# fortransitions with LFortran

fortransitions builds and passes all tests with LFortran
(version 0.68.0-243-gf97c7d364b), in the debug and the release profile:

```sh
fpm test --profile release --compiler lfortran
fpm test --profile debug --compiler lfortran
```

The results agree with gfortran to better than 1e-6 for most values; the
largest difference (1e-3, double well at c = 0.25) is the known spread of
the path deformation, also seen between gfortran, ifx and cosmoTransitions.
The example program gives the same actions as with gfortran.

This file lists the problems found on the way, with their causes and
what was done about them.

## Summary

| # | Problem | Cause | Solution |
| - | ------- | ----- | -------- |
| 1 | Link error `_lfortran_string_read_f128` | LFortran runtime lacks the function | csv-fortran dependency removed |
| 2 | "Array is not allocated" at run time | No reallocation on assignment by default | `--realloc-lhs-arrays` in the profiles |
| 3 | Profile integration fails, all actions 0 | **LFortran bug** (see below) | workaround in `sfi_integrate_profile` |
| 4 | `--fast` breaks the NaN checks | fast math | release profile uses `--fast --no-fast-math` |
| 5 | Minimizer status in `test/double_well.f90` | test tolerance below rounding noise | test accepts status 3 |
| 6 | Debug builds fail without `llvm-dwarfdump` | `-g` needs `llvm-dwarfdump` on the PATH | environment (LFortran's build environment has it) |

## 1. Link error `_lfortran_string_read_f128`

Linking the example program failed with an undefined reference to
`_lfortran_string_read_f128`, called from csv-fortran's `to_real_qp`
(an internal `read(str, *) val` with `real(real128) :: val`). LFortran
emits a call to this function, but its runtime only defines
`_lfortran_string_read_f32` and `_lfortran_string_read_f64` (and the
`_array` variants); the function is missing from the runtime library and
from `src/libasr/runtime/lfortran_intrinsics.c`.

The example program was the only user of csv-fortran. It now writes its
CSV files with plain `write` statements, and csv-fortran is no longer a
dependency. Worth reporting upstream as a missing runtime function.

## 2. Arrays not allocated at run time

The code relies on the Fortran 2003 rule that assignment to an
unallocated (or differently shaped) allocatable array allocates it, as do
odeint and the test programs. LFortran does this only with
`--realloc-lhs-arrays`, which is in both LFortran profiles of `fpm.toml`.
fpm passes the profile flags to the dependencies as well.

## 3. LFortran bug: value lost through an argument (worked around)

**Symptom:** `single_field_instanton%find_profile` returned status 1
(integration error), so `full_tunneling` failed for every case. odeint's
`rkqs_vec` returned the step size `hdid = 0` to `sfi_integrate_profile`,
although it assigns `hdid = h` (all other outputs were correct), and the
check `dr_did < drmin` stopped the integration.

**Cause:** when a procedure is called with an internal procedure as an
actual argument, a value that the callee assigns to an output argument is
lost, if that variable is also used by another internal procedure of the
same host. In `sfi_integrate_profile`, `dr_did` is passed to `rkqs_vec`
together with the internal procedure `deriv`, and is used by the internal
procedures `interp_dphi` and `interp_phi_diff`. It looks as if LFortran
copies the variables used by internal procedures into separate storage
for the call and writes them back afterwards, overwriting the output.

**Minimal reproducer** (expected output: `42.0` everywhere; LFortran
prints `1.0` in case A):

```fortran
module m
  implicit none
  abstract interface
    function fun(x) result(y)
      real, intent(in) :: x
      real :: y
    end function fun
  end interface
contains
  subroutine set_it(f, x)
    procedure(fun) :: f
    real, intent(out) :: x
    x = 42.0 + 0.0*f(1.0)
  end subroutine set_it

  function plain(x) result(y)
    real, intent(in) :: x
    real :: y
    y = x
  end function plain

  subroutine case_a()
    ! internal procedure passed; another internal procedure uses `a`
    real :: a
    a = 1.0
    call set_it(g, a)
    print *, "A internal g passed, h uses a :", a, h()
  contains
    function g(x) result(y)
      real, intent(in) :: x
      real :: y
      y = x
    end function g
    function h() result(y)
      real :: y
      y = a
    end function h
  end subroutine case_a

  subroutine case_b()
    ! module procedure passed; an internal procedure uses `a`
    real :: a
    a = 1.0
    call set_it(plain, a)
    print *, "B module proc passed, h uses a:", a, h()
  contains
    function h() result(y)
      real :: y
      y = a
    end function h
  end subroutine case_b

  subroutine case_c()
    ! internal procedure passed; no internal procedure uses `a`
    real :: a
    a = 1.0
    call set_it(g, a)
    print *, "C internal g passed, a unused :", a
  contains
    function g(x) result(y)
      real, intent(in) :: x
      real :: y
      y = x
    end function g
  end subroutine case_c
end module m

program min
  use m
  implicit none
  call case_a()
  call case_b()
  call case_c()
end program min
```

```
$ gfortran min.f90 && ./a.out
 A internal g passed, h uses a :   42.0000000       42.0000000
 B module proc passed, h uses a:   42.0000000       42.0000000
 C internal g passed, a unused :   42.0000000
$ lfortran min.f90 && ./a.out
A internal g passed, h uses a : 1.00000000 1.00000000
B module proc passed, h uses a: 42.0000000 42.0000000
C internal g passed, a unused : 42.0000000
```

**Workaround:** in `sfi_integrate_profile`, `rkqs_vec` returns the step
size in a separate local variable (`h_did`), which no internal procedure
uses, and it is copied to `dr_did` afterwards. The other calls in `src/`
that pass internal procedures (to `brentq`, `minimize_bounded`,
`minimize_unbounded`, and `rkqs_vec` in `sfi_integrate_and_save_profile`)
do not assign to variables used by internal procedures, so they are not
affected. Code changes in this area must keep this pattern in mind.

## 4. `--fast` and fast math

`--fast` makes LFortran builds about 7-8 times faster, similar to
gfortran. But it enables fast math, under which `ieee_is_finite` reports
NaN as finite and comparisons with NaN come out true. The guards that
keep the code from hanging on non-finite values (e.g. from a potential,
or after an exception in a Python callback) rely on these checks. The
release profile therefore uses `--fast --no-fast-math`, which keeps the
NaN semantics at the same speed, much like ifx needs `-fp-model=precise`.

(Before the workaround for problem 3, `--fast` also made `check` and
`c_interface` crash; that was a consequence of problem 3, since without
bounds checking the code then accessed the never allocated profile.)

## 5. Minimizer status in `test/double_well.f90`

With `--fast --no-fast-math`, gradmin's `minimize` returned status 3
(maximum number of iterations) for the true minimum at c = 0.25. The test
asks for `|grad V| < 1e-12`, below the rounding noise of gradmin's
finite-difference gradient (step 1e-8, noise ~1e-9), so reaching it
depends on the last bits of the arithmetic. The minimum itself is
correct. The test now accepts status 3; the minimum is still checked
against the reference values, and gfortran results are unchanged.

## 6. Debug builds need `llvm-dwarfdump`

LFortran's `-g` (part of the general debug flags in `fpm.toml`) calls
`llvm-dwarfdump`. Without it on the PATH, every compilation in the debug
profile fails ("Error in creating the files used to generate the debug
information"). LFortran's own build environment provides it.
