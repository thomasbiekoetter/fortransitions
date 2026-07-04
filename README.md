# cosmotransitions

A Fortran port of the Python package
[CosmoTransitions](https://github.com/clwainwright/CosmoTransitions)
(C.L. Wainwright, Comput. Phys. Commun. 183 (2012) 2006,
[arXiv:1109.4189](http://arxiv.org/abs/1109.4189)).

This port was written entirely by
[Claude Fable 5](https://www.anthropic.com/news/claude-fable-5-mythos-5)
(Anthropic), including the source code, the test programs, the reference
scripts in `reference/`, and this README.

## Ported so far

- `pathDeformation.fullTunneling` → `full_tunneling`: multi-field
  instanton calculation via iterative path deformation.
- `tunneling1D.SingleFieldInstanton` → `single_field_instanton`:
  one-dimensional instantons via the overshoot/undershoot method.
- `pathDeformation.SplinePath` → `spline_path` (spline mode only, i.e.
  `V_spline_samples > 0`).
- `pathDeformation.Deformation_Spline` → `deformation_spline`
  (without `save_all_steps`).
- The required pieces of `helper_functions` (rkqs/rkck, Nbspld2, finite
  differences, ...) plus replacements for the scipy routines that the
  Python package relies on (brentq, fminbound/fmin, splprep/splrep via
  bspline-fortran, lstsq via normal equations, fractional-order Bessel
  functions).

Errors that the Python package raises as exceptions
(`PotentialError`, `IntegrationError`, `DeformationError`) are reported
through integer `status` arguments; the codes are defined in
`cosmotransitions__config` and re-exported by the top-level module.

## Usage

Potentials are supplied as derived types. For `full_tunneling`, extend
`potential_nd` and implement `v(x)` and `grad(x)`:

```fortran
use cosmotransitions

type, extends(potential_nd) :: my_potential
contains
  procedure :: v => my_v
  procedure :: grad => my_grad
end type
```

Then, with the first path point at (or near) the stable minimum and the
last one at the metastable minimum:

```fortran
type(my_potential), target :: pot
type(full_tunneling_result) :: res
real(wp) :: path_pts(2, 2)
integer :: status

path_pts(1, :) = [1.0_wp, 1.0_wp]  ! stable (deeper) minimum
path_pts(2, :) = [0.0_wp, 0.0_wp]  ! metastable minimum
call full_tunneling(path_pts, pot, res, status)
print *, res%action, res%fratio
```

`res%profile` holds the one-dimensional bubble profile (r, phi, dphi)
and `res%phi` the corresponding points of the deformed path in field
space. As additional diagnostics, `res%action_pot` and `res%action_kin`
give the action computed from only the potential or only the kinetic
term of the integrand, rescaled by Derrick's theorem; for an exact
solution of the bounce equation both equal `res%action`, so their
deviation measures the quality of the solution (like `res%fratio`). See `app/full_tunneling_example.f90` for a complete example (the two-field example
of the Python package, `examples/fullTunneling.py`) including CSV
output, and `test/check.f90` for validation against reference values
computed with the Python package.

For purely one-dimensional problems, use `single_field_instanton` with a
`potential_1d` (or the ready-made `pot1d_func` wrapper for plain
functions; derivatives default to finite differences).

## Building and testing

```sh
fpm build --profile release
fpm test --profile release   # validates against Python reference values
fpm run --profile release    # runs the two-field example, writes CSVs
```

Always pass `--profile debug` or `--profile release` explicitly: without
it, fpm compiles with no flags at all, so none of the flags defined in
`fpm.toml` are applied. In particular, ifx then falls back to its
value-unsafe floating-point default (`-fp-model=fast`), which breaks the
thin-wall shooting; both profiles therefore set `-fp-model=precise` for
ifx.

The reference values asserted in the test programs were produced with the
original Python package; the scripts that generate them are collected in
`reference/` (see `reference/README.md`) for reproducibility.

## License

MIT, the same license as the original CosmoTransitions Python package;
see [LICENSE](LICENSE). As a derivative work, the license retains the
original package's copyright notice.
