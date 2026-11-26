# pyfortransitions

Python interface to **fortransitions**, the Fortran port of
[CosmoTransitions](https://github.com/clwainwright/CosmoTransitions). It
provides `fullTunneling` and `SingleFieldInstanton` with the same
arguments and defaults as cosmoTransitions. The potential is written in
Python; the instanton calculation runs in Fortran.

The Fortran code is compiled into a shared library
(`pyfortransitions/lib/libfortransitions.so`) and called with `ctypes`
through `bind(C)` wrappers (`../src/c_*.f90`). Only Linux is supported
for now.

## Installation

### With pixi (recommended)

[pixi](https://pixi.sh) provides the whole toolchain (gfortran, fpm,
Python, numpy) from conda-forge. Nothing has to be installed on the
system. From this directory:

```sh
pixi run install   # build the shared library, install the package
pixi run test      # build, install and run the tests
pixi shell         # a shell in the environment, with pyfortransitions
```

### Manual build

This requires gfortran (or ifx), [fpm](https://fpm.fortran-lang.org),
make, and a Python environment with pip and numpy. From this directory,
with that Python environment active:

```sh
make wrapper           # gfortran
make wrapper FC=ifx    # Intel compiler
python test/run_tests.py
```

`make wrapper` builds the library with fpm (release profile) and installs
the package in editable mode (`pip install -e .`). The libraries find
each other at run time through their runpath, so no `patchelf` step is
needed. `make clean` uninstalls the package and removes all build
products. Use `make wrapper PC=python3` if your interpreter is not called
`python`.

## Usage

### Multi-field instantons

```python
import numpy as np
from pyfortransitions import fullTunneling

def V(X):
    x, y = X[..., 0], X[..., 1]
    r1 = x*x + 5*y*y
    r2 = 5*(x - 1)**2 + (y - 1)**2
    return r1*r2 + 80*(0.25*y**4 - y**3/3)

def dV(X):
    x, y = X[..., 0], X[..., 1]
    r1 = x*x + 5*y*y
    r2 = 5*(x - 1)**2 + (y - 1)**2
    return np.array([r1*10*(x - 1) + 2*x*r2,
                     r1*2*(y - 1) + 10*y*r2 + 80*y*y*(y - 1)])

# First point at the stable minimum, last point at the metastable one.
res = fullTunneling([[1.0, 1.0], [0.0, 0.0]], V, dV)
print(res.action, res.fRatio)
```

`V` and `dV` are called with a 1d array `X` of shape `(ndim,)`. Potentials
written for cosmoTransitions (indexing with `X[..., i]`) therefore work
unchanged. The result has the fields of cosmoTransitions'
`fullTunneling_rval` (`profile1D`, `Phi`, `action`, `fRatio`,
`saved_steps`). It also has the attributes `action_pot` and `action_kin`
(the action from only the potential or only the kinetic term, rescaled by
Derrick's theorem; their spread from `action` measures the accuracy),
`num_iters` and `converged`.

Use `V_spline_samples=None` to evaluate the potential directly along the
path instead of through a spline of 100 samples. This is needed if the
barrier is narrow compared to the path length, e.g. very close to the
metastable minimum.

### One-field instantons

```python
from pyfortransitions import SingleFieldInstanton

V = lambda phi: 0.25*phi**4 - 0.49*phi**3 + 0.235*phi**2
dV = lambda phi: phi*(phi - 0.47)*(phi - 1)

inst = SingleFieldInstanton(1.0, 0.0, V, dV)  # dV, d2V are optional
profile = inst.findProfile()
print(inst.findAction(profile))
```

### Errors

Failures raise exceptions with the cosmoTransitions names:
`PotentialError` (with `args[1]` either `"no barrier"` or
`"stable, not metastable"`), `IntegrationError` and `DeformationError`.
All of them derive from `FortransitionsError`, which is also raised
directly for other numerical failures. An exception raised inside `V` or
`dV` aborts the calculation and is re-raised unchanged (this includes
`KeyboardInterrupt`).

## Differences from cosmoTransitions

- `fullTunneling` supports these parameter keys: `alpha` in
  `tunneling_init_params`; `npoints`, `xtol`, `phitol` in
  `tunneling_findProfile_params`; `nb`, `kb`, `fix_start`, `fix_end` in
  `deformation_init_params`; and `verbose` in
  `deformation_deform_params`. The options `save_all_steps` and
  `callback`, other tunneling or deformation classes, and any other keys
  raise `NotImplementedError`. `fixEndCutoff` is ignored, as it is in
  cosmoTransitions.
- The deformation messages are printed only if `verbose` (or
  `deformation_deform_params["verbose"]`) is true. cosmoTransitions always
  prints them.
- The path fit has no limit on the number of fields. cosmoTransitions
  stops at 10, a limit of `scipy.interpolate.splprep`.
- `SingleFieldInstanton.findActionDerrick(profile)` returns the action
  together with its two Derrick estimates.
- Only the methods listed above are available on `SingleFieldInstanton`,
  plus the attributes `phi_bar`, `rscale`, `phi_eps`, `alpha` and `V`.

## Tests

`test/run_tests.py` repeats the Fortran tests (`../test/check.f90`,
`double_well.f90` and `many_fields.f90`) with the potentials written in
Python. It asserts the same cosmoTransitions reference values with the
same tolerances. It also tests error handling, argument validation and
memory use. `python test/run_tests.py --quick` skips the many-field scan
above 10 fields.

`test/compare_cosmoTransitions.py` runs the same problems with
cosmoTransitions and pyfortransitions, using the same potential functions
and keyword arguments. It prints the results, their relative differences
and the run times of both packages. Run it in the pixi environment
`compare`, which adds cosmoTransitions and scipy and pins Python, numpy
and scipy to the versions used for the reference values:

```sh
pixi run -e compare compare            # or: ... compare --quick
```

## Performance

The Fortran code calls `V` and `dV` once per point, while cosmoTransitions
evaluates them on whole arrays of points. For one-field problems and
moderate path deformations, pyfortransitions is several times to ~30x
faster. When the path deformation needs many steps, the per-point Python
calls dominate the run time, and pyfortransitions can be slower than
cosmoTransitions: about 3x slower for the two-field example of
cosmoTransitions. Most of the time then goes into evaluating `dV` in
numpy for a single point. Writing `V` and `dV` with plain floats for one
point, instead of numpy array operations, reduces this. For the
two-field example, this brings the time from 2.6 s to 1.0 s, with an
identical result (cosmoTransitions: 0.7 s):

```python
def dV(X):
    x, y = float(X[0]), float(X[1])
    ...
    return [dVdx, dVdy]
```
