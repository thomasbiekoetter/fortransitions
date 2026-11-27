# fortransitions for Python

The Python package `fortransitions` is the Python interface to the
Fortran package of the same name, the Fortran port of
[CosmoTransitions](https://github.com/clwainwright/CosmoTransitions). It
provides `fullTunneling` and `SingleFieldInstanton` with the same
arguments and defaults as cosmoTransitions. The potential is written in
Python; the instanton calculation runs in Fortran.

The Fortran code is compiled into a shared library
(`fortransitions/lib/libfortransitions.so`) and called with `ctypes`
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
pixi shell         # a shell in the environment, with fortransitions
```

With the Intel compiler ifx (pixi environment `ifx`):

```sh
pixi run -e ifx install
pixi run -e ifx test
pixi run -e compare-ifx compare    # comparison with cosmoTransitions
```

The ifx build links Intel's runtime libraries from that environment, so
it only works while the environment exists.

With LFortran (pixi environment `lfortran`):

```sh
pixi run -e lfortran install
pixi run -e lfortran test
pixi run -e compare-lfortran compare    # comparison with cosmoTransitions
```

### Manual build

This requires gfortran (or ifx), [fpm](https://fpm.fortran-lang.org),
make, and a Python environment with pip and numpy. From this directory,
with that Python environment active:

```sh
make wrapper           # gfortran
make wrapper FC=ifx    # Intel compiler
make wrapper FC=lfortran
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
from fortransitions import fullTunneling

def V(X):
    x, y = X[..., 0], X[..., 1]
    r1 = x*x + 5*y*y
    r2 = 5*(x - 1)**2 + (y - 1)**2
    return r1*r2 + 80*(0.25*y**4 - y**3/3)

def dV(X):
    x, y = X[..., 0], X[..., 1]
    r1 = x*x + 5*y*y
    r2 = 5*(x - 1)**2 + (y - 1)**2
    g = np.empty_like(X)
    g[..., 0] = r1*10*(x - 1) + 2*x*r2
    g[..., 1] = r1*2*(y - 1) + 10*y*r2 + 80*y*y*(y - 1)
    return g

# First point at the stable minimum, last point at the metastable one.
res = fullTunneling([[1.0, 1.0], [0.0, 0.0]], V, dV)
print(res.action, res.fRatio)
```

`V` and `dV` are called with a single point, a 1d array `X` of shape
`(ndim,)`. In the path deformation, `dV` is also called with all points
of the path at once, an array of shape `(npts, ndim)`, and must then
return an array of the same shape. This is the same convention as in
cosmoTransitions, so potentials written for it (indexing with `X[..., i]`)
work unchanged. A `dV` that only handles single points is detected and
evaluated point by point instead, which is slower (see
[Performance](#performance)). The option `vectorized=True` or `False`
skips the detection. The result has the fields of cosmoTransitions'
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
from fortransitions import SingleFieldInstanton

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
cosmoTransitions and fortransitions, using the same potential functions
and keyword arguments. It prints the results, their relative differences
and the run times of both packages. Run it in the pixi environment
`compare`, which adds cosmoTransitions and scipy and pins Python, numpy
and scipy to the versions used for the reference values:

```sh
pixi run -e compare compare            # or: ... compare --quick
```

## Performance

The 1d tunneling runs entirely in Fortran and calls `V` and `dV` one
point at a time. The path deformation evaluates `dV` at all points of the
path, about 750 of them, at every deformation step. It does this in a
single call per step, as cosmoTransitions does, so numpy evaluates all
points in one vectorized operation.

For this to work, `dV` must accept an array of shape `(npts, ndim)`, as
in the example above. Otherwise fortransitions falls back to calling
`dV` once per point, which costs about 50 µs per point. For the
two-field example of cosmoTransitions (`examples/fullTunneling.py`):

| `dV` | fortransitions | cosmoTransitions |
| ---- | ---------------- | ---------------- |
| vectorized (`X[..., i]`, `np.empty_like`) | 0.05 s | 0.66 s |
| single points only, plain floats | 0.30 s | (fails) |
| vectorized, but `vectorized=False` | 2.6 s | |

For all test problems of `test/compare_cosmoTransitions.py`, with
vectorized potentials, fortransitions is 10–40x faster for one-field
problems, about 10x for the two-field example, about 60x for the double
well and 40–70x for 2–10 fields. The narrow barrier with
`V_spline_samples=None` is 5x faster: there, the 1d tunneling has to call
`V` and `dV` one point at a time. Over all problems it is 27x faster
(1.9 s vs. 52 s).
