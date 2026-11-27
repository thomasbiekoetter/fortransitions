# fortransitions

A Fortran port of the Python package
[CosmoTransitions](https://github.com/clwainwright/CosmoTransitions)
(C.L. Wainwright, Comput. Phys. Commun. 183 (2012) 2006,
[arXiv:1109.4189](http://arxiv.org/abs/1109.4189)).

This port was written with AI models by Anthropic. The initial port was
written by
[Claude Fable 5](https://www.anthropic.com/news/claude-fable-5-mythos-5),
including the source code, the test programs, the reference scripts in
`reference/`, and this README. Later additions, among them the Python
interface in `python/`, were written with Claude Opus 5.5.

## Citation

If you use this code, please cite the original CosmoTransitions paper
(C.L. Wainwright, Comput. Phys. Commun. 183 (2012) 2006,
[arXiv:1109.4189](http://arxiv.org/abs/1109.4189)) and, in addition,
the paper for which this Fortran port was developed
([INSPIRE](https://inspirehep.net/literature/3231213)):

```bibtex
@article{Biekotter:2026oiv,
    author = {Biek{\"o}tter, Thomas and Dashko, Andrii and Radchenko, Kateryna},
    title = "{Simulation-based inference for gravitational-wave cosmology: parameter reconstruction of a first-order electroweak phase transition}",
    eprint = "2611.25036",
    archivePrefix = "arXiv",
    primaryClass = "astro-ph.CO",
    reportNumber = "IFT-UAM/CSIC-26-169",
    month = "11",
    year = "2026"
}
```

## Ported so far

- `pathDeformation.fullTunneling` → `full_tunneling`: multi-field
  instanton calculation via iterative path deformation.
- `tunneling1D.SingleFieldInstanton` → `single_field_instanton`:
  one-dimensional instantons via the overshoot/undershoot method.
- `pathDeformation.SplinePath` → `spline_path`, with two modes selected
  by the optional argument `v_mode` of `full_tunneling` (and
  `spline_path%init`): `"spline"` (default, like Python) interpolates
  `v_spline_samples` samples of V along the path (default 100), and
  `"direct"` (`V_spline_samples=None` in Python) evaluates V and its
  gradient directly along the path; d²V/dx² is taken from a central
  difference of the exact dV/dx (Python: finite differences of V). Use
  `v_mode = "direct"` if the barrier is narrow compared to the path
  length, e.g. very close to the metastable minimum: 100 samples of V
  then do not resolve it, which can make the one-dimensional tunneling
  fail ("no barrier") or give inconsistent profiles (`test/check.f90`
  has such a case). Unlike the Python class, which inherits a hard limit
  of 10 field dimensions from `scipy.interpolate.splprep`
  (FITPACK's `parcur`), the port fits each field dimension as its own
  1d B-spline and works for any number of fields
  (`test/many_fields.f90` exercises up to 20).
- `pathDeformation.Deformation_Spline` → `deformation_spline`
  (without `save_all_steps`).
- The required pieces of `helper_functions` (Nbspld2, finite
  differences, ...) plus replacements for the scipy routines that the
  Python package relies on (brentq, fminbound/fmin, splprep/splrep via
  bspline-fortran, lstsq via normal equations, fractional-order Bessel
  functions). The adaptive Cash-Karp Runge-Kutta stepper
  (`helper_functions.rkqs`/`_rkck` in the Python package) is provided by
  the [odeint](https://github.com/thomasbiekoetter/odeint) package
  (`rkqs_vec` in `odeint__rkck`, which supports the per-component
  tolerances and status-code error reporting that the shooting method
  needs).

Errors that the Python package raises as exceptions
(`PotentialError`, `IntegrationError`, `DeformationError`) are reported
through integer `status` arguments; the codes are defined in
`fortransitions__config` and re-exported by the top-level module.

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

The package builds and passes all tests with gfortran, ifx and LFortran:

```sh
fpm test --profile release --compiler ifx
fpm test --profile release --compiler lfortran
```

For LFortran, the profiles add the required options, among them
`--realloc-lhs-arrays`. The release profile uses `--fast --no-fast-math`,
since fast math would break the checks for non-finite values. See
`LFORTRAN_ISSUES.md` for the LFortran problems found and their
workarounds.

The reference values asserted in the test programs were produced with the
original Python package; the scripts that generate them are collected in
`reference/` (see `reference/README.md`) for reproducibility.

## Using fortransitions in your fpm project

Add fortransitions to the `[dependencies]` of your `fpm.toml`:

```toml
[dependencies]
fortransitions = { git = "https://github.com/thomasbiekoetter/fortransitions.git", profile = "release" }
```

`profile = "release"` compiles fortransitions with its own release flags,
which for ifx include the required `-fp-model=precise` (see
[Building and testing](#building-and-testing)).

## Usage

Potentials are supplied as derived types. For `full_tunneling`, extend
`potential_nd` and implement `v(x)` and `grad(x)`:

```fortran
use fortransitions

type, extends(potential_nd) :: my_potential
contains
  procedure :: v => my_v
  procedure :: grad => my_grad
end type
```

Optionally, also override `grad_many(x, dv)`. It returns the gradient at
all rows of an `(npts, ndim)` array at once, and the path deformation uses
it for the forces on the whole path. By default it calls `grad` for each
point.

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

## Python interface

The directory `python/` contains the Python package `fortransitions`.
It exposes `fullTunneling` and `SingleFieldInstanton` with the same
arguments as cosmoTransitions, with the potential written in Python.

### Installation

The Python package is built and installed from the `python/` folder.
Both installations use the same compiled library in
`python/fortransitions/lib/`, from the last `pixi run install` or
`make wrapper`.

#### Recommended: with pixi

[pixi](https://pixi.sh) provides gfortran, fpm, make, Python and numpy.
The package is installed into the pixi environment (`python/.pixi/`).

Build the library and install the package:

```sh
cd python
pixi run install
```

Run the tests, and the comparison with cosmoTransitions:

```sh
pixi run test
pixi run -e compare compare
```

Build and test with the Intel compiler ifx instead of gfortran:

```sh
pixi run -e ifx install
pixi run -e ifx test
pixi run -e compare-ifx compare
```

The ifx build links Intel's runtime libraries from the pixi environment
`ifx`, so it only works while that environment exists.

The same with LFortran (from conda-forge):

```sh
pixi run -e lfortran install
pixi run -e lfortran test
pixi run -e compare-lfortran compare
```

Uninstall and remove the build products:

```sh
pixi run clean
```

#### Manual: with make

This installs the package into the active Python environment. It
requires gfortran or ifx, fpm, make, pip and numpy.

Build the library and install the package (with gfortran):

```sh
cd python
make wrapper
```

With the Intel compiler or LFortran, or with a Python interpreter that is
not called `python`:

```sh
make wrapper FC=ifx
make wrapper FC=lfortran
make wrapper PC=python3
```

Run the tests:

```sh
python test/run_tests.py
```

Uninstall and remove the build products:

```sh
make clean
```

### Usage

After the installation, for example:

```python
from fortransitions import fullTunneling
res = fullTunneling([[1.0, 1.0], [0.0, 0.0]], V, dV)
print(res.action)
```

The package calls the Fortran library through `bind(C)` wrappers in
`src/c_*.f90`. See `python/README.md` for the full documentation of the
Python interface.

## License

MIT, the same license as the original CosmoTransitions Python package;
see [LICENSE](LICENSE). As a derivative work, the license retains the
original package's copyright notice.
