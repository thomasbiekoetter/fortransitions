# Reference scripts

These Python scripts call the original Python package
[cosmoTransitions](https://github.com/clwainwright/CosmoTransitions) to
produce the reference values (bounce actions, minima locations, barrier
positions, radial scales) that the Fortran test programs assert against:

| Script               | Produces reference values for |
| -------------------- | ----------------------------- |
| `ref_check.py`       | `test/check.f90`              |
| `ref_double_well.py` | `test/double_well.f90`        |

Run them with

```sh
python3 ref_check.py
python3 ref_double_well.py
```

The scripts require the `cosmoTransitions` package to be importable
(e.g. installed with pip), along with numpy and scipy.

The values hard-coded in the Fortran tests were generated with
Python 3.13, numpy 2.2.6 and scipy 1.16.1. The Fortran port and the
Python package agree to about 1e-4 (relative) or better on the actions;
the residual differences come from different spline knot placement
(bspline-fortran vs. FITPACK) and slightly different quadrature details,
which is why the tests use a relative tolerance of 5e-3.
