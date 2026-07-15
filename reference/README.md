# Reference scripts

These Python scripts call the original Python package
[cosmoTransitions](https://github.com/clwainwright/CosmoTransitions) to
produce the reference values (bounce actions, minima locations, barrier
positions, radial scales) that the Fortran test programs assert against:

| Script               | Produces reference values for |
| -------------------- | ----------------------------- |
| `ref_check.py`       | `test/check.f90`              |
| `ref_double_well.py` | `test/double_well.f90`        |
| `ref_many_fields.py` | `test/many_fields.f90`        |

Run them with

```sh
python3 ref_check.py
python3 ref_double_well.py
python3 ref_many_fields.py
```

The scripts require the `cosmoTransitions` package to be importable
(e.g. installed with pip), along with numpy and scipy.

The values hard-coded in the Fortran tests were generated with
Python 3.13, numpy 2.2.6 and scipy 1.16.1. With the default settings
both the Fortran port and the Python package compute the 1d bounce
profile with the findProfile tolerances xtol = phitol = 1e-4, which
leaves a relative error of up to a few 1e-3 in the action (visible in
the Fortran result as the spread between `action`, `action_pot` and
`action_kin`); this is why `test/check.f90`, which runs with the
default tolerances, asserts at a relative tolerance of 5e-3.
`test/double_well.f90` and `ref_double_well.py` instead tighten the
tolerances to xtol = phitol = 1e-9 (exposed as optional arguments of
`full_tunneling` in the Fortran port and via
`tunneling_findProfile_params` in Python); the actions of the two
codes then agree to about 1e-7 on most points, with occasional
differences up to ~1e-3 from slightly different path-deformation
trajectories (which depend on the implementation, the compiler and the
optimization level), and the test asserts at a relative tolerance of
2e-3.
