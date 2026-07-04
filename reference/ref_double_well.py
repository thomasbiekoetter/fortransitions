"""Reference values for test/double_well.f90.

Runs the Python cosmoTransitions package on the two-field potential

    V(x, y) = (x^2 + y^2) * ( a*(x-1)^2 + b*(y-1)^2 - c )

with a = 1.8, b = 0.2, scanning the coefficient c between 0.01 and 0.3.
For each c the true minimum near (1, 1) is located with scipy (BFGS with
analytic gradient), and the O(3)-symmetric (finite-temperature) bounce
action S_3 is computed with fullTunneling. The findProfile tolerances
are tightened to xtol = phitol = 1e-9 (the defaults of 1e-4 leave a
relative error of up to a few 1e-3 in the action); the Fortran test
uses the same tolerances and asserts at a relative tolerance of 2e-3.
The residual differences of up to ~1e-3 come from slightly different
path-deformation trajectories, which depend on the implementation, the
compiler and the optimization level.
The printed minima and actions are the reference values asserted in
test/double_well.f90.

Usage: python3 ref_double_well.py

Requires the cosmoTransitions package to be importable (e.g. installed
with pip).
"""

import numpy as np
from scipy import optimize

from cosmoTransitions import pathDeformation as pd

a, b = 1.8, 0.2


def make_pot(c):
    def V(X):
        x, y = X[..., 0], X[..., 1]
        return (x**2 + y**2) * (a*(x - 1.0)**2 + b*(y - 1.0)**2 - c)

    def dV(X):
        x, y = X[..., 0], X[..., 1]
        f = a*(x - 1.0)**2 + b*(y - 1.0)**2 - c
        r = x**2 + y**2
        rval = np.empty_like(X)
        rval[..., 0] = 2.0*x*f + r*2.0*a*(x - 1.0)
        rval[..., 1] = 2.0*y*f + r*2.0*b*(y - 1.0)
        return rval

    return V, dV


for c in [0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3]:
    V, dV = make_pot(c)
    res = optimize.minimize(lambda x: V(x), [1.0, 1.0], jac=lambda x: dV(x),
                            method="BFGS", tol=1e-14)
    xmin = res.x
    Y = pd.fullTunneling([xmin, [0.0, 0.0]], V, dV, verbose=False,
                         tunneling_findProfile_params={"xtol": 1e-9,
                                                       "phitol": 1e-9})
    print("c = %.2f  min = (%.12f, %.12f)  Vmin = %.6e  S3 = %.10e"
          % (c, xmin[0], xmin[1], V(xmin), Y.action))
