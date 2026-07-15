"""Reference values for test/many_fields.f90.

Runs the Python cosmoTransitions package on a family of potentials with
an increasing number of field dimensions ndim,

    V(x) = 0.4*phi^2 - 1.6*phi^3 + phi^4
           + 1/2 * sum_{i=2}^{ndim} m2_i * (x_i - a_i*phi^2)^2,

with phi = x_1, a_i = 0.04*(i-1) and m2_i = 2 + 0.2*(i-1). The minima
are known exactly by construction: the metastable minimum is the origin
(V = 0) and the true minimum is x = (1, a_2, ..., a_ndim) (V = -0.2),
with the barrier top along phi at phi = 0.2. The transverse fields want
to sit on the curve x_i = a_i*phi^2, so the tunneling path bends in
every transverse direction and the path deformation is genuinely
exercised. The initial guess is the straight line between the minima.

The scan stops at ndim = 10: for ndim > 10 the Python package fails
with "TypeError: 0 < idim < 11 must hold", raised by
scipy.interpolate.splprep (a hard limit of the underlying FITPACK
parcur routine, which SplinePath uses to fit the path spline). The
Fortran port fits each field dimension as its own 1d B-spline and has
no such limit; test/many_fields.f90 therefore scans ndim = 2..20,
asserting the actions printed here for ndim = 2..10 and internal
consistency checks beyond.

The findProfile tolerances are tightened to xtol = phitol = 1e-9 (the
defaults of 1e-4 leave a relative error of up to a few 1e-3 in the
action); the Fortran test uses the same tolerances and asserts at a
relative tolerance of 2e-3 (see reference/README.md).

Usage: python3 ref_many_fields.py

Requires the cosmoTransitions package to be importable (e.g. installed
with pip).
"""

import numpy as np
from cosmoTransitions import pathDeformation as pd


def make_pot(ndim):
    a = 0.04 * np.arange(1, ndim)         # a_i, i = 2..ndim
    m2 = 2.0 + 0.2 * np.arange(1, ndim)   # m2_i, i = 2..ndim

    def V(X):
        X = np.asanyarray(X)
        phi = X[..., 0]
        d = X[..., 1:] - a * phi[..., None]**2
        return (0.4 * phi**2 - 1.6 * phi**3 + phi**4
                + 0.5 * np.sum(m2 * d**2, axis=-1))

    def dV(X):
        X = np.asanyarray(X)
        phi = X[..., 0]
        d = X[..., 1:] - a * phi[..., None]**2
        rval = np.empty_like(X)
        rval[..., 0] = (0.8 * phi - 4.8 * phi**2 + 4.0 * phi**3
                        - 2.0 * phi * np.sum(m2 * a * d, axis=-1))
        rval[..., 1:] = m2 * d
        return rval

    x_true = np.concatenate(([1.0], a))
    return V, dV, x_true


npath = 41

for ndim in range(2, 11):
    V, dV, x_true = make_pot(ndim)

    # exact minima by construction
    assert abs(V(np.zeros(ndim))) < 1e-15
    assert abs(V(x_true) + 0.2) < 1e-15
    assert np.max(np.abs(dV(np.zeros(ndim)))) < 1e-15
    assert np.max(np.abs(dV(x_true))) < 1e-15

    # straight line from the true to the metastable minimum
    t = np.linspace(0.0, 1.0, npath)[:, None]
    path_pts = x_true * (1.0 - t)

    Y = pd.fullTunneling(path_pts, V, dV, verbose=False,
                         tunneling_findProfile_params={"xtol": 1e-9,
                                                       "phitol": 1e-9})
    print("ndim = %2d  S_3 = %.10e  fRatio = %.6e"
          % (ndim, Y.action, Y.fRatio))
