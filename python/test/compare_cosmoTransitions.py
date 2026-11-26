"""Side-by-side comparison of fortransitions with cosmoTransitions.

Runs the same problems with both packages, using the very same Python
potential functions and the same keyword arguments, and prints the
results, their relative differences and the run times. The problems are
those of the Fortran test programs (test/check.f90, test/double_well.f90,
test/many_fields.f90); the tolerances of the comparison are the ones used
there. With the default findProfile tolerances (1e-4) the actions of the
two codes agree to a few 1e-3; with 1e-9 they agree to ~1e-7 on most
points, with occasional differences up to ~1e-3 from slightly different
path-deformation trajectories (see reference/README.md).

Run with ``pixi run -e compare compare`` (from python/), or with
``python test/compare_cosmoTransitions.py`` in an environment where both
packages are installed. ``--quick`` skips the slow cases. The exit code
is 1 if any difference exceeds its tolerance.
"""

import contextlib
import io
import sys
import time

import numpy as np
from scipy import optimize

import cosmoTransitions
from cosmoTransitions import pathDeformation as ct_pd
from cosmoTransitions import tunneling1D as ct_t1d

import fortransitions as ft

# ---------------------------------------------------------------------------
# Potentials, written in the vectorized cosmoTransitions style (X[..., i]),
# as cosmoTransitions evaluates them on arrays of points. fortransitions
# calls the same functions with one point of shape (ndim,).
# ---------------------------------------------------------------------------


def v1_thin(phi):
    return 0.25 * phi**4 - 0.49 * phi**3 + 0.235 * phi**2


def dv1_thin(phi):
    return phi * (phi - 0.47) * (phi - 1.0)


def v1_thick(phi):
    return 0.25 * phi**4 - 0.4 * phi**3 + 0.1 * phi**2


def dv1_thick(phi):
    return phi * (phi - 0.2) * (phi - 1.0)


def example_2d(c=5.0, fx=0.0, fy=2.0):
    """The two-field potential of examples/fullTunneling.py."""

    def V(X):
        x, y = X[..., 0], X[..., 1]
        r1 = x * x + c * y * y
        r2 = c * (x - 1)**2 + (y - 1)**2
        r3 = fx * (0.25 * x**4 - x**3 / 3.) + fy * (0.25 * y**4 - y**3 / 3.)
        return r1 * r2 + r3

    def dV(X):
        x, y = X[..., 0], X[..., 1]
        r1 = x * x + c * y * y
        r2 = c * (x - 1)**2 + (y - 1)**2
        rval = np.empty_like(X, dtype=float)
        rval[..., 0] = r1 * 2 * c * (x - 1) + 2 * x * r2 + fx * x * x * (x - 1)
        rval[..., 1] = r1 * 2 * (y - 1) + 2 * c * y * r2 + fy * y * y * (y - 1)
        return rval

    return V, dV


def narrow_barrier(a=1.0, b=0.004, g=0.3):
    def V(X):
        x, y = X[..., 0], X[..., 1]
        return (a * (0.25 * x**4 - (1 + b) * x**3 / 3 + 0.5 * b * x**2)
                + (y - g * x**2)**2)

    def dV(X):
        x, y = X[..., 0], X[..., 1]
        rval = np.empty_like(X, dtype=float)
        rval[..., 0] = a * x * (x - b) * (x - 1) - 4 * g * x * (y - g * x**2)
        rval[..., 1] = 2 * (y - g * x**2)
        return rval

    return V, dV


def double_well(c, a=1.8, b=0.2):
    def V(X):
        x, y = X[..., 0], X[..., 1]
        return (x**2 + y**2) * (a * (x - 1)**2 + b * (y - 1)**2 - c)

    def dV(X):
        x, y = X[..., 0], X[..., 1]
        f = a * (x - 1)**2 + b * (y - 1)**2 - c
        r = x**2 + y**2
        rval = np.empty_like(X, dtype=float)
        rval[..., 0] = 2 * x * f + r * 2 * a * (x - 1)
        rval[..., 1] = 2 * y * f + r * 2 * b * (y - 1)
        return rval

    return V, dV


def many_fields(ndim):
    a = 0.04 * np.arange(1, ndim)
    m2 = 2.0 + 0.2 * np.arange(1, ndim)

    def V(X):
        phi = X[..., 0]
        d = X[..., 1:] - a * phi[..., None]**2
        return (0.4 * phi**2 - 1.6 * phi**3 + phi**4
                + 0.5 * np.sum(m2 * d * d, axis=-1))

    def dV(X):
        phi = X[..., 0]
        d = X[..., 1:] - a * phi[..., None]**2
        rval = np.empty_like(X, dtype=float)
        rval[..., 0] = (0.8 * phi - 4.8 * phi**2 + 4.0 * phi**3
                        - 2.0 * phi * np.sum(m2 * a * d, axis=-1))
        rval[..., 1:] = m2 * d
        return rval

    return V, dV, np.concatenate([[1.0], a])


# ---------------------------------------------------------------------------
# Comparison machinery
# ---------------------------------------------------------------------------

ROWS = []      # (case, quantity, ct value, ft value, rel. diff, rtol)
TIMES = []     # (case, t_ct, t_ft)


def timed(func):
    t0 = time.perf_counter()
    out = func()
    return out, time.perf_counter() - t0


def quiet(func):
    """Run func with stdout suppressed: cosmoTransitions always prints
    the progress of the path deformation."""
    with contextlib.redirect_stdout(io.StringIO()):
        return func()


def compare(case, quantity, ct_val, ft_val, rtol):
    ref = abs(ct_val) if ct_val != 0 else 1.0
    ROWS.append((case, quantity, ct_val, ft_val, abs(ft_val - ct_val) / ref,
                 rtol))


def run_1d(case, V, dV, rtol, **profile_params):
    def run_cosmo():
        inst = ct_t1d.SingleFieldInstanton(1.0, 0.0, V, dV)
        prof = inst.findProfile(**profile_params)
        return inst, prof, inst.findAction(prof)

    def run_fortran():
        inst = ft.SingleFieldInstanton(1.0, 0.0, V, dV)
        prof = inst.findProfile(**profile_params)
        return inst, prof, inst.findAction(prof)

    (ci, cp, ca), t_ct = timed(run_cosmo)
    (fi, fp, fa), t_ft = timed(run_fortran)
    TIMES.append((case, t_ct, t_ft))
    compare(case, "phi_bar", ci.phi_bar, fi.phi_bar, 1e-4)
    compare(case, "rscale", ci.rscale, fi.rscale, 1e-4)
    compare(case, "phi(r=0)", cp.Phi[0], fp.Phi[0], rtol)
    compare(case, "action", ca, fa, rtol)


def run_full(case, path, V, dV, rtol, **kwargs):
    ct_res, t_ct = timed(lambda: quiet(
        lambda: ct_pd.fullTunneling(path, V, dV, **kwargs)))
    ft_res, t_ft = timed(lambda: ft.fullTunneling(path, V, dV, **kwargs))
    TIMES.append((case, t_ct, t_ft))
    compare(case, "action", ct_res.action, ft_res.action, rtol)
    # The final paths are sampled at different points along the path, so
    # compare their end points (the two minima after extend_to_minima).
    for i, label in [(0, "start"), (-1, "end")]:
        dist = np.linalg.norm(ct_res.Phi[i] - ft_res.Phi[i])
        scale = np.linalg.norm(ct_res.Phi[0] - ct_res.Phi[-1])
        ROWS.append((case, "path {} (|d|/L)".format(label), 0.0, dist / scale,
                     dist / scale, 1e-2))
    return ct_res, ft_res


def print_report():
    print()
    print("{:<26} {:<18} {:>17} {:>17} {:>9} {:>7}  ".format(
        "case", "quantity", "cosmoTransitions", "fortransitions", "rel.diff",
        "rtol"))
    print("-" * 102)
    nfail = 0
    for case, q, cv, fv, rd, rtol in ROWS:
        ok = rd <= rtol
        nfail += not ok
        if q.startswith("path"):
            values = "{:>17} {:>17.3e}".format("", fv)
        else:
            values = "{:>17.10e} {:>17.10e}".format(cv, fv)
        print("{:<26} {:<18} {} {:>9.1e} {:>7.0e}  {}".format(
            case, q, values, rd, rtol, "ok" if ok else "FAIL"))
    print()
    print("{:<26} {:>12} {:>12} {:>9}".format(
        "case", "t_cosmo [s]", "t_fortr [s]", "speedup"))
    print("-" * 62)
    for case, t_ct, t_ft in TIMES:
        print("{:<26} {:>12.3f} {:>12.3f} {:>8.1f}x".format(
            case, t_ct, t_ft, t_ct / t_ft))
    t_ct = sum(t[1] for t in TIMES)
    t_ft = sum(t[2] for t in TIMES)
    print("{:<26} {:>12.3f} {:>12.3f} {:>8.1f}x".format(
        "total", t_ct, t_ft, t_ct / t_ft))
    return nfail


def main():
    quick = "--quick" in sys.argv[1:]
    print("cosmoTransitions {} ({})".format(
        getattr(cosmoTransitions, "__version__", "?"),
        cosmoTransitions.__file__))
    print("numpy {}, Python {}".format(np.__version__,
                                       sys.version.split()[0]))
    tight = dict(xtol=1e-9, phitol=1e-9)

    # test/check.f90 (default tolerances: agreement to a few 1e-3)
    run_1d("1d thin", v1_thin, dv1_thin, 5e-3)
    run_1d("1d thick", v1_thick, dv1_thick, 5e-3)
    run_1d("1d thin, tight tol.", v1_thin, dv1_thin, 1e-3, **tight)
    run_1d("1d thick, no dV", v1_thick, None, 5e-3)
    V, dV = example_2d(fy=2.0)
    run_full("2d thin", [[1., 1.], [0., 0.]], V, dV, 5e-3)
    V, dV = example_2d(fy=80.0)
    run_full("2d thick", [[1., 1.], [0., 0.]], V, dV, 5e-3)
    V, dV = narrow_barrier()
    run_full("narrow barrier, direct", [[1., 0.3], [0., 0.]], V, dV, 1e-3,
             V_spline_samples=None, tunneling_findProfile_params=tight)

    # test/double_well.f90 (tight tolerances)
    for c in ([0.05, 0.2] if quick else
              [0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3]):
        V, dV = double_well(c)
        xmin = optimize.minimize(V, [1.0, 1.0], jac=dV, method="BFGS",
                                 tol=1e-14).x
        run_full("double well c={:.2f}".format(c), [xmin, [0., 0.]], V, dV,
                 2e-3, tunneling_findProfile_params=tight)

    # test/many_fields.f90 (tight tolerances; cosmoTransitions is limited
    # to 10 fields)
    for ndim in ([2, 5] if quick else range(2, 11)):
        V, dV, x_true = many_fields(ndim)
        t = np.linspace(0.0, 1.0, 41)[:, None]
        run_full("many fields ndim={}".format(ndim), x_true * (1.0 - t), V,
                 dV, 2e-3, tunneling_findProfile_params=tight)

    nfail = print_report()
    print()
    if nfail:
        print("FAILED: {} of {} comparisons exceed their tolerance.".format(
            nfail, len(ROWS)))
        sys.exit(1)
    print("All {} comparisons within tolerance.".format(len(ROWS)))


if __name__ == "__main__":
    main()
