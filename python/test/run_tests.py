"""Tests of the Python interface pyfortransitions.

The physics tests repeat the Fortran test programs (test/check.f90,
test/double_well.f90, test/many_fields.f90) through the Python interface,
with the potentials written in Python, and assert against the same
reference values (computed with cosmoTransitions, see reference/) with the
same relative tolerances. The remaining tests cover the error handling,
argument validation and memory management of the interface.

Run with ``pixi run test`` or ``python test/run_tests.py`` (optionally with
``--quick`` to skip the slow many-field scan above 10 fields).
"""

import resource
import sys
import time

import numpy as np

import pyfortransitions as pf

NFAIL = 0
NPASS = 0


def _fail(msg):
    global NFAIL
    NFAIL += 1
    print("FAIL " + msg)


def _pass(msg):
    global NPASS
    NPASS += 1
    print("ok   " + msg)


def assert_close(name, val, ref, rtol):
    if abs(val - ref) > rtol * abs(ref):
        _fail("{}: got {:.8e}, expected {:.8e}, rtol {:.1e}".format(
            name, val, ref, rtol))
    else:
        _pass("{}: got {:.8e} (ref {:.8e})".format(name, val, ref))


def assert_true(name, cond):
    if cond:
        _pass(name)
    else:
        _fail(name)


def assert_raises(name, exc_type, func, check=None):
    try:
        func()
    except exc_type as err:
        if check is not None and not check(err):
            _fail("{}: {} raised with unexpected content: {!r}".format(
                name, exc_type.__name__, err))
        else:
            _pass("{}: {}: {}".format(name, type(err).__name__, err))
    except BaseException as err:  # noqa: B902
        _fail("{}: expected {}, got {}: {}".format(
            name, exc_type.__name__, type(err).__name__, err))
    else:
        _fail("{}: expected {}, nothing raised".format(name, exc_type.__name__))


# ---------------------------------------------------------------------------
# Potentials (as in the Fortran tests)
# ---------------------------------------------------------------------------

def v1_thin(phi):
    return 0.25 * phi**4 - 0.49 * phi**3 + 0.235 * phi**2


def dv1_thin(phi):
    return phi * (phi - 0.47) * (phi - 1.0)


def v2_thick(phi):
    return 0.25 * phi**4 - 0.4 * phi**3 + 0.1 * phi**2


def dv2_thick(phi):
    return phi * (phi - 0.2) * (phi - 1.0)


def example_potential(c=5.0, fx=10.0, fy=10.0):
    """The two-field potential of examples/fullTunneling.py, written in the
    vectorized cosmoTransitions style (X[..., i]), so that dV can be
    evaluated on all points of the path at once."""

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


def narrow_barrier_potential(a=1.0, b=0.004, g=0.3):
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


def double_well_potential(c, a=1.8, b=0.2):
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


def many_fields_potential(ndim):
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
        g = np.empty_like(X, dtype=float)
        g[..., 0] = (0.8 * phi - 4.8 * phi**2 + 4.0 * phi**3
                     - 2.0 * phi * np.sum(m2 * a * d, axis=-1))
        g[..., 1:] = m2 * d
        return g

    return V, dV, np.concatenate([[1.0], a])


# ---------------------------------------------------------------------------
# Tests repeating test/check.f90
# ---------------------------------------------------------------------------

def test_1d():
    inst = pf.SingleFieldInstanton(1.0, 0.0, v1_thin, dv1_thin)
    assert_close("1d thin: phi_bar", inst.phi_bar, 8.37171431e-1, 1e-4)
    assert_close("1d thin: rscale", inst.rscale, 1.66770894, 1e-4)
    prof = inst.findProfile()
    action = inst.findAction(prof)
    assert_close("1d thin: action", action, 1.0927555625e3, 5e-3)
    assert_true("1d thin: profile arrays consistent",
                len(prof.R) == len(prof.Phi) == len(prof.dPhi) > 500)
    assert_true("1d thin: Rerr is None", prof.Rerr is None)
    s, s_pot, s_kin = inst.findActionDerrick(prof)
    assert_true("1d thin: findActionDerrick action", s == action)
    # With the default tolerances (1e-4), the Derrick estimates of this
    # thin-walled profile are off by up to ~1%; with tight tolerances, all
    # three agree.
    prof_tight = inst.findProfile(xtol=1e-9, phitol=1e-9)
    s, s_pot, s_kin = inst.findActionDerrick(prof_tight)
    assert_close("1d thin, tight tolerances: action", s, 1.0927555625e3,
                 5e-3)
    assert_close("1d thin, tight tolerances: Derrick (potential term)",
                 s_pot, s, 1e-3)
    assert_close("1d thin, tight tolerances: Derrick (kinetic term)",
                 s_kin, s, 1e-3)

    # Finite-difference derivatives (no dV given)
    inst_fd = pf.SingleFieldInstanton(1.0, 0.0, v1_thin)
    action_fd = inst_fd.findAction(inst_fd.findProfile())
    assert_close("1d thin, FD derivatives: action", action_fd,
                 1.0927555625e3, 5e-3)

    # findAction on a modified profile (a plain namedtuple)
    prof2 = pf.Profile1D(prof.R, prof.Phi, 2 * prof.dPhi, prof.Rerr)
    assert_true("1d thin: findAction on modified profile differs",
                inst.findAction(prof2) != action)

    inst = pf.SingleFieldInstanton(1.0, 0.0, v2_thick, dv2_thick)
    prof = inst.findProfile()
    assert_close("1d thick: action", inst.findAction(prof),
                 6.6489883968, 5e-3)
    assert_close("1d thick: phi(0)", prof.Phi[0], 7.41882997e-1, 5e-3)

    # Given phi_bar and rscale are used as such
    inst = pf.SingleFieldInstanton(1.0, 0.0, v1_thin, dv1_thin,
                                   phi_bar=0.83717, rscale=1.6677)
    assert_true("1d thin: given phi_bar/rscale kept",
                inst.phi_bar == 0.83717 and inst.rscale == 1.6677)
    assert_close("1d thin, given phi_bar/rscale: action",
                 inst.findAction(inst.findProfile()), 1.0927555625e3, 5e-3)


def test_2d():
    V, dV = example_potential(c=5.0, fx=0.0, fy=2.0)
    res = pf.fullTunneling([[1.0, 1.0], [0.0, 0.0]], V, dV, verbose=True)
    assert_close("2d thin: action", res.action, 1.7669376756e3, 5e-3)
    assert_close("2d thin: path start x", res.Phi[0, 0], 1.0, 1e-2)
    assert_close("2d thin: path start y", res.Phi[0, 1], 1.0, 1e-2)
    assert_true("2d thin: path end x", abs(res.Phi[-1, 0]) < 1e-2)
    assert_true("2d thin: path end y", abs(res.Phi[-1, 1]) < 1e-2)
    assert_true("2d thin: Phi shape",
                res.Phi.shape == (len(res.profile1D.R), 2))
    profile1D, Phi, action, fRatio, saved_steps = res
    assert_true("2d thin: unpacks into 5 cosmoTransitions fields",
                action == res.action and saved_steps == [])

    V, dV = example_potential(c=5.0, fx=0.0, fy=80.0)
    res = pf.fullTunneling([[1.0, 1.0], [0.0, 0.0]], V, dV, verbose=True)
    assert_close("2d thick: action", res.action, 4.5036661308, 5e-3)


def test_narrow_barrier():
    V, dV = narrow_barrier_potential()
    path = [[1.0, 0.3], [0.0, 0.0]]
    params = dict(xtol=1e-9, phitol=1e-9)
    try:
        res = pf.fullTunneling(path, V, dV,
                               tunneling_findProfile_params=params)
        print("narrow barrier, spline (100 samples): action {:.8e}, "
              "Derrick {:.8e} {:.8e}".format(res.action, res.action_pot,
                                            res.action_kin))
    except pf.FortransitionsError as err:
        print("narrow barrier, spline (100 samples): {}: {}".format(
            type(err).__name__, err))
    res = pf.fullTunneling(path, V, dV, V_spline_samples=None,
                           tunneling_findProfile_params=params)
    assert_true("narrow barrier, direct: converged", res.converged)
    assert_close("narrow barrier, direct: action", res.action,
                 1.11271076e-2, 1e-3)
    assert_close("narrow barrier, direct: Derrick (potential term)",
                 res.action_pot, res.action, 1e-3)
    assert_close("narrow barrier, direct: Derrick (kinetic term)",
                 res.action_kin, res.action, 1e-3)


# ---------------------------------------------------------------------------
# Tests repeating test/double_well.f90 and test/many_fields.f90
# ---------------------------------------------------------------------------

def test_double_well():
    c_vals = [0.01, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30]
    # True minima from test/double_well.f90 (found there with gradmin)
    ref_xmin = [1.002674189573, 1.011663435734, 1.020199200731,
                1.026807410988, 1.032131164050, 1.036548270078,
                1.040296717098]
    ref_ymin = [1.024593855202, 1.115773451273, 1.216831369462,
                1.307134557354, 1.389232396274, 1.464850194264,
                1.535207198375]
    ref_action = [7.8805837752e3, 3.8913237834e2, 1.1963573348e2,
                  6.2230671791e1, 3.9636887635e1, 2.7922475045e1,
                  2.0887209824e1]
    for c, x, y, ref in zip(c_vals, ref_xmin, ref_ymin, ref_action):
        V, dV = double_well_potential(c)
        assert_true("[c={:.2f}] grad V(true minimum) ~ 0".format(c),
                    np.max(np.abs(dV(np.array([x, y])))) < 1e-9)
        res = pf.fullTunneling(
            [[x, y], [0.0, 0.0]], V, dV,
            tunneling_findProfile_params=dict(xtol=1e-9, phitol=1e-9))
        assert_close("[c={:.2f}] bounce action S_3".format(c), res.action,
                     ref, 2e-3)


def test_many_fields(ndim_max):
    ref_action = {
        2: 3.3399851500, 3: 3.3265446619, 4: 3.3478390266,
        5: 3.3864753646, 6: 3.4488214698, 7: 3.5359964603,
        8: 3.6595890944, 9: 3.8232885005, 10: 4.0330848003}
    action_prev = 0.0
    for ndim in range(2, ndim_max + 1):
        V, dV, x_true = many_fields_potential(ndim)
        t = np.linspace(0.0, 1.0, 41)[:, None]
        path = x_true * (1.0 - t)
        res = pf.fullTunneling(
            path, V, dV,
            tunneling_findProfile_params=dict(xtol=1e-9, phitol=1e-9))
        label = "[ndim={}] ".format(ndim)
        assert_true(label + "deformation converged", res.converged)
        if ndim in ref_action:
            assert_close(label + "bounce action S_3", res.action,
                         ref_action[ndim], 2e-3)
        assert_close(label + "Derrick action (potential term)",
                     res.action_pot, res.action, 1e-3)
        assert_close(label + "Derrick action (kinetic term)",
                     res.action_kin, res.action, 1e-3)
        assert_true(label + "action grows with ndim",
                    ndim <= 3 or res.action > action_prev)
        action_prev = res.action


# ---------------------------------------------------------------------------
# Error handling, validation, memory
# ---------------------------------------------------------------------------

def test_errors():
    assert_raises(
        "stable minimum", pf.PotentialError,
        lambda: pf.SingleFieldInstanton(0.0, 1.0, v1_thin, dv1_thin),
        check=lambda e: e.args[1] == "stable, not metastable")
    assert_raises(
        "no barrier", pf.PotentialError,
        lambda: pf.SingleFieldInstanton(1.0, 0.0, lambda p: -p),
        check=lambda e: e.args[1] == "no barrier")
    assert_true("PotentialError is a FortransitionsError",
                issubclass(pf.PotentialError, pf.FortransitionsError))

    # Exceptions raised inside the potential propagate unchanged, also when
    # they occur in the middle of a calculation.
    class Boom(Exception):
        pass

    def failing_after(n, f, exc=Boom):
        count = [0]

        def g(*args):
            count[0] += 1
            if count[0] > n:
                raise exc("raised in call {}".format(count[0]))
            return f(*args)

        return g, count

    v, _ = failing_after(0, v1_thin)
    assert_raises("exception in V during init", Boom,
                  lambda: pf.SingleFieldInstanton(1.0, 0.0, v, dv1_thin))

    v, count = failing_after(200, v1_thin)
    inst = pf.SingleFieldInstanton(1.0, 0.0, v, dv1_thin)
    assert_raises("exception in V during findProfile", Boom,
                  inst.findProfile)
    calls = count[0]
    assert_true("user V not called again after the exception "
                "({} calls)".format(calls), calls == 201)

    V, dV = example_potential(c=5.0, fx=0.0, fy=80.0)
    for vectorized, k in [(None, 10), (True, 10), (False, 50)]:
        dv, count = failing_after(k, dV)
        t0 = time.time()
        assert_raises(
            "exception in dV during fullTunneling (vectorized={})".format(
                vectorized), Boom,
            lambda: pf.fullTunneling([[1, 1], [0, 0]], V, dv,
                                     vectorized=vectorized))
        assert_true("fullTunneling stops quickly after the exception "
                    "(vectorized={})".format(vectorized),
                    time.time() - t0 < 5.0 and count[0] == k + 1)

    v, _ = failing_after(10, V, exc=KeyboardInterrupt)
    assert_raises("KeyboardInterrupt in V", KeyboardInterrupt,
                  lambda: pf.fullTunneling([[1, 1], [0, 0]], v, dV))

    assert_raises("dV with wrong number of components", ValueError,
                  lambda: pf.fullTunneling([[1, 1], [0, 0]], V,
                                           lambda X: np.zeros(3)))

    # The instanton still works after a failed call
    inst = pf.SingleFieldInstanton(1.0, 0.0, v1_thin, dv1_thin)
    assert_close("instanton usable after errors",
                 inst.findAction(inst.findProfile()), 1.0927555625e3, 5e-3)


def test_vectorization():
    from pyfortransitions._callbacks import _Vectorization

    V, dV = example_potential(c=5.0, fx=0.0, fy=2.0)
    path = [[1.0, 1.0], [0.0, 0.0]]
    calls = {"n": 0, "batched": 0}

    def counting_dV(X):
        calls["n"] += 1
        calls["batched"] += np.ndim(X) == 2
        return dV(X)

    res_auto = pf.fullTunneling(path, V, counting_dV)
    n_auto, n_batched = calls["n"], calls["batched"]
    calls.update(n=0, batched=0)
    res_off = pf.fullTunneling(path, V, counting_dV, vectorized=False)
    n_off = calls["n"]
    res_on = pf.fullTunneling(path, V, dV, vectorized=True)
    # numpy may round differently in the last bit when evaluating an array
    # instead of a single point, so the results agree to rounding only.
    assert_true("vectorized auto and True identical",
                res_auto.action == res_on.action)
    assert_close("vectorized auto vs. False (point by point)",
                 res_auto.action, res_off.action, 1e-10)
    assert_true("vectorized auto: dV batched ({} of {} calls, vs. {} "
                "calls point by point)".format(n_batched, n_auto, n_off),
                n_batched > 10 and n_auto < n_off / 100)
    assert_close("2d thin, vectorized: action", res_auto.action,
                 1.7669376756e3, 5e-3)

    # A potential written for single points only (plain floats) is
    # detected and evaluated point by point.
    c, fy = 5.0, 2.0

    def V_float(X):
        x, y = float(X[0]), float(X[1])
        return (x * x + c * y * y) * (c * (x - 1)**2 + (y - 1)**2) \
            + fy * (0.25 * y**4 - y**3 / 3.)

    def dV_float(X):
        x, y = float(X[0]), float(X[1])
        r1 = x * x + c * y * y
        r2 = c * (x - 1)**2 + (y - 1)**2
        return [r1 * 2 * c * (x - 1) + 2 * x * r2,
                r1 * 2 * (y - 1) + 2 * c * y * r2 + fy * y * y * (y - 1)]

    res_float = pf.fullTunneling(path, V_float, dV_float)
    assert_close("single-point dV detected, same action",
                 res_float.action, res_off.action, 1e-10)
    assert_raises("single-point dV with vectorized=True",
                  (TypeError, ValueError),
                  lambda: pf.fullTunneling(path, V_float, dV_float,
                                           vectorized=True))
    assert_raises("invalid vectorized value", ValueError,
                  lambda: pf.fullTunneling(path, V, dV, vectorized="yes"))

    # A non-vectorized dV that happens to return the right shape (here for
    # npts == ndim) must not be mistaken for a vectorized one.
    def swap(X):
        x, y = X  # unpacks the first axis: wrong for X of shape (npts, 2)
        return np.array([y, x])

    vec = _Vectorization(None)
    X = np.array([[1.0, 2.0], [3.0, 4.0]])
    out = vec.gradients(swap, X)
    assert_true("same-shape non-vectorized dV detected",
                vec.mode is False
                and np.array_equal(out, [[2.0, 1.0], [4.0, 3.0]]))
    vec = _Vectorization(None)
    out = vec.gradients(lambda X: 2 * X, X)
    assert_true("vectorized dV detected",
                vec.mode is True and np.array_equal(out, 2 * X))


def test_validation():
    V, dV = example_potential(c=5.0, fx=0.0, fy=80.0)
    path = [[1, 1], [0, 0]]
    assert_raises("maxiter=0", ValueError,
                  lambda: pf.fullTunneling(path, V, dV, maxiter=0))
    assert_raises("V_spline_samples=0", ValueError,
                  lambda: pf.fullTunneling(path, V, dV, V_spline_samples=0))
    assert_raises("bad path shape", ValueError,
                  lambda: pf.fullTunneling([1.0, 0.0], V, dV))
    assert_raises("save_all_steps", NotImplementedError,
                  lambda: pf.fullTunneling(path, V, dV, save_all_steps=True))
    assert_raises("callback", NotImplementedError,
                  lambda: pf.fullTunneling(path, V, dV,
                                           callback=lambda *a: True))
    assert_raises("unsupported findProfile key", NotImplementedError,
                  lambda: pf.fullTunneling(
                      path, V, dV,
                      tunneling_findProfile_params=dict(thinCutoff=0.1)))
    assert_raises("negative alpha", ValueError,
                  lambda: pf.SingleFieldInstanton(1.0, 0.0, v1_thin,
                                                  alpha=-1))
    inst = pf.SingleFieldInstanton(1.0, 0.0, v1_thin, dv1_thin)
    assert_raises("npoints=1", ValueError,
                  lambda: inst.findProfile(npoints=1))

    # Supported parameters are passed through
    res = pf.fullTunneling(
        path, V, dV, maxiter=1,
        tunneling_init_params=dict(alpha=3),
        tunneling_findProfile_params=dict(npoints=200),
        deformation_init_params=dict(nb=8, kb=3),
        deformation_deform_params=dict(verbose=False))
    assert_true("maxiter=1: one iteration", res.num_iters == 1)
    assert_true("npoints=200 used", 200 <= len(res.profile1D.R) <= 300)
    res4 = pf.fullTunneling(path, V, dV,
                            tunneling_init_params=dict(alpha=3))
    assert_true("alpha=3 (O(4)) action differs from alpha=2",
                abs(res4.action - 4.5036661308) > 0.1)


def test_memory():
    """Repeated calls must not leak memory (handles are freed)."""

    def maxrss_mb():
        return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0

    V, dV = example_potential(c=5.0, fx=0.0, fy=80.0)
    for _ in range(3):  # warm-up
        pf.fullTunneling([[1, 1], [0, 0]], V, dV)
        inst = pf.SingleFieldInstanton(1.0, 0.0, v1_thin, dv1_thin)
        inst.findProfile()
    rss0 = maxrss_mb()
    for _ in range(300):
        inst = pf.SingleFieldInstanton(1.0, 0.0, v1_thin, dv1_thin)
        prof = inst.findProfile()
        inst.findAction(prof)
    for _ in range(20):
        pf.fullTunneling([[1, 1], [0, 0]], V, dV)
    growth = maxrss_mb() - rss0
    assert_true("memory growth over 300 instantons and 20 fullTunneling "
                "calls: {:.1f} MB".format(growth), growth < 20.0)


def main():
    quick = "--quick" in sys.argv[1:]
    tests = [
        ("check: 1d", test_1d),
        ("check: 2d", test_2d),
        ("check: narrow barrier", test_narrow_barrier),
        ("double_well", test_double_well),
        ("many_fields", lambda: test_many_fields(10 if quick else 20)),
        ("vectorization", test_vectorization),
        ("errors", test_errors),
        ("validation", test_validation),
        ("memory", test_memory),
    ]
    for name, test in tests:
        print("--- {} ---".format(name))
        t0 = time.time()
        try:
            test()
        except BaseException as err:  # noqa: B902
            _fail("{}: unexpected {}: {}".format(name, type(err).__name__,
                                                 err))
        print("    ({:.1f} s)".format(time.time() - t0))
    print()
    if NFAIL:
        print("FAILED: {} of {} test(s) failed.".format(NFAIL, NFAIL + NPASS))
        sys.exit(1)
    print("All {} tests passed.".format(NPASS))


if __name__ == "__main__":
    main()
