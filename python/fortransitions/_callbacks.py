"""Wrapping of Python potentials as C callbacks for the Fortran code.

ctypes cannot propagate exceptions out of a callback: it prints them and
returns 0, so the Fortran code would silently continue with wrong values.
The wrappers therefore catch every exception (including
KeyboardInterrupt), store it in a guard object and return NaN. After the
first exception, all further callbacks return NaN immediately without
calling the user function, so that the Fortran code finishes quickly. The
caller re-raises the stored exception with ``guard.check()`` after the
Fortran call returns.
"""

import ctypes
import math

import numpy as np

from . import _lib

_NAN = math.nan


class CallbackGuard:
    """Collects the first exception raised inside a callback."""

    def __init__(self):
        self.error = None

    def reset(self):
        self.error = None

    def check(self):
        """Re-raise a stored exception, if any."""
        if self.error is not None:
            err, self.error = self.error, None
            raise err


def _grad_per_point(dV, X):
    """dV evaluated point by point at the rows of X, shape (npts, ndim)."""
    npts, n = X.shape
    out = np.empty((npts, n))
    for i in range(npts):
        dv = np.asarray(dV(X[i].copy()), dtype=np.float64)
        if dv.size != n:
            raise ValueError(
                "dV(X) returned {} values for {} fields".format(dv.size, n))
        out[i] = dv.reshape(n)
    return out


class _Vectorization:
    """Decides whether dV can be called with many points at once.

    mode None (auto): the first batched call checks that dV(X) with X of
    shape (npts, ndim) returns an array of the same shape whose rows agree
    with point-by-point evaluations at two sample points. If not (or if the
    call raises an ordinary exception), dV is evaluated point by point from
    then on. mode True: always batched (a wrong shape is an error). mode
    False: never batched (handled by passing a null callback).
    """

    def __init__(self, mode):
        self.mode = mode

    def gradients(self, dV, X):
        npts, n = X.shape
        if self.mode is True:
            dv = np.asarray(dV(X), dtype=np.float64)
            if dv.shape != (npts, n):
                raise ValueError(
                    "vectorized dV(X) with X of shape {} returned shape {}; "
                    "expected {}".format(X.shape, dv.shape, (npts, n)))
            return dv
        if self.mode is False:
            return _grad_per_point(dV, X)
        # Auto-detection on the first call.
        try:
            dv = np.asarray(dV(X.copy()), dtype=np.float64)
            ok = dv.shape == (npts, n)
            if ok:
                for i in sorted({0, npts // 2}):
                    ref = np.asarray(dV(X[i].copy()), dtype=np.float64)
                    if ref.size != n or not np.allclose(
                            dv[i], ref.reshape(n), rtol=1e-10,
                            atol=1e-12 * max(1.0, np.max(np.abs(ref)))):
                        ok = False
                        break
        except Exception:
            ok = False
        self.mode = ok
        return dv if ok else _grad_per_point(dV, X)


def nd_callbacks(V, dV, guard, vectorized=None):
    """C callbacks for a multi-field potential V(X) and its gradient dV(X).

    For single points, X is passed as a fresh 1d array of shape (ndim,), so
    potentials written for cosmoTransitions (using e.g. ``X[..., 0]``) work
    unchanged. The third callback evaluates the gradient on all points of
    the path at once (X of shape (npts, ndim)); see _Vectorization for the
    meaning of `vectorized`. With vectorized=False it is a null pointer, and
    the Fortran code calls the single-point callback for each point.
    """

    def v_cb(n, x_ptr):
        if guard.error is not None:
            return _NAN
        try:
            x = np.ctypeslib.as_array(x_ptr, shape=(n,)).copy()
            return float(V(x))
        except BaseException as err:  # noqa: B902 - see module docstring
            guard.error = err
            return _NAN

    def grad_cb(n, x_ptr, g_ptr):
        g = np.ctypeslib.as_array(g_ptr, shape=(n,))
        if guard.error is not None:
            g[:] = _NAN
            return
        try:
            x = np.ctypeslib.as_array(x_ptr, shape=(n,)).copy()
            dv = np.asarray(dV(x), dtype=np.float64)
            if dv.size != n:
                raise ValueError(
                    "dV(X) returned {} values for {} fields".format(dv.size, n))
            g[:] = dv.reshape(n)
        except BaseException as err:  # noqa: B902
            guard.error = err
            g[:] = _NAN

    vec = _Vectorization(vectorized)

    def grad_many_cb(npts, n, x_ptr, g_ptr):
        # Column-major (npts, n) arrays, viewed as (npts, n) numpy arrays.
        g = np.ctypeslib.as_array(g_ptr, shape=(n, npts)).T
        if guard.error is not None:
            g[:] = _NAN
            return
        try:
            X = np.ascontiguousarray(
                np.ctypeslib.as_array(x_ptr, shape=(n, npts)).T)
            g[:] = vec.gradients(dV, X)
        except BaseException as err:  # noqa: B902
            guard.error = err
            g[:] = _NAN

    grad_many = (_lib.GRAD_MANY_CB() if vectorized is False
                 else _lib.GRAD_MANY_CB(grad_many_cb))
    return _lib.V_ND_CB(v_cb), _lib.GRAD_ND_CB(grad_cb), grad_many


def f1d_callback(f, guard):
    """C callback for a function f(phi) of one field (V, dV or d2V).

    Returns a null function pointer if f is None, which the Fortran side
    interprets as "use finite differences".
    """
    if f is None:
        return _lib.F_1D_CB()

    def cb(phi):
        if guard.error is not None:
            return _NAN
        try:
            return float(f(phi))
        except BaseException as err:  # noqa: B902
            guard.error = err
            return _NAN

    return _lib.F_1D_CB(cb)


def out_double():
    return ctypes.c_double()


def out_int():
    return ctypes.c_int()


def out_handle():
    return ctypes.c_void_p()
