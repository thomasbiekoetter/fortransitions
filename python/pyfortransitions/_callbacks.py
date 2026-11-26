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


def nd_callbacks(V, dV, guard):
    """C callbacks for a multi-field potential V(X) and its gradient dV(X).

    X is passed as a fresh 1d array of shape (ndim,), so potentials written
    for cosmoTransitions (using e.g. ``X[..., 0]``) work unchanged.
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

    return _lib.V_ND_CB(v_cb), _lib.GRAD_ND_CB(grad_cb)


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
