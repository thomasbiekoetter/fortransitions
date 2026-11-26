"""One-dimensional instantons, mirroring cosmoTransitions.tunneling1D."""

import ctypes
from collections import namedtuple

import numpy as np

from . import _lib
from ._callbacks import CallbackGuard
from ._callbacks import f1d_callback
from .exceptions import PotentialError  # noqa: F401 (re-exported)
from .exceptions import raise_for_status

Profile1D = namedtuple("Profile1D", "R Phi dPhi Rerr")
Profile1D.__doc__ = """The bubble profile.

R, Phi, dPhi : 1d arrays
    Radii, field values and field derivatives of the profile.
Rerr : float or None
    The first radius at which the integration step size fell below its
    minimum, or None if that never happened.
"""


def _profile_from_handle(prof_handle):
    """Copy a Fortran profile into a Profile1D and free the handle."""
    try:
        n = ctypes.c_int()
        _lib.profile_size(prof_handle, ctypes.byref(n))
        r = np.empty(n.value)
        phi = np.empty(n.value)
        dphi = np.empty(n.value)
        rerr = ctypes.c_double()
        _lib.profile_arrays(prof_handle, n.value, r, phi, dphi,
                            ctypes.byref(rerr))
    finally:
        _lib.profile_free(prof_handle)
    return Profile1D(r, phi, dphi, rerr.value if rerr.value >= 0 else None)


class SingleFieldInstanton:
    """Instanton of a single scalar field (without gravity), found with the
    overshoot/undershoot method.

    Python interface to the Fortran port of
    cosmoTransitions.tunneling1D.SingleFieldInstanton, with the same
    arguments and defaults.

    Parameters
    ----------
    phi_absMin : float
        The field value at the stable vacuum to which the instanton
        tunnels.
    phi_metaMin : float
        The field value at the metastable vacuum.
    V : callable
        The potential V(phi), called with a float.
    dV, d2V : callable, optional
        First and second derivatives of the potential. If not given, they
        are computed with fourth-order finite differences with step
        ``phi_eps * abs(phi_absMin - phi_metaMin)``.
    phi_eps : float, optional
        Relative step for the finite differences (default 1e-3).
    alpha : float, optional
        Coefficient of the friction term in the ODE, i.e. the number of
        spacetime dimensions minus 1 (default 2).
    phi_bar : float, optional
        Field value at the edge of the barrier. Found automatically if not
        given.
    rscale : float, optional
        Approximate radial scale of the instanton. Found automatically if
        not given.

    Raises
    ------
    PotentialError
        If ``V(phi_metaMin) <= V(phi_absMin)`` ("stable, not metastable")
        or if there is no barrier between the minima ("no barrier").
    """

    def __init__(self, phi_absMin, phi_metaMin, V, dV=None, d2V=None,
                 phi_eps=1e-3, alpha=2, phi_bar=None, rscale=None):
        self._handle = None
        if not alpha >= 0:
            raise ValueError("alpha must be non-negative")
        self.phi_absMin = float(phi_absMin)
        self.phi_metaMin = float(phi_metaMin)
        self.alpha = float(alpha)
        self.V = V
        self._guard = CallbackGuard()
        # Keep the ctypes callbacks alive as long as the Fortran object can
        # call them.
        self._callbacks = (f1d_callback(V, self._guard),
                           f1d_callback(dV, self._guard),
                           f1d_callback(d2V, self._guard))
        handle = ctypes.c_void_p()
        status = ctypes.c_int()
        _lib.sfi_new(
            self.phi_absMin, self.phi_metaMin, *self._callbacks,
            self.alpha, float(phi_eps),
            int(phi_bar is not None),
            float(phi_bar) if phi_bar is not None else 0.0,
            int(rscale is not None),
            float(rscale) if rscale is not None else 0.0,
            ctypes.byref(handle), ctypes.byref(status))
        self._handle = handle.value
        if self._handle is not None and (self._guard.error is not None
                                         or status.value != 0):
            self._free()
        self._guard.check()
        raise_for_status(status.value, "SingleFieldInstanton")
        self.phi_eps = float(phi_eps) * abs(self.phi_absMin - self.phi_metaMin)
        phi_bar_out = ctypes.c_double()
        rscale_out = ctypes.c_double()
        _lib.sfi_params(self._handle, ctypes.byref(phi_bar_out),
                        ctypes.byref(rscale_out))
        self.phi_bar = phi_bar_out.value
        self.rscale = rscale_out.value

    def findProfile(self, xguess=None, xtol=1e-4, phitol=1e-4,
                    thinCutoff=.01, npoints=500, rmin=1e-4, rmax=1e4,
                    max_interior_pts=None):
        """Calculate the bubble profile by iteratively over/undershooting.

        The arguments are those of
        cosmoTransitions.tunneling1D.SingleFieldInstanton.findProfile;
        `rmin` and `rmax` are relative to ``self.rscale``.

        Returns
        -------
        Profile1D
        """
        if int(npoints) < 2:
            raise ValueError("npoints must be at least 2")
        if max_interior_pts is not None and int(max_interior_pts) < 0:
            raise ValueError("max_interior_pts must be non-negative")
        prof_handle = ctypes.c_void_p()
        status = ctypes.c_int()
        self._guard.reset()
        _lib.sfi_find_profile(
            self._handle,
            int(xguess is not None),
            float(xguess) if xguess is not None else 0.0,
            float(xtol), float(phitol), float(thinCutoff), int(npoints),
            float(rmin), float(rmax),
            int(max_interior_pts is not None),
            int(max_interior_pts) if max_interior_pts is not None else 0,
            ctypes.byref(prof_handle), ctypes.byref(status))
        if prof_handle.value is not None and self._guard.error is not None:
            _lib.profile_free(prof_handle.value)
        self._guard.check()
        raise_for_status(status.value, "findProfile")
        return _profile_from_handle(prof_handle.value)

    def _action_terms(self, profile):
        r = np.ascontiguousarray(profile.R, dtype=np.float64)
        phi = np.ascontiguousarray(profile.Phi, dtype=np.float64)
        dphi = np.ascontiguousarray(profile.dPhi, dtype=np.float64)
        if not (r.ndim == phi.ndim == dphi.ndim == 1
                and r.size == phi.size == dphi.size):
            raise ValueError("profile.R, profile.Phi and profile.dPhi must be "
                             "1d arrays of equal length")
        rerr = getattr(profile, "Rerr", None)
        action = ctypes.c_double()
        action_pot = ctypes.c_double()
        action_kin = ctypes.c_double()
        self._guard.reset()
        _lib.sfi_find_action(
            self._handle, r.size, r, phi, dphi,
            float(rerr) if rerr is not None else -1.0,
            ctypes.byref(action), ctypes.byref(action_pot),
            ctypes.byref(action_kin))
        self._guard.check()
        return action.value, action_pot.value, action_kin.value

    def findAction(self, profile):
        """The Euclidean action of the instanton for a given profile.

        S = int [ (dphi/dr)^2/2 + V(phi) - V(phi_metaMin) ] r^alpha dr dOmega
        """
        return self._action_terms(profile)[0]

    def findActionDerrick(self, profile):
        """The action together with its Derrick-theorem estimates.

        Returns ``(action, action_pot, action_kin)``, where `action_pot` and
        `action_kin` are the action computed from only the potential or only
        the kinetic term of the integrand, rescaled by Derrick's theorem.
        For an exact solution of the bounce equation all three agree, so
        their spread measures the accuracy of the profile. (Not part of
        cosmoTransitions.)
        """
        return self._action_terms(profile)

    def _free(self):
        if self._handle is not None:
            _lib.sfi_free(self._handle)
            self._handle = None

    def __del__(self):
        try:
            self._free()
        except Exception:  # interpreter shutdown
            pass
