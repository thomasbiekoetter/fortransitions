"""Multi-field instantons via path deformation, mirroring
cosmoTransitions.pathDeformation.fullTunneling."""

import ctypes
from collections import namedtuple

import numpy as np

from . import _lib
from ._callbacks import CallbackGuard
from ._callbacks import nd_callbacks
from .exceptions import DeformationError  # noqa: F401 (re-exported)
from .exceptions import raise_for_status
from .tunneling1D import Profile1D

_rval_base = namedtuple("fullTunneling_rval",
                        "profile1D Phi action fRatio saved_steps")


class fullTunneling_rval(_rval_base):
    """Result of fullTunneling.

    The tuple fields are those of cosmoTransitions:

    profile1D : Profile1D
        The one-dimensional bubble profile along the final path. Its `Phi`
        is the distance along the path.
    Phi : array, shape (len(profile1D.R), ndim)
        The points of the final path in field space, one per profile point.
    action : float
        The Euclidean action of the instanton.
    fRatio : float
        The largest force perpendicular to the final path relative to the
        largest potential gradient (zero for a perfect solution).
    saved_steps : list
        Always empty (``save_all_steps`` is not supported).

    Extra attributes (not part of the tuple, so that unpacking into five
    values keeps working):

    action_pot, action_kin : float
        The action computed from only the potential or only the kinetic
        term of the integrand, rescaled by Derrick's theorem. For an exact
        solution both equal `action`; their spread measures the accuracy.
    num_iters : int
        Number of tunneling/deformation iterations performed.
    converged : bool
        Whether the deformation converged before `maxiter` was reached.
    """

    def __new__(cls, profile1D, Phi, action, fRatio, saved_steps,
                action_pot, action_kin, num_iters, converged):
        self = super().__new__(cls, profile1D, Phi, action, fRatio,
                               saved_steps)
        self.action_pot = action_pot
        self.action_kin = action_kin
        self.num_iters = num_iters
        self.converged = converged
        return self


def _take(params, name, allowed):
    """Copy of a parameter dict, rejecting keys the Fortran port lacks."""
    params = dict(params)
    unknown = set(params) - set(allowed)
    if unknown:
        raise NotImplementedError(
            "fullTunneling: {} {} not supported by fortransitions "
            "(supported: {})".format(name, sorted(unknown),
                                     sorted(allowed) or "none"))
    return params


def _check_class(cls, name, expected):
    if cls is not None and getattr(cls, "__name__", None) != expected:
        raise NotImplementedError(
            "fullTunneling: only the default {} ({}) is supported".format(
                name, expected))


def fullTunneling(path_pts, V, dV, maxiter=20, fixEndCutoff=.03,
                  save_all_steps=False, verbose=False,
                  callback=None, callback_data=None,
                  V_spline_samples=100,
                  tunneling_class=None,
                  tunneling_init_params={},
                  tunneling_findProfile_params={},
                  deformation_class=None,
                  deformation_init_params={},
                  deformation_deform_params={},
                  vectorized=None):
    """Calculate the instanton solution in multiple field dimensions.

    Python interface to the Fortran port of
    cosmoTransitions.pathDeformation.fullTunneling, with the same
    arguments and defaults (with the default tunneling_class
    SingleFieldInstanton and deformation_class Deformation_Spline).

    Parameters
    ----------
    path_pts : array_like, shape (num_points, ndim)
        Initial guess for the tunneling path. The first point should be at
        (or near) the stable minimum to which the field tunnels, the last
        point at the metastable minimum.
    V, dV : callable
        The potential V(X) and its gradient dV(X), called with a 1d array X
        of shape (ndim,); dV must return ndim values. In the path
        deformation, dV is also called with all points of the path at once,
        X of shape (npts, ndim), and must then return an array of the same
        shape, like in cosmoTransitions (see `vectorized`).
    maxiter : int, optional
        Maximum number of tunneling/deformation iterations (default 20).
    fixEndCutoff : float, optional
        Ignored (it is also unused by cosmoTransitions).
    verbose : bool, optional
        Print a message at the start of each iteration.
    V_spline_samples : int or None, optional
        Number of samples of V along the path used for the potential
        spline (default 100). None evaluates V and dV directly along the
        path, which resolves barriers that are narrow compared to the path
        length.
    tunneling_init_params : dict, optional
        Supported key: ``alpha``.
    tunneling_findProfile_params : dict, optional
        Supported keys: ``npoints``, ``xtol``, ``phitol``.
    deformation_init_params : dict, optional
        Supported keys: ``nb``, ``kb``, ``fix_start``, ``fix_end``.
    deformation_deform_params : dict, optional
        Supported key: ``verbose`` (default: the value of `verbose`).
    vectorized : bool or None, optional
        Whether dV can be called with many points at once (not part of
        cosmoTransitions). None (default) detects it on the first call: dV
        is used with X of shape (npts, ndim) if it returns an array of that
        shape that agrees with point-by-point evaluations; otherwise dV is
        evaluated point by point (slower). True always passes all points
        at once; False always evaluates point by point.

    The options `save_all_steps`, `callback`, other tunneling and
    deformation classes and other parameter keys are not supported and
    raise NotImplementedError.

    Returns
    -------
    fullTunneling_rval

    Raises
    ------
    PotentialError, IntegrationError, DeformationError, FortransitionsError
    """
    if save_all_steps:
        raise NotImplementedError(
            "fullTunneling: save_all_steps is not supported")
    if callback is not None:
        raise NotImplementedError("fullTunneling: callback is not supported")
    _check_class(tunneling_class, "tunneling_class", "SingleFieldInstanton")
    _check_class(deformation_class, "deformation_class", "Deformation_Spline")
    tip = _take(tunneling_init_params, "tunneling_init_params",
                ["alpha"])
    tfp = _take(tunneling_findProfile_params, "tunneling_findProfile_params",
                ["npoints", "xtol", "phitol"])
    dip = _take(deformation_init_params, "deformation_init_params",
                ["nb", "kb", "fix_start", "fix_end", "save_all_steps",
                       "v2min"])
    if dip.pop("save_all_steps", False):
        raise NotImplementedError(
            "fullTunneling: save_all_steps is not supported")
    if dip.pop("v2min", 0.0) != 0.0:
        raise NotImplementedError(
            "fullTunneling: deformation_init_params v2min is not supported")
    ddp = _take(deformation_deform_params, "deformation_deform_params",
                ["verbose"])

    pts = np.asfortranarray(path_pts, dtype=np.float64)
    if pts.ndim != 2 or pts.shape[0] < 2 or pts.shape[1] < 1:
        raise ValueError("path_pts must have shape (num_points, ndim) with "
                         "num_points >= 2")
    maxiter = int(maxiter)
    alpha = float(tip.get("alpha", 2.0))
    npoints = int(tfp.get("npoints", 500))
    nb = int(dip.get("nb", 10))
    kb = int(dip.get("kb", 3))
    if maxiter < 1:
        raise ValueError("maxiter must be positive")
    if V_spline_samples is not None and int(V_spline_samples) < 1:
        raise ValueError("V_spline_samples must be positive or None")
    if not alpha >= 0:
        raise ValueError("alpha must be non-negative")
    if npoints < 2:
        raise ValueError("npoints must be at least 2")
    if nb < 1 or kb < 1 or nb + 3 - kb < 2:
        raise ValueError("deformation spline requires nb >= 1, kb >= 1 and "
                         "nb >= kb - 1")

    guard = CallbackGuard()
    if vectorized not in (None, True, False):
        raise ValueError("vectorized must be None, True or False")
    v_cb, grad_cb, grad_many_cb = nd_callbacks(V, dV, guard, vectorized)
    handle = ctypes.c_void_p()
    status = ctypes.c_int()
    _lib.full_tunneling(
        pts.shape[0], pts.shape[1], pts, v_cb, grad_cb, grad_many_cb,
        maxiter, int(bool(verbose)),
        int(bool(ddp.get("verbose", verbose))),
        100 if V_spline_samples is None else int(V_spline_samples),
        1 if V_spline_samples is None else 0,
        alpha, npoints, nb, kb,
        int(bool(dip.get("fix_start", False))),
        int(bool(dip.get("fix_end", False))),
        float(tfp.get("xtol", 1e-4)), float(tfp.get("phitol", 1e-4)),
        ctypes.byref(handle), ctypes.byref(status))
    handle = handle.value
    try:
        guard.check()
        raise_for_status(status.value, "fullTunneling")
        nprof = ctypes.c_int()
        ndim = ctypes.c_int()
        _lib.ft_result_sizes(handle, ctypes.byref(nprof), ctypes.byref(ndim))
        action = ctypes.c_double()
        action_pot = ctypes.c_double()
        action_kin = ctypes.c_double()
        fratio = ctypes.c_double()
        rerr = ctypes.c_double()
        num_iters = ctypes.c_int()
        converged = ctypes.c_int()
        _lib.ft_result_scalars(
            handle, ctypes.byref(action), ctypes.byref(action_pot),
            ctypes.byref(action_kin), ctypes.byref(fratio),
            ctypes.byref(rerr), ctypes.byref(num_iters),
            ctypes.byref(converged))
        r = np.empty(nprof.value)
        phi = np.empty(nprof.value)
        dphi = np.empty(nprof.value)
        phi_path = np.empty((nprof.value, ndim.value), order="F")
        _lib.ft_result_arrays(handle, nprof.value, ndim.value, r, phi, dphi,
                              phi_path)
    finally:
        if handle is not None:
            _lib.ft_result_free(handle)

    profile = Profile1D(r, phi, dphi, rerr.value if rerr.value >= 0 else None)
    return fullTunneling_rval(
        profile, np.ascontiguousarray(phi_path), action.value, fratio.value,
        [], action_pot.value, action_kin.value, num_iters.value,
        bool(converged.value))
