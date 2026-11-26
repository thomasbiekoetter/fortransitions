"""Loading of the shared library and ctypes declarations of its C interface.

The C interface is defined by the bind(C) routines in the Fortran modules
fortransitions__c_potentials, fortransitions__c_pathdeformation and
fortransitions__c_tunneling1d (../src/c_*.f90).
"""

import ctypes
import os

import numpy as np

_lib_path = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "lib", "libfortransitions.so")

try:
    lib = ctypes.CDLL(_lib_path)
except OSError as err:
    raise ImportError(
        "Could not load the fortransitions shared library at {}. Build it "
        "first with `pixi run install` or `make wrapper` in the python/ "
        "directory. ({})".format(_lib_path, err)) from err

c_int = ctypes.c_int
c_double = ctypes.c_double
c_void_p = ctypes.c_void_p
p_int = ctypes.POINTER(ctypes.c_int)
p_double = ctypes.POINTER(ctypes.c_double)
p_void_p = ctypes.POINTER(ctypes.c_void_p)

# Callback prototypes (abstract interfaces in fortransitions__c_potentials)
V_ND_CB = ctypes.CFUNCTYPE(c_double, c_int, p_double)
GRAD_ND_CB = ctypes.CFUNCTYPE(None, c_int, p_double, p_double)
F_1D_CB = ctypes.CFUNCTYPE(c_double, c_double)


def _arr(ndim=1, order="C"):
    flags = "F_CONTIGUOUS" if order == "F" else "C_CONTIGUOUS"
    return np.ctypeslib.ndpointer(dtype=np.float64, ndim=ndim, flags=flags)


def _declare(name, argtypes):
    func = getattr(lib, name)
    func.argtypes = argtypes
    func.restype = None
    return func


# full_tunneling
full_tunneling = _declare("fortransitions_full_tunneling", [
    c_int, c_int, _arr(2, "F"),           # npts, ndim, path_pts
    V_ND_CB, GRAD_ND_CB,                  # v_cb, grad_cb
    c_int, c_int, c_int,                  # maxiter, verbose, deform_verbose
    c_int, c_int,                         # v_spline_samples, v_mode
    c_double, c_int, c_int, c_int,        # alpha, npoints, nb, kb
    c_int, c_int,                         # fix_start, fix_end
    c_double, c_double,                   # xtol, phitol
    p_void_p, p_int,                      # handle, status
])
ft_result_sizes = _declare("fortransitions_ft_result_sizes", [
    c_void_p, p_int, p_int])
ft_result_scalars = _declare("fortransitions_ft_result_scalars", [
    c_void_p, p_double, p_double, p_double, p_double, p_double, p_int,
    p_int])
ft_result_arrays = _declare("fortransitions_ft_result_arrays", [
    c_void_p, c_int, c_int, _arr(), _arr(), _arr(), _arr(2, "F")])
ft_result_free = _declare("fortransitions_ft_result_free", [c_void_p])

# single_field_instanton
sfi_new = _declare("fortransitions_sfi_new", [
    c_double, c_double,                   # phi_absmin, phi_metamin
    F_1D_CB, F_1D_CB, F_1D_CB,            # v_cb, dv_cb, d2v_cb
    c_double, c_double,                   # alpha, phi_eps_rel
    c_int, c_double, c_int, c_double,     # has_phi_bar, phi_bar, has_rscale, rscale
    p_void_p, p_int,                      # handle, status
])
sfi_params = _declare("fortransitions_sfi_params", [
    c_void_p, p_double, p_double])
sfi_find_profile = _declare("fortransitions_sfi_find_profile", [
    c_void_p,
    c_int, c_double,                      # has_xguess, xguess
    c_double, c_double, c_double,         # xtol, phitol, thin_cutoff
    c_int, c_double, c_double,            # npoints, rmin_rel, rmax_rel
    c_int, c_int,                         # has_max_interior, max_interior_pts
    p_void_p, p_int,                      # prof_handle, status
])
sfi_find_action = _declare("fortransitions_sfi_find_action", [
    c_void_p, c_int, _arr(), _arr(), _arr(), c_double,
    p_double, p_double, p_double])
sfi_free = _declare("fortransitions_sfi_free", [c_void_p])
profile_size = _declare("fortransitions_profile_size", [c_void_p, p_int])
profile_arrays = _declare("fortransitions_profile_arrays", [
    c_void_p, c_int, _arr(), _arr(), _arr(), p_double])
profile_free = _declare("fortransitions_profile_free", [c_void_p])
