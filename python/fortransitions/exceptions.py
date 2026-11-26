"""Exceptions, mirroring those of cosmoTransitions.

The Fortran code reports errors through integer status codes (defined in
fortransitions__config); they are translated into these exceptions.
"""


class FortransitionsError(Exception):
    """Base class of all errors raised by fortransitions.

    Raised directly for generic numerical failures (root finding, linear
    solve, spline setup) and failed spline evaluations, which have no
    counterpart in cosmoTransitions.
    """


class PotentialError(FortransitionsError):
    """The potential does not have the expected characteristics.

    As in cosmoTransitions, the second item of ``args`` is one of
    ``"no barrier"`` or ``"stable, not metastable"``.
    """


class IntegrationError(FortransitionsError):
    """The integration of the bubble profile failed."""


class DeformationError(FortransitionsError):
    """The path deformation failed."""


STATUS_OK = 0

_STATUS = {
    1: (IntegrationError, ("Integration of the bubble profile failed.",)),
    2: (PotentialError, ("No barrier between the minima; tunneling cannot "
                         "occur.", "no barrier")),
    3: (PotentialError, ("V(phi_metaMin) <= V(phi_absMin); tunneling cannot "
                         "occur.", "stable, not metastable")),
    4: (DeformationError, ("Path deformation failed.",)),
    5: (FortransitionsError, ("Numerical failure (root finding, linear solve "
                              "or spline setup).",)),
    6: (FortransitionsError, ("A spline evaluation failed (e.g. extreme "
                              "extrapolation of the path spline).",)),
}


def raise_for_status(status, where):
    """Raise the exception corresponding to a nonzero Fortran status code."""
    if status == STATUS_OK:
        return
    cls, args = _STATUS.get(
        status, (FortransitionsError, ("Unknown status code.",)))
    args = ("{}: {}".format(where, args[0]),) + args[1:]
    raise cls(*args)
