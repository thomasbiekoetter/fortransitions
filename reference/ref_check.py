"""Reference values for test/check.f90.

Runs the Python cosmoTransitions package on the two single-field potentials
from the tunneling1D docstring and on the two-field potential from
examples/fullTunneling.py, and prints the values that test/check.f90
asserts against (actions, barrier location, radial scale, profile
endpoints).

Usage: python3 ref_check.py

Requires the cosmoTransitions package to be importable (e.g. installed
with pip).
"""

import numpy as np

from cosmoTransitions import pathDeformation as pd
from cosmoTransitions import tunneling1D


# --- 1D potentials (from the tunneling1D docstring) ---

def V1_thin(phi):
    return 0.25*phi**4 - 0.49*phi**3 + 0.235*phi**2


def dV1_thin(phi):
    return phi*(phi - 0.47)*(phi - 1.0)


def V2_thick(phi):
    return 0.25*phi**4 - 0.4*phi**3 + 0.1*phi**2


def dV2_thick(phi):
    return phi*(phi - 0.2)*(phi - 1.0)


for name, V, dV in [("1d thin", V1_thin, dV1_thin),
                    ("1d thick", V2_thick, dV2_thick)]:
    inst = tunneling1D.SingleFieldInstanton(1.0, 0.0, V, dV)
    profile = inst.findProfile()
    action = inst.findAction(profile)
    print("=== %s ===" % name)
    print("phi_bar = %.8e" % inst.phi_bar)
    print("rscale  = %.8e" % inst.rscale)
    print("action  = %.10e" % action)
    print("phi(0)  = %.8e" % profile.Phi[0])


# --- 2D potential (from examples/fullTunneling.py) ---

class Potential:
    def __init__(self, c=5., fx=10., fy=10.):
        self.params = c, fx, fy

    def V(self, X):
        x, y = X[..., 0], X[..., 1]
        c, fx, fy = self.params
        r1 = x*x + c*y*y
        r2 = c*(x-1)**2 + (y-1)**2
        r3 = fx*(0.25*x**4 - x**3/3.)
        r3 += fy*(0.25*y**4 - y**3/3.)
        return r1*r2 + r3

    def dV(self, X):
        x, y = X[..., 0], X[..., 1]
        c, fx, fy = self.params
        r1 = x*x + c*y*y
        r2 = c*(x-1)**2 + (y-1)**2
        dr1dx = 2*x
        dr1dy = 2*c*y
        dr2dx = 2*c*(x-1)
        dr2dy = 2*(y-1)
        dVdx = r1*dr2dx + dr1dx*r2 + fx*x*x*(x-1)
        dVdy = r1*dr2dy + dr1dy*r2 + fy*y*y*(y-1)
        rval = np.empty_like(X)
        rval[..., 0] = dVdx
        rval[..., 1] = dVdy
        return rval


for name, fy in [("2d thin", 2.), ("2d thick", 80.)]:
    m = Potential(c=5., fx=0., fy=fy)
    Y = pd.fullTunneling([[1., 1.], [0., 0.]], m.V, m.dV, verbose=False)
    p = Y.profile1D
    print("=== %s ===" % name)
    print("action  = %.10e" % Y.action)
    print("fRatio  = %.6e" % Y.fRatio)
    print("npoints = %d" % len(p.R))
    print("Phi2D[0]  = (%.8e, %.8e)" % tuple(Y.Phi[0]))
    print("Phi2D[-1] = (%.8e, %.8e)" % tuple(Y.Phi[-1]))
