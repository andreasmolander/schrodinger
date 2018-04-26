# -*- coding: utf-8 -*-

import numpy as np


class SchrodingerEquation1D(object):
    
    def __init__(self, x, v, dx, dt, psi0, m=1.0, hbar=1.0, t0=0):
        if any(len(lst) != len(x) for lst in [v, psi0]):
            raise ValueError("Creating SchrodingerEquaion1D object " 
            + "failed, different lengths of x, v or psi0 given.")
        self.x = x
        self.v = v
        self.dx = dx
        self.dt = dt
        self.m = m
        self.psi = psi0
        self.hbar = hbar
        self.t = t0

    def prob(self):
        return self.psi * np.conjugate(self.psi)
