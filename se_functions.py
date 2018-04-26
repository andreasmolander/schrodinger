# -*- coding: utf-8 -*-

import numpy as np


def const_v(x, potential=0):
    v = np.full(len(x), potential)
    return v


def harmonic_v(x, k):
    return 0.5 * k * x**2


def square_v(x, square_potential, square_start, square_width, potential=0):
    v = np.full(len(x), potential)
    start = int(v.size * square_start)
    width = int(v.size * square_width)
    v[start:(start + width)] = square_potential
    return v


def gauss_wave_init(x, sigma, x0, k0):
    psi = (((2 * np.pi * sigma ** 2) ** (-0.25))
          * np.exp(1j * k0 * x - (((x - x0)/(2 * sigma)) ** 2)))
    return psi


def box_wave(x, t, edges_at=0.25, m=1.0, hbar=1.0, n=1):
    if not isinstance(n, int):
        raise ValueError("The parameter n has to be of type int.")

    lowerlimit = int(x.size * edges_at)
    upperlimit = x.size - lowerlimit

    pi = np.pi
    l = x[upperlimit] - x[lowerlimit]
    a = np.sqrt(2/l)                            # Normalization constant
    e = (n**2 * pi**2 * hbar**2)/(2*m*l**2)     # Energy
    fi = np.exp(-1j*t*e/hbar)                   # Time dependence

    if (n%2 == 0):
        psi = fi * a * np.sin((n*pi*x)/l)
    else:
        psi = fi * a * np.cos((n*pi*x)/l)

    psi[:lowerlimit] = 0
    psi[upperlimit:] = 0

    return psi
