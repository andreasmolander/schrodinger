#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib import rc

import se_functions as se_functions
import se_equations as se_equations
import se_solvers as se_solvers


def simulate():
    # Set up parameters for the simulation
    xlim = 10                               # The x-dimension goes from -xlim to xlim
    n = 500                                 # The x-dimension is discretized in n parts
    x = np.linspace(-xlim, xlim, n)         # The x-dimension
    dx = x[1] - x[0]                        # The size of a step in the x-dimension
    dt = 0.05                               # The timestep
    m = 1                                   # The mass
    hbar = 1                                # The reduced Plank's constant
    t0 = 0                                  # The initial time value

    # The harmonic potential
    k = 0.5
    v = se_functions.harmonic_v(x, k)

    # # The barrier potential
    # potential = 0
    # barrier = 50
    # v = se_functions.square_v(x, barrier, 0.5, 0.01, potential)

    # # The well potential
    # potential = 0
    # barrier = -10
    # v = se_functions.square_v(x, barrier, 0.5, 0.01, potential)

    # The wave packet for the harmonic potential
    sigma = 0.5
    x0 = -7
    k0 = 0
    psi0 = se_functions.gauss_wave_init(x, sigma, x0, k0)

    # # The wave packet for barrier/well potentials
    # sigma = 1
    # x0 = -7
    # k0 = 7
    # psi0 = se_functions.gauss_wave_init(x, sigma, x0, k0)

    # The Schrödinger equation to be solved
    se = se_equations.SchrodingerEquation1D(x, v, dx, dt, psi0, m, hbar, t0)

    # The solver used for the simulation
    se_solver = se_solvers.SeSolver1DCN(se)

    # The harmonic potential to be plotted
    vplot = se.v / 30 - 0.7

    # # The barrier potential to be plotted
    # vplot = se.v / 50

    # # The well potential to be plotted
    # vplot = se.v / 10

    # Set up the plot
    rc('text', usetex=True)
    fig = plt.figure()
    ylim = np.max([np.max(se.prob()), np.max(np.abs(vplot))]) * 1.1
    ax = fig.add_subplot(111, xlim=(-xlim, xlim), ylim=(-ylim, ylim))

    potential_plot_data, = ax.plot(se.x, vplot, 'grey', label='$V(x)$')
    probability_plot_data, = plt.plot([], [], 'xkcd:pastel orange', label=r'$\Psi^* \Psi $')
    psi_plot_data, = ax.plot([], [], 'xkcd:medium blue', label=r'$\Psi(x, t)$')
    
    plt.legend(handles=[potential_plot_data, probability_plot_data, psi_plot_data], loc=1, prop={'size': 20})

    ax.set_xlabel("$x$", fontsize=24)
    ax.set_ylabel(r"$\Psi$", fontsize=24)
    ax.set_xticklabels([])
    ax.set_yticklabels([])

    plt.tight_layout()
    plt.subplots_adjust(top=0.9)

    def init():
        psi_plot_data.set_data(se.x, se.psi)
        probability_plot_data.set_data(se.x, se.prob())
        return psi_plot_data, probability_plot_data


    def update(frame):
        se_solver.iterate()
        psi_plot_data.set_data(se.x, se.psi)
        probability_plot_data.set_data(se.x, se.prob())
        return psi_plot_data, probability_plot_data

    ani = FuncAnimation(fig, update, init_func=init, interval=1, blit=False)
    plt.show()


if __name__ == "__main__":
    simulate()
