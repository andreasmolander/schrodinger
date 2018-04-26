# -*- coding: utf-8 -*-

import numpy as np
from scipy.linalg import solve


class SeSolver(object):

    def __init__(self, se):
        self.se = se

    def iterate(self):
        raise NotImplementedError


class SeSolver1DCN(SeSolver):

    def __init__(self, se):
        super(SeSolver1DCN, self).__init__(se)
        self.__cn_current = self.__init_cn_matrix_current()
        self.__cn_next = self.__init_cn_matrix_next()

    def __init_cn_matrix_current(self):

        def __get_matrix_element(i, j):
            if i == j:
                return (1 - ((1j * self.se.dt) / (2 * self.se.hbar)
                        * ((self.se.hbar ** 2) / (self.se.m * self.se.dx ** 2)
                        + self.se.v[i])))
            elif abs(i-j) == 1:
                return (1 + (1j * self.se.dt) / (2 * self.se.hbar)
                        * (self.se.hbar ** 2)
                        / (2 * self.se.m * self.se.dx ** 2))
            else:
                return 0

        cn_matrix = np.zeros([len(self.se.x), len(self.se.x)], dtype=complex)

        for i in range(len(self.se.x)):
            for j in range(len(self.se.x)):
                cn_matrix[i][j] = __get_matrix_element(i, j)

        return cn_matrix

    def __init_cn_matrix_next(self):
        
        def __get_matrix_element(i, j):
            if i == j:
                return (1 + ((1j * self.se.dt) / (2 * self.se.hbar)
                        * ((self.se.hbar ** 2) / (self.se.m * self.se.dx ** 2) 
                        + self.se.v[i])))
            elif abs(i-j) == 1:
                return (1 - (1j * self.se.dt) / (2 * self.se.hbar)
                        * (self.se.hbar ** 2)
                        / (2 * self.se.m * self.se.dx ** 2))
            else:
                return 0

        cn_matrix = np.empty([len(self.se.x), len(self.se.x)], dtype=complex)

        for i in range(len(self.se.x)):
            for j in range(len(self.se.x)):
                cn_matrix[i][j] = __get_matrix_element(i, j)

        return cn_matrix

    def crank_nicholson(self):
        a = np.dot(self.__cn_current, self.se.psi)
        self.se.psi = solve(self.__cn_next, a)

    def iterate(self):
        self.crank_nicholson()
        self.se.t = self.se.t + self.se.dt


class SeSolver1DFE(SeSolver):

    def __init__(self, se):
        super(SeSolver1DFE, self).__init__(se)
        self.__h_matrix = self.__init_h_matrix()

    def __get_h_matrix_element(self, i, j):
        if i == j:
            return ((self.se.hbar ** 2 / self.se.m * self.se.dx ** 2)
                     + self.se.v[i])
        elif abs(i-j) == 1:
            return - self.se.hbar ** 2 / 2 * self.se.m * self.se.dx ** 2
        else:
            return 0

    def __init_h_matrix(self):
        h_matrix = np.empty([len(self.se.x), len(self.se.x)])
        for i in range(len(self.se.x)):
            for j in range(len(self.se.x)):
                h_matrix[i][j] = self.__get_h_matrix_element(i, j)
        return h_matrix

    def forward_euler(self):
        self.se.psi = (self.se.psi - (1j * self.se.dt) / self.se.hbar
                        * np.dot(self.__h_matrix, self.se.psi))

    def iterate(self):
        self.forward_euler()
        self.se.t = self.se.t + self.se.dt
