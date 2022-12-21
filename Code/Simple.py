# Solving matrix equation using simple Gauss method

from math import pi
from Parameters import V, R, L, C, mu_0, a1, R_coil, L_coil, Rings, a
from scipy import linalg

import numpy as np
from numpy import real, imag


# Some lists of data to plot

# Modeling current in rings for this frequency range

w_max = 90 * 10 ** 6    # Maximum integer frequency in MGz
w_min = 30 * 10 ** 6    # Minimum integer frequency in MGz

Omega = np.linspace(w_min, w_max, 50)                            # Frequency for each solving
Impedance_real = []                                              # Real part of impedance on responding ring
Impedance_imag = []                                              # Imaginary part of impedance on responding ring

# Calculating responding impedance normalized on ring current

def Z_coil(omega, matrix_M, matrix_I):
    return sum([1j*omega*float(matrix_M[len(matrix_M)-1][i]) * matrix_I[i]/matrix_I[len(matrix_M)-1] for i in range(len(matrix_M)-1)]),

# Reading data file with geometric matrix

with open("DATA/Data.txt", "r") as res:
    RES = res.read()
    M = [[k for k in x.split(" ")] for x in RES.split("\n")]
M = [[float(M[i][k])*(a1/a) ** 1 * mu_0/(4*pi) for k in range(len(M[i]))] for i in range(len(M)-1)]


for omega in Omega:

    # Solving equation for matrix M instead of impedance

    M_0 = (R/(1j*omega) + L - 1/(omega ** 2 * C))    # Self-impedance of each ring
    M_0coil = R_coil/(1j*omega) + L_coil           # Self-impedance of responding ring
    Mi = M + np.eye(len(M)) * M_0                  # Add self-impedance of diagonal
    Mi[len(V)-1][len(V)-1] = M_0coil
    V = V[:len(V)-1] + [1/(1j*omega)]              # Correcting voltage to new matrix equation

    # Current list
    I = linalg.solve(Mi, V)

    Impedance_real.append(real(Z_coil(omega, M, I)))
    Impedance_imag.append(imag(Z_coil(omega, M, I)))

