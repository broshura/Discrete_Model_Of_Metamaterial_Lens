# Solving matrix equation using simple Gauss method

from math import pi
from scipy.sparse.linalg import gmres

from random import random
import numpy as np
from numpy import real, imag, sqrt
from scipy.linalg import solve
from Parameters_anisotropic import *
from math import ceil

# Choose set of parameters to calculating

w_min = 1
w_max = 10
n = 100
# Calculating responding impedance normalized on ring current

def Z_coil(omega, matrix_M, matrix_I):
    return sum([1j*omega*complex(matrix_M[-1][i]) * matrix_I[i]/matrix_I[-1] for i in range(len(matrix_M))])

def Mu(I_n):
    return 1 - (pi * Radius1 ** 2 * I_n/a1**2/b1)/H_0z

# Modeling current in rings for this frequency range

Omega = np.logspace(w_min, w_max, n)                           # Frequency for each solving
Impedance_real = []                                              # Real part of impedance on responding ring
Impedance_imag = []                                              # Imaginary part of impedance on responding ring

mu = np.zeros(n, dtype=complex)
# Reading data file with geometric matrix

print("Reading data..")
with open(f"DATA/Data-{name}.txt", "r") as res:
   RES = res.read()
   M = [[k for k in x.split(" ")] for x in RES.split("\n")]
M = np.array([[complex(M[i][k])*(a1/a) ** 1 * mu_0 for k in range(len(M[i]))] for i in range(len(M)-1)])
print(f"Modeling responding ring for {name} set of parameters\n Number of rings: {len(M)}")

N_middle = ceil(Nz/2) * Nx * Ny + ceil(Ny/2)*Nx + ceil(Nx/2)
for i in range(n):
    omega = Omega[i]
    # Solving equation for matrix M instead of impedance

    print(f"Frequency: {(omega)} MGz")

    M_0 = R/1j/omega + L
    Mi = M + np.eye(Number) * M_0

    I = solve(Mi, E_w)
    print(min(I), max(I))
    mu[i] = Mu(I[N_middle])



# Add result to Data folder
print("Saving data...")
with open(f"DATA/Responding-{name}.txt", 'w') as res:
    for i in range(len(Omega)):
        res.write(f"{Omega[i]} {real(mu[i])} {imag(mu[i])}\n")
print("Ended")
