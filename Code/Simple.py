# Solving matrix equation using simple Gauss method

from math import pi
from scipy import linalg

import numpy as np
from numpy import real, imag, sqrt
#from numpy import linalg
# Choose set of parameters to calculating

from Parameters import R, L, C, mu_0, a1, a, R_coil, L_coil, V, name

# Calculating responding impedance normalized on ring current

def Z_coil(omega, matrix_M, matrix_I):
    return sum([1j*omega*complex(matrix_M[-1][i]) * matrix_I[i]/matrix_I[-1] for i in range(len(matrix_M))])

# Modeling current in rings for this frequency range

w_max = 75*10 ** 6 * 2 * pi     # Maximum frequency in MGz
w_min = 55*10 ** 6 * 2 * pi   # Minimum frequency in MGz

Omega = np.linspace(w_min, w_max, 100)                           # Frequency for each solving
Impedance_real = []                                              # Real part of impedance on responding ring
Impedance_imag = []                                              # Imaginary part of impedance on responding ring

# Reading data file with geometric matrix

print("Reading data..")
with open(f"DATA/Data-{name}.txt", "r") as res:
   RES = res.read()
   M = [[k for k in x.split(" ")] for x in RES.split("\n")]
M = [[float(M[i][k])*(a1/a) ** 1 * mu_0/(2*pi) for k in range(len(M[i]))] for i in range(len(M)-1)]
print(f"Modeling responding ring for {name} set of parameters\n Number of rings: {len(M)}")

for omega in Omega:
    # Solving equation for matrix M instead of impedance

    print(omega)
    M_0 = (R/(1j*omega) + L - 1/(omega ** 2 * C))    # Self-impedance of each ring
    M_0coil = R_coil/(1j*omega) + L_coil           # Self-impedance of responding ring
    Mi = M + np.eye(len(M)) * M_0                  # Add self-impedance of diagonal
    Mi[len(V)-1][len(V)-1] = M_0coil
    Vi = [v/1j/omega for v in V]                    # Correcting voltage to new matrix equation

    # Current list

    I = linalg.solve(Mi, Vi)
    print(I)
    print(M[-1])
    print(max(M[-1]))
    Impedance_real.append(real(Z_coil(omega, M, I)))
    Impedance_imag.append(imag(Z_coil(omega, M, I)))

Impedance_real = np.asarray(Impedance_real)
Impedance_imag = np.asarray(Impedance_imag)

# Add result to Data folder
print("Saving data...")
with open(f"DATA/Responding-{name}.txt", 'w') as res:
    for i in range(len(Omega)):
        res.write(f"{Omega[i]} {Impedance_real[i]} {Impedance_imag[i]}\n")
print("Ended")
