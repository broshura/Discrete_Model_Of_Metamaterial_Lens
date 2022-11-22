from Parameters import V, R, L, C, mu_0, a1, R_coil, L_coil
from Parameters import Number as N
import numpy as np
from scipy import linalg

from matplotlib import rcParams
import matplotlib.pyplot as plt
import numpy as np
from numpy import real, imag
# Customizing plot and fonts
rcParams['font.family'] = 'Times New Roman'

# Some lists of data to plot
Omega = []
Impedance_real = []
Impedance_imag = []
def Z_coil(omega, matrix_M, matrix_I):
    return sum([1j*omega*float(matrix_M[len(matrix_M)-1][i]) * matrix_I[i]/matrix_I[len(matrix_M)-1] for i in range(len(matrix_M)-1)]),

with open("DATA/Data.txt", "r") as res:
    RES = res.read()
    M = [[k for k in x.split(" ")] for x in RES.split("\n")]
M = [[float(M[i][k])*(a1/2) ** 1 for k in range(len(M[i]))] for i in range(len(M)-1)]
for omega_0 in range(30, 90, 2):
    omega = omega_0 * 10 ** 5
    print(omega_0)
    #Zi = [[1j*omega*M[i][k] for k in range(len(M[i]))] for i in range(len(M))]
    M_0 = R/(1j*omega) + L - 1/(omega ** 2 * C)
    M_0coil = R_coil/(1j*omega) + L_coil
    #Z_0coil = R_coil + 1j*omega * L_coil
    #Z_0 = R + 1j*omega * L + 1/(1j *omega* C)
    #Zi = Zi + np.eye(len(Zi)) * Z_0
    Zi = M + np.eye(len(M)) * M_0
    Zi[len(V)-1][len(V)-1] = M_0coil
    #Zi[len(V) - 1][len(V) - 1] = Z_0coil
    V = V[:len(V)-1] + [1/(1j*omega)]
    I = linalg.solve(Zi, V)
    print(I[-1], I[0], I[-2])
    print(M_0coil * 1j*omega, Zi[len(Zi)-1][-2] * 1j*omega, Zi[len(Zi)-1][-1]*1j*omega)
    Omega.append(omega)
    Impedance_real.append(real(Z_coil(omega, M, I)))
    Impedance_imag.append(imag(Z_coil(omega, M, I)))


# Making subplots and figure
fig, ax = plt.subplots(figsize = (10, 6))
ax.set_xlabel("Frequency, Gz")
ax.set_ylabel("Real part of impedance, Ohm")
plt.grid(True)
#Plotting dots of data
plt.scatter(Omega, Impedance_real, label=r'Dots label', color='blue')
plt.show()

fig, ax = plt.subplots(figsize = (10, 6))
ax.set_xlabel("Frequency, Gz")
ax.set_ylabel("Imaginary part of impedance, Ohm")
plt.grid(True)
#Plotting dots of data
plt.scatter(Omega, Impedance_imag, label=r'Dots label', color='blue')
plt.show()
