from Parameters import V, R, L, C, mu_0, a
from Parameters import Number as N
import numpy as np
from scipy import linalg

from matplotlib import rcParams
import matplotlib.pyplot as plt
import numpy as np
from numpy import real
# Customizing plot and fonts
rcParams['font.family'] = 'Times New Roman'

# Some lists of data to plot
x = []
y = []

def Z_coil(omega, matrix_M, matrix_I):
    return sum([1j*omega*float(matrix_M[len(matrix_M)-1][i]) * matrix_I[i]/matrix_I[len(matrix_M)-1] for i in range(len(matrix_M)-1)]),

with open("DATA/Data.txt", "r") as res:
    RES = res.read()
    M = [[k for k in x.split(" ")] for x in RES.split("\n")]
M = [[mu_0*float(M[i][k])*(a/2) ** 1 for k in range(len(M[i]))] for i in range(len(M)-1)]
for omega_0 in range(30, 90, 2):
    omega = omega_0 * 10 ** 6
    print(omega_0, omega/10**9)
    #Zi = [[1j*omega*M[i][k] for k in range(len(M[i]))] for i in range(len(M))]
    Z_0 = R/(1j*omega) + L - 1/(omega ** 2 * C)
    #Z_0 = R + 1j*omega * L + 1/(1j *omega* C)
    #Zi = Zi + np.eye(len(Zi)) * Z_0
    Zi = M + np.eye(len(M)) * Z_0
    V = V[:len(V)-1] + [1/(1j*omega)]
    I = linalg.solve(Zi, V)
    x.append(omega)
    y.append(real(Z_coil(omega, M, I)))


# Making subplots and figure
fig, ax = plt.subplots(figsize = (10, 6))
ax.set_title('Plot title', fontsize=14)

#Plotting dots of data
plt.scatter(x, y, label=r'Dots label', color='blue')

plt.show()
