# Visualization of each part of modeling

import matplotlib.pyplot as plt

from matplotlib import rcParams
import numpy as np
from numpy import sqrt
from Permeability_anisotropic import Omega, MuReal, MuImag

# Function for drawing circles in 3D

def plot_circle(ring):
    n = 250
    #if ring.r > 1:
    #    Y = np.array([ring.y for i in range(n)])
    #    Z = np.linspace(ring.z - ring.r, ring.z + ring.r, n // 2)
    #    X = ring.x + sqrt(ring.r ** 2 - (Z - ring.z) ** 2 + 0.001)
    #    X = np.append(X, 2 * ring.x - X)
    #    Z = np.append(Z, Z[::-1])
    #    ax.plot(
    #        X, Y, Z, color="green", linewidth = 1
    #    )
    if ring.pos == "x":
        X = np.array([ring.x for i in range(n)])
        Y = np.linspace(ring.y - ring.r, ring.y + ring.r, n//2)
        Z = ring.z + sqrt(ring.r ** 2 - (Y-ring.y) ** 2+0.001)
        Z = np.append(Z, 2*ring.z-Z)
        Y = np.append(Y, Y[::-1])
        ax.plot(
            X, Y, Z, color = "blue", linewidth =1
        )
    elif ring.pos == "y":
        Y = np.array([ring.y for i in range(n)])
        Z = np.linspace(ring.z - ring.r, ring.z + ring.r, n//2)
        X = ring.x + sqrt(ring.r ** 2 - (Z-ring.z) ** 2+0.001)
        X = np.append(X, 2*ring.x-X)
        Z = np.append(Z, Z[::-1])
        ax.plot(
            X, Y, Z, color = "red",  linewidth =1
        )
    elif ring.pos == "z" and (ring.x == 9*14 or ring.y == 0 or ring.z >= 2*32):
        Z = np.array([ring.z for i in range(n)])
        X = np.linspace(ring.x - ring.r, ring.x + ring.r, n//2)
        Y = ring.y + sqrt(ring.r ** 2 - (X - ring.x) ** 2)
        Y = np.append(Y, (2 * ring.y - Y))
        X = np.append(X, X[::-1])
        ax.plot(
            X, Y, Z, color="black",  linewidth =1
        )

# Customizing plot and fonts
rcParams['font.family'] = 'Times New Roman'

# Repeating plots for MRI lenz of 3 x 18 x 18 lenz

#from Simple import Impedance_real, Impedance_imag, Omega

# Making subplots and figure for Real part of Impedance on responding ring

from Simple import Omega, Impedance_imag, Impedance_real

fig, ax = plt.subplots(figsize = (10, 6))
ax.set_xlabel("Frequency, Gz")
ax.set_ylabel("Real part impedance, Ohm")
plt.grid(True)

plt.plot(Omega, Impedance_real, label=r'Real part', color='blue')
plt.savefig(f"Plots/Responding_Impedance_Real")
plt.show()

#Making subplots and figure of Imaginary part of Impedance on responding ring

fig, ax = plt.subplots(figsize = (10, 6))
ax.set_xlabel("Frequency, Gz")
ax.set_ylabel("Imaginary part of impedance, ohm")
plt.grid(True)

plt.plot(Omega, Impedance_imag, label=r'Imaginary part', color='red')
plt.savefig(f"Plots/Responding_Impedance_Im")
plt.show()

# Plots of magnetic permeability of system

# Cj

#from Permeability_anistotropic import Omega, MuReal, MuImag, name

#with open(f"DATA/SumM-{name}.txt", "r") as res:
#    SumM = int(res.read())
#
# fig, ax = plt.subplots(figsize = (10, 6))
# ax.set_xlabel("Frequency, Gz")
# ax.set_xscale("log")
# ax.set_ylabel("Real part of magnetic permeability")
# plt.grid(True)
#
# plt.plot(Omega, MuReal, label=r'Real part', color='blue')
# plt.savefig(f"Plots/Anisotropic-mureal")
# plt.show()
#
#
# fig, ax = plt.subplots(figsize = (10, 6))
# ax.set_xlabel("Frequency, Gz")
# ax.set_xscale("log")
# ax.set_ylabel("Imaginary part of magnetic permeability")
# plt.grid(True)
#
# plt.scatter(Omega, MuImag, label=r'Imaginary part', color='blue')
# plt.savefig(f"Plots/Anisotropic-muim")
# plt.show()

# 3D plot of structure with dots in middles of each ring

# Choosing set of parameters

##from Parameters import Rings

##from Parameters_anisotropic import Rings

# fig = plt.figure(figsize = (20, 20))
# ax = fig.add_subplot(1, 1, 1, projection = '3d')
#
# Orientation = "z"                                               # Orientation of ring that will be shown
#
# for ring in Rings[:len(Rings)]:
#     if ring.pos == Orientation:
#         plot_circle(ring)
#
# # Other parameters of plot
# plt.title(f"Anisotropic structure", size = 36) #Size
# ax.set_ylim(0, 9*12)
# ax.set_xlabel("x", size = 36)
# ax.set_ylabel("y", size = 36)
# ax.set_zlabel("z", size = 36)
#
# plt.savefig(f"Plots/Anisotropic2-rings")
# plt.show()

