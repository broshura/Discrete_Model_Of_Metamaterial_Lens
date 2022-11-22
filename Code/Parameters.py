# This file contains all parameters for modeling and rings

from math import pi
from Ring_Class import Ring
import numpy as np
import matplotlib.pyplot as plt

def Rectangle_packing(nx, ny, nz, r,  type = "closed"):
    rings = []
    if type == "closed":
        for i in range(nx):
            for j in range(ny):
                for k in range(1, nz):
                    rings.append(Ring(i * 2 + 1 - nx, j * 2 + 1 - ny, k * 2 - nz, "xy", r))
        for i in range(nx):
            for j in range(ny + 1):
                for k in range(nz):
                    rings.append(Ring(i * 2 + 1 - nx, j * 2 - ny, k * 2 + 1 - nz, "xz", r))
        for i in range(1, nx):
            for j in range(ny):
                for k in range(nz):
                    rings.append(Ring(i * 2 - nx, j * 2 + 1 - ny, k * 2 + 1 - nz, "zy", r))
    elif type == "open":
        for i in range(nx):
            for j in range(ny):
                for k in range(1, nz):
                    rings.append(Ring(i * 2 + 1, j * 2 + 1, k * 2, "xy"))
        for i in range(nx):
            for j in range(1, ny):
                for k in range(nz):
                    rings.append(Ring(i * 2 + 1, j * 2, k * 2 + 1, "xz"))
        for i in range(1,nx):
            for j in range(ny):
                for k in range(nz):
                    rings.append(Ring(i * 2, j * 2 + 1, k * 2 + 1, "zy"))
        return np.array(rings, dtype=Ring)
    else:
        return "Unidentified type of packing"
    return np.array(rings, dtype=Ring)
def Sphere_Packing(Nr, type = "closed"):
    if type == "closed":
        pass
    else:
        pass




L = 13.5 * 10 ** -9                     # Self-inductance
C = 470 * 10 ** -12                     # Capacitance
R = 0.0465                              # Resistance
omega_0 = 63.28 * 10 ** 6               # Frequency of resonance in free space


Nx, Ny, Nz = 18, 2, 18                    # Number of cell on each row
a1 = 15 * 10 ** -3                       # Length of cell
a = 2
Radius = a / 3                    # Mean radius of rings
w = 0.7*0.15 * a                        # Width of strip
mu_0 = 4 * pi * 10 ** -7               # Permeability of vacuum

R_coil = 1.5
L_coil = 1.8*10**-7
#L_coil = 2.6*10**-7

Rings = Rectangle_packing(Nx, Ny, Nz, Radius)   # List of Rings with their coordinates
#Rings = [Ring(0, 0, 0, "xy", Radius), Ring(0, 0, a, "xy", Radius)]#, Ring(0, 0, 2*a, "xy", Radius)]
Number  = len(Rings)+1                     # Number of Rings
Rings = np.append(Rings, Ring(0, -4, 0, "xz", 7.62/1.5/3*5))
V = [0 for x in range(Number-1)] + [1]      # Voltage on each ring


#fig = plt.figure(figsize = (20, 20))

#ax = fig.add_subplot(1, 1, 1, projection = '3d')

#surf = ax.scatter(
#[ring.x for ring in Rings if ring.pos=="xy"], [ring.y for ring in Rings if ring.pos=="xy"], [ring.z for ring in Rings if ring.pos=="xy"],
#color = "b")
#plt.title("Centers of all xy oriented rings", size = 36)
#ax.set_ylim(-18, 18)
#ax.set_xlabel("x", size = 36)
#ax.set_ylabel("y", size = 36)
#ax.set_zlabel("z", size = 36)
#plt.show()
#print(len(Rings))
