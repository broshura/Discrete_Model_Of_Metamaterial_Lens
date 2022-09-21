# This file contains all parameters for modeling and rings

from math import pi
from Ring_Class import Ring
import numpy as np

def Rectangle_packing(nx, ny, nz, type = "closed"):
    rings = []
    if type == "closed":
        for i in range(nx):
            for j in range(ny):
                for k in range(nz + 1):
                    rings.append(Ring(i * 2 + 1, j * 2 + 1, k * 2, "xy"))
        for i in range(nx):
            for j in range(ny + 1):
                for k in range(nz):
                    rings.append(Ring(i * 2 + 1, j * 2, k * 2 + 1, "xz"))
        for i in range(nx + 1):
            for j in range(ny):
                for k in range(nz):
                    rings.append(Ring(i * 2, j * 2 + 1, k * 2 + 1, "zy"))
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

Nx, Ny, Nz = 2, 2, 2                    # Number of cell on each row
a = 15 * 10 ** -3                       # Length of cell
Rings = Rectangle_packing(Nx, Ny, Nz)   # List of Rings with their coordinates
Number  = len(Rings)                    # Number of Rings
Radius  = 4.9 * 10 ** -3 / a * 2        # Mean radius of rings
w = 2.2 * 10 ** -3 / a * 2  * 0         # Width of strip
mu_0 = 4 * pi * 10 ** -7                # Permeability of vacuum
L = 13.5 * 10 ** -9                     # Self-inductance
C = 470 * 10 ** -12                     # Capacitance
R = 0.0465                              # Resistance
omega = 63.28 * 10 ** 9                 # Frequency of resonance in free space

Thickness = 0
Z_0 = R + 1j * omega * L + 1/(1j * omega * C) # Self-impedance
Z_0 = round(Z_0.real) + 1j*round(Z_0.imag)
V = np.array([1 for x in range(Number)])      # Voltage on each ring


