# This file contains all parameters for modeling and rings

from math import pi
from Ring_Class import Ring
import numpy as np

def Rectangle_packing(d, nx, ny, nz, type = "closed"):
    if type == "closed":
        rings = []
        for i in range(nz + 1):
            for j in range(ny):
                for k in range(nz):
                    rings.append(Ring(i*d - 2, j*d - 2, k*d - 2, "xy"))
        for i in range(ny + 1):
            for j in range(ny):
                for k in range(nz):
                    rings.append(Ring(i*d - 2, j*d - d/2 - 2, k*d + d/2 - 2, "xz"))
        for i in range(nx + 1):
            for j in range(ny):
                for k in range(nz):
                    rings.append(Ring(i*d - d/2 - 2, j*d - 2, k*d+d/2 - 2, "zy"))
    else:
        pass
    return np.array(rings, dtype=Ring)
def Sphere_Packing(Nr, type = "closed"):
    if type == "closed":
        pass
    else:
        pass

Thickness = 0
Nx = 1
Ny = 1
Nz = 1
Number  = (Nx + 1) * Ny * Nz + Nx * (Ny + 1) * Nz + Nx * Ny * (Nz + 1)
Radius  = 1
Error = 1.1
mu_0 = 4 * pi * 10 ** -7
L = 1
C = 1
R = 10
j = 1
omega = 1
Z_0 = R + j * omega * L + 1/(j * omega * C)

V = np.array([1 for x in range(Number)])
Rings = Rectangle_packing(2 * Radius * Error, Nx, Ny, Nz)

