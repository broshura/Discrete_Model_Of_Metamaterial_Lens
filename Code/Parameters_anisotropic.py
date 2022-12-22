# This file contains all parameters for modeling and geometry of rings
# Special for anistropic structures

from math import pi, log
from Ring_Class import Ring
import numpy as np

# To simplify calculations distance between rings
# are normalized to cell sizes thar equals 2 units
# Also middle of structure have zero coordinates


# Geometry for our system for MRI lenz - parallelogram without borders at x and z surfaces.
def Rectangle_packing(nx, ny, nz, r, a, b, orientation = "z"):
    rings = []

    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                rings.append(Ring(i*a, j*a, k*b, orientation, r))
    return np.array(rings, dtype=Ring)

def Hexagonal_packing(nx, ny, nz, r, d, b):
    pass

# Parameters for system used in modeling

mu_0 = 4 * pi * 10 ** -7               # Permeability of vacuum

Nx, Ny, Nz = 15, 12, 34                # Number of cell on each row
a1 = 4.5 * 10 ** -3                     # Length of cell
b1 = 1 * 10 ** -3
a = 9                                  # Length of cell in packing units
b = 2
Radius1 = a1/3
Radius = a/3                           # Mean radius of rings
w1 = Radius/3
w = 0.7*Radius/3                       # Width of strip

L = 4*10**-9                 # Self-inductance
C = "Infinity"                          # Capacitance
R = 926*10**-6                              # Resistance



Rings = Rectangle_packing(Nx, Ny, Nz, Radius, a, b)   # List of Rings with their coordinates

# Adding responding ring to identify resonance frequency

#Rings = np.append(Rings, Ring(0, -4, 0, "y", 2.54/1.5*5))

Number = len(Rings)                             # Number of Rings
