# This file contains all parameters for modeling and rings

from math import pi
Thickness = 0                # Thickness of flat rings
Number  = 8                  # Number of rings
Radius  = 1                  # Middle radius of flat rings
mu_0 = 4 * pi * 10 ** -7     # Magnetic-field constant
L = 1                        # Self - Inductance
C = 1                        # Self-capacitance
R = 10                       # Resistance
j = 1
omega = 1
Z_0 = R + j * omega * L + 1/(j * omega * C)
