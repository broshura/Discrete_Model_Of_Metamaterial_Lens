# Combination of different parameters and method for modeling

# Solving and plotting frequency-response characteristics for MRI lenz

from Parameters_MRI import *

# Geometry part

from Impedance_matrix import *
Impedance_matrix(Rings, name, w)    # Making Data-files with Impedance matrix

# Part with solving matrix equation in the frequency range
from Simple_Method import *

w_max = 75 * 10 ** 6 * 2 * pi       # Maximum frequency in MGz
w_min = 55 * 10 ** 6 * 2 * pi       # Minimum frequency in MGz
n = 100                             # number of dots

Simple_Method(name, w_min, w_max, n, R, L, C, mu_0, a1, a, R_coil, L_coil, V)

#Final Plotting to compare results

import Plotting
