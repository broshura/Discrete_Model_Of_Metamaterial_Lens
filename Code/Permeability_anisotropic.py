#Calculating permeability of vacuum for structure

from math import pi
from Parameters_anisotropic import R, L, mu_0, a1, b1, Radius1
from scipy import linalg

import numpy as np
from numpy import real, imag
import numpy as np
def mu(R, L, SumM, omega, r_0, a, b):
    n_0 = 1/a**2/b
    Z_0 = R + 1j*omega*L
    Z = SumM * 1j*omega
    const = 1j*omega* pi ** 2*r_0**4*mu_0*n_0
    return (Z_0 + Z + 2/3*const)/(Z_0 + Z - 1/3 * const)


Omega = np.logspace(-1, 10, 1000)
SumM = 8.724352040907938e-06

MuReal = []
MuImag = []
for omega in Omega:
    MuReal.append(real(mu(R, L, SumM, omega, Radius1, a1, b1)))
    MuImag.append(imag(mu(R, L, SumM, omega, Radius1, a1, b1)))
print(MuReal[999])
