#Calculating permeability of vacuum for structure

from math import pi
from Parameters_anisotropic import R, L, mu_0, a1, b1, Radius1, name
import matplotlib.pyplot as plt
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

# Using data for calculation

#with open(f"DATA/SumM-{name}.txt", "r") as res:
#    SumM = float(res.read())
SumM = 8 * 10 ** -6
MuReal = []
MuImag = []
for omega in Omega:
    MuReal.append(real(mu(R, L, SumM, omega, Radius1, a1, b1)))
    MuImag.append(imag(mu(R, L, SumM, omega, Radius1, a1, b1)))

# Plots of magnetic permeability of system

with open(f"DATA/SumM-{name}.txt", "r") as res:
   SumM = int(res.read())

fig, ax = plt.subplots(figsize = (10, 6))
ax.set_xlabel("Frequency, Gz")
ax.set_xscale("log")
ax.set_ylabel("Real part of magnetic permeability")
plt.grid(True)

plt.plot(Omega, MuReal, label=r'Real part', color='blue')
plt.savefig(f"Plots/Anisotropic-mureal")
plt.show()


fig, ax = plt.subplots(figsize = (10, 6))
ax.set_xlabel("Frequency, Gz")
ax.set_xscale("log")
ax.set_ylabel("Imaginary part of magnetic permeability")
plt.grid(True)

plt.scatter(Omega, MuImag, label=r'Imaginary part', color='blue')
plt.savefig(f"Plots/Anisotropic-muim")
plt.show()

