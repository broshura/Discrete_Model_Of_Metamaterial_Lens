#Calculating permeability of vacuum for structure

from math import pi
from Parameters_anisotropic import R, L, mu_0, a1, a, b1, b, Radius1, name
import matplotlib.pyplot as plt
import numpy as np
from numpy import real, imag
import numpy as np

def mu(R, L, SumM, omega, r_0, a, b):
    Z_0 = R + 1j*omega*L
    Sigma = SumM/r_0
    Const = a ** 3 / pi ** 2 / r_0 ** 3
    print(Z_0)
    print(Z_0/omega/mu_0/r_0)
    print(SumM)
    print(Const * (1j * Z_0/omega/mu_0/r_0 + Sigma))
    return 1 - (Const * (1j * Z_0/omega/mu_0/r_0 + Sigma) + 1 / 3) ** -1

Omega = np.logspace(-1, 10, 1000)

# Using data for calculation

with open(f"DATA/SumM-{name}.txt", "r") as res:
    SumM = complex(res.read()) * a1/a * mu_0

MuReal = []
MuImag = []
for omega in Omega:
    MuReal.append(real(mu(R, L, SumM, omega, Radius1, a1, b1)))
    MuImag.append(imag(mu(R, L, SumM, omega, Radius1, a1, b1)))
print(MuReal)
# Plots of magnetic permeability of system


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

plt.plot(Omega, MuImag, label=r'Imaginary part', color='blue')
plt.savefig(f"Plots/Anisotropic-muim")
plt.show()

