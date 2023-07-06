#Calculating permeability of vacuum for structure

from math import pi
from Parameters_anisotropic import R, L, mu_0, a1, a, Radius1, name, n_0, Number
import matplotlib.pyplot as plt
import numpy as np
from numpy import real, imag
import numpy as np

def mu(omega):
    Z_0 = R + 1j * omega * L
    SumZ = SumM * 1j * omega
    C = 1j * omega * pi ** 2 * Radius1 ** 4 * n_0 * mu_0
    print(Z_0, SumZ, C)
    return (Z_0 + SumZ + 2/3 * C)/(Z_0 + SumZ - 1/3 * C)



Omega = np.logspace(-1, 10, 1000)

# Using data for calculation

with open(f"DATA/SumM-{name}.txt", "r") as res:
    SumM = complex(res.read()) * a1/a * mu_0

#with open(f"DATA/Data-{name}.txt", "r") as res:
#    RES = res.read()
#    M = [[k for k in x.split(" ")] for x in RES.split("\n")]

#M = np.array([[complex(M[i][k])*(a1/a) ** 1 * mu_0 for k in range(len(M[i]))] for i in range(len(M)-1)])


MuReal = []
MuImag = []
for omega in Omega:
    MuReal.append(real(mu(omega)))
    MuImag.append(imag(mu(omega)))
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

