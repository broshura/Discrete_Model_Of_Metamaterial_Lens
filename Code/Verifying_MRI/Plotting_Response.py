
import matplotlib.pyplot as plt

from matplotlib import rcParams
import numpy as np
from numpy import sqrt, pi
from Parameters_anisotropic import *



with open(f"DATA/Responding-{name}.txt", "r") as res:
    XY = res.read().split("\n")
    XY = XY[:len(XY)-1]
    Omega = np.array([float((xy.split()[0])) for xy in XY])
    mu_real = np.array([float(xy.split()[1]) for xy in XY])
    mu_imag = np.array([float(xy.split()[2]) for xy in XY])

# with open(f"DATA/Responding-{name}-fast.txt", "r") as res:
#     XY = res.read().split("\n")
#     XY = XY[:len(XY)-1]
#     Omega_f = np.array([float((xy.split()[0])) for xy in XY])
#     Impedance_real_f = np.array([float(xy.split()[1]) for xy in XY])
#     Impedance_imag_f = np.array([float(xy.split()[2]) for xy in XY])
n = len(Omega)
mu_min = list(mu_imag).index(min(mu_imag))
Omega_id = f'$\omega = {round(Omega[mu_min]/10**6/3, 2)}-{round(Omega[mu_min]/10**6*3, 2)}$МГц'

fig, ax = plt.subplots(figsize = (10, 6))
ax.set_xlabel("Частота, Гц", fontsize= 16)
ax.set_ylabel("Магнитная проницаемость", fontsize = 16)
plt.xscale('log')
plt.title(f'Зависимость $\mu$ от частоты \n', fontsize = 22)
plt.plot(Omega, [0.68 for omega in Omega], linestyle = '--', color = 'blue')
#plt.grid(True)
plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)


plt.plot(Omega, mu_real, label=r'Реальная часть', color='blue', linewidth = 2)
plt.plot(Omega, mu_imag, label=r'Мнимая часть', color='red', linewidth = 2)

plt.plot([Omega[mu_min]/3 for omega in Omega],np.linspace(min(mu_imag), max(mu_real), 100), color = 'red', linestyle = '--' )
plt.plot([Omega[mu_min]*3 for omega in Omega],np.linspace(min(mu_imag), max(mu_real), 100), color = 'red', linestyle = '--' )

plt.text(10**1, 0.74, f'$\mu = 0.68$',color = 'blue', fontsize = 20)
plt.text(10 ** 1, 0.4,
         Omega_id,
         color = 'red',
         fontsize = 20
         )
#plt.xlim(10**4, 10**8)



plt.legend(fontsize = 16)
plt.savefig(f"Plots/MRI/Responding_Re(z)_{name}", dpi = 300 )
plt.show()

#Making subplots and figure of Imaginary part of Impedance on responding ring


plt.show()
