# Visualization of each part of modeling

import matplotlib.pyplot as plt

from matplotlib import rcParams
import numpy as np
from numpy import sqrt, pi
from Permeability_anisotropic import Omega, MuReal, MuImag

# Customizing plot and fonts
rcParams['font.family'] = 'Times New Roman'

# Repeating plots for MRI lenz of 3 x 18 x 18 lenz

# Making subplots and figure for Real part of Impedance on responding ring

# Using data from modeling and published plots

with open("DATA/RespondingRingReal.txt", "r") as res:
    XY = res.read().split("\n")
    Omega_real_exp = np.array([float(xy.split()[0]) for xy in XY])
    Impedance_real_exp = np.array([float(xy.split()[1]) for xy in XY])

with open("DATA/RespondingRingIm.txt", "r") as res:
    XY = res.read().split("\n")
    Omega_im_exp = np.array([float(xy.split()[0]) for xy in XY])
    Impedance_imag_exp = np.array([float(xy.split()[1]) for xy in XY])

with open("DATA/Responding-MRI.txt", "r") as res:
    XY = res.read().split("\n")
    XY = XY[:len(XY)-1]
    Omega = np.array([float((xy.split()[0])) for xy in XY])
    Impedance_real = np.array([float(xy.split()[1]) for xy in XY])
    Impedance_imag = np.array([float(xy.split()[2]) for xy in XY])


Impedance_real_exp = Impedance_real_exp/max(Impedance_real_exp)
Impedance_imag_exp = Impedance_imag_exp/max(Impedance_imag_exp)
Impedance_real = Impedance_real/max(Impedance_real)
Impedance_imag = Impedance_imag/max(Impedance_imag)

fig, ax = plt.subplots(figsize = (10, 6))
ax.set_xlabel("Frequency, MGz")
ax.set_ylabel("Real part impedance, Ohm")
plt.grid(True)

plt.plot(Omega/2/pi/10 ** 6, Impedance_real, label=r'Modeling', color='blue')
plt.plot(Omega_real_exp, Impedance_real_exp, label = "Experiment", color = "black", linestyle = "-")
plt.savefig(f"Plots/Responding_Impedance_Real")
plt.show()

#Making subplots and figure of Imaginary part of Impedance on responding ring

fig, ax = plt.subplots(figsize = (10, 6))
ax.set_xlabel("Frequency, MGz")
ax.set_ylabel("Imaginary part of impedance, ohm")
plt.grid(True)

plt.plot(Omega/2/pi/10 ** 6, Impedance_imag, label=r'Imaginary part', color='red')
plt.plot(Omega_im_exp, Impedance_imag_exp, label = "Experiment", color = "black", linestyle = "-")
plt.savefig(f"Plots/Responding_Impedance_Im")
plt.show()



