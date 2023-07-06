# Visualization of frequency-response characteristics to Compare with verified results

import matplotlib.pyplot as plt

from matplotlib import rcParams
import numpy as np
from numpy import sqrt, pi



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

with open("DATA/RespondingRingRealModel.txt", "r") as res:
    XY = res.read().split("\n")
    Omega_real_mod = np.array([float(xy.split()[0]) for xy in XY])
    Impedance_real_mod = np.array([float(xy.split()[1]) for xy in XY])

with open("DATA/RespondingRingImModel.txt", "r") as res:
    XY = res.read().split("\n")
    Omega_im_mod = np.array([float(xy.split()[0]) for xy in XY])
    Impedance_imag_mod = np.array([float(xy.split()[1]) for xy in XY])

with open("DATA/Responding-MRI.txt", "r") as res:
    XY = res.read().split("\n")
    XY = XY[:len(XY)-1]
    Omega = np.array([float((xy.split()[0])) for xy in XY])
    Impedance_real = np.array([float(xy.split()[1]) for xy in XY])
    Impedance_imag = np.array([float(xy.split()[2]) for xy in XY])


fig, ax = plt.subplots(figsize = (10, 6))
ax.set_xlabel("Frequency, MGz")
ax.set_ylabel("Real part impedance, Ohm")
plt.grid(True)

plt.plot(Omega/2/pi/10 ** 6, Impedance_real, label=r'Modeling', color='blue')
plt.plot(Omega_real_exp, Impedance_real_exp, label = "Experiment", color = "black", linestyle = "dotted")
plt.plot(Omega_real_mod, Impedance_real_mod, label = "Verified modeling", color = "green", linestyle = "-")

plt.legend()
plt.savefig(f"Plots/MRI/Responding_Impedance_Real")
plt.show()

#Making subplots and figure of Imaginary part of Impedance on responding ring

fig, ax = plt.subplots(figsize = (10, 6))
ax.set_xlabel("Frequency, MGz")
ax.set_ylabel("Imaginary part of impedance, ohm")
ax.set_xscale('log')

plt.grid(True)

plt.plot(Omega/2/pi/10 ** 6, Impedance_imag, label=r'Modeling', color='red')
plt.plot(Omega_im_exp, Impedance_imag_exp, label = "Experiment", color = "black", linestyle = "dotted")
plt.plot(Omega_im_mod, Impedance_imag_mod, label = "Verified modeling", color = "green", linestyle = "-")

plt.legend()
plt.savefig(f"Plots/MRI/Responding_Impedance_Im")
plt.show()



