# Visualization of each part of modeling

import matplotlib.pyplot as plt
from Parameters import Rings
from matplotlib import rcParams
import numpy as np
from numpy import sqrt
from Simple import Impedance_real, Impedance_imag, Omega

# Function for drawing circles in 3D

def plot_circle(ring):
    n = 250
    if ring.r > 1:
        Y = np.array([ring.y for i in range(n)])
        Z = np.linspace(ring.z - ring.r, ring.z + ring.r, n // 2)
        X = ring.x + sqrt(ring.r ** 2 - (Z - ring.z) ** 2 + 0.001)
        X = np.append(X, 2 * ring.x - X)
        Z = np.append(Z, Z[::-1])
        ax.plot(
            X, Y, Z, color="green", linewidth = 20
        )
    elif ring.pos == "x":
        X = np.array([ring.x for i in range(n)])
        Y = np.linspace(ring.y - ring.r, ring.y + ring.r, n//2)
        Z = ring.z + sqrt(ring.r ** 2 - (Y-ring.y) ** 2+0.001)
        Z = np.append(Z, 2*ring.z-Z)
        Y = np.append(Y, Y[::-1])
        ax.plot(
            X, Y, Z, color = "blue"
        )
    elif ring.pos == "y":
        Y = np.array([ring.y for i in range(n)])
        Z = np.linspace(ring.z - ring.r, ring.z + ring.r, n//2)
        X = ring.x + sqrt(ring.r ** 2 - (Z-ring.z) ** 2+0.001)
        X = np.append(X, 2*ring.x-X)
        Z = np.append(Z, Z[::-1])
        ax.plot(
            X, Y, Z, color = "red"
        )
    elif ring.pos == "z":
        Z = np.array([ring.z for i in range(n)])
        X = np.linspace(ring.x - ring.r, ring.x + ring.r, n//2)
        Y = ring.y + sqrt(ring.r ** 2 - (X - ring.x) ** 2)
        Y = np.append(Y, (2 * ring.y - Y))
        X = np.append(X, X[::-1])
        ax.plot(
            X, Y, Z, color="black"
        )

# Customizing plot and fonts
rcParams['font.family'] = 'Times New Roman'

# Repeating plots for MRI lenz of 3 x 18 x 18 lenz

# Making subplots and figure of Real part of Impedance on responding ring

fig, ax = plt.subplots(figsize = (10, 6))
ax.set_xlabel("Frequency, Gz")
ax.set_ylabel("Real part of impedance, Ohm")
plt.grid(True)

#Plotting dots of data

plt.scatter(Omega, Impedance_real, label=r'Dots label', color='blue')
plt.show()

#Making subplots and figure of Imaginary part of Impedance on responding ring

fig, ax = plt.subplots(figsize = (10, 6))
ax.set_xlabel("Frequency, Gz")
ax.set_ylabel("Imaginary part of impedance, Ohm")
plt.grid(True)

#Plotting dots of data
plt.scatter(Omega, Impedance_imag, label=r'Dots label', color='blue')
plt.show()




# 3D plot of structure with dots in middles of each ring

fig = plt.figure(figsize = (20, 20))

ax = fig.add_subplot(1, 1, 1, projection = '3d')

Orientation = "z"                                               # Orientation of ring that will be shown


for ring in Rings[:len(Rings)]:                                              # Coloring and so one
    if ring.pos == "z":
        plot_circle(ring)

# Other parameters of plot
plt.title(f"Centers of all {Orientation}-oriented rings", size = 36) #Size
ax.set_ylim(-18, 18)
ax.set_xlabel("x", size = 36)
ax.set_ylabel("y", size = 36)
ax.set_zlabel("z", size = 36)

plt.savefig("Plots/Z-rings")
plt.show()

