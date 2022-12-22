# 3D visualization of structure

import matplotlib.pyplot as plt

from matplotlib import rcParams
import numpy as np
from numpy import sqrt, pi

# Function for drawing circles in 3D

def plot_circle(ring):
    n = 250
    #if ring.r > 1:
    #    Y = np.array([ring.y for i in range(n)])
    #    Z = np.linspace(ring.z - ring.r, ring.z + ring.r, n // 2)
    #    X = ring.x + sqrt(ring.r ** 2 - (Z - ring.z) ** 2 + 0.001)
    #    X = np.append(X, 2 * ring.x - X)
    #    Z = np.append(Z, Z[::-1])
    #    ax.plot(
    #        X, Y, Z, color="green", linewidth = 1
    #    )
    if ring.pos == "x":
        X = np.array([ring.x for i in range(n)])
        Y = np.linspace(ring.y - ring.r, ring.y + ring.r, n//2)
        Z = ring.z + sqrt(ring.r ** 2 - (Y-ring.y) ** 2+0.001)
        Z = np.append(Z, 2*ring.z-Z)
        Y = np.append(Y, Y[::-1])
        ax.plot(
            X, Y, Z, color = "blue", linewidth =1
        )
    elif ring.pos == "y":
        Y = np.array([ring.y for i in range(n)])
        Z = np.linspace(ring.z - ring.r, ring.z + ring.r, n//2)
        X = ring.x + sqrt(ring.r ** 2 - (Z-ring.z) ** 2+0.001)
        X = np.append(X, 2*ring.x-X)
        Z = np.append(Z, Z[::-1])
        ax.plot(
            X, Y, Z, color = "red",  linewidth =1
        )
    elif ring.pos == "z":
        Z = np.array([ring.z for i in range(n)])
        X = np.linspace(ring.x - ring.r, ring.x + ring.r, n//2)
        Y = ring.y + sqrt(ring.r ** 2 - (X - ring.x) ** 2+0.001)
        Y = np.append(Y, (2 * ring.y - Y))
        X = np.append(X, X[::-1])
        ax.plot(
            X, Y, Z, color="black",  linewidth =1
        )


# Choosing set of parameters

#from Parameters import Rings

from Parameters_anisotropic import Rings

fig = plt.figure(figsize = (20, 20))
ax = fig.add_subplot(1, 1, 1, projection = '3d')

Orientation = "z"                                               # Orientation of ring that will be shown

for ring in Rings[:len(Rings)]:
    if ring.pos == Orientation and ring.z < 5:
        plot_circle(ring)

# Other parameters of plot
plt.title(f"Anisotropic structure", size = 36) #Size
ax.set_ylim(0, 9*12)
ax.set_xlabel("x", size = 36)
ax.set_ylabel("y", size = 36)
ax.set_zlabel("z", size = 36)

plt.savefig(f"Plots/Anisotropic2-rings")
plt.show()
