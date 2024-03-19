# 3D visualization of structure

import matplotlib.pyplot as plt

from matplotlib import rcParams
import numpy as np
from numpy import sqrt, pi

import plotly.graph_objects as go

# Function for drawing circles in 3D

def plot_circle(ring):
    n = 250
    if ring.r == 2.54/1.5*3:
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
            X, Y, Z, color = "blue", linewidth = 2
        )
    elif ring.pos == "y":
        Y = np.array([ring.y for i in range(n)])
        Z = np.linspace(ring.z - ring.r, ring.z + ring.r, n//2)
        X = ring.x + sqrt(ring.r ** 2 - (Z-ring.z) ** 2+0.001)
        X = np.append(X, 2*ring.x-X)
        Z = np.append(Z, Z[::-1])
        ax.plot(
            X, Y, Z, color = "red",  linewidth = 2
        )
    elif ring.pos == "z":
        Z = np.array([ring.z for i in range(n)])
        X = np.linspace(ring.x - ring.r, ring.x + ring.r, n//2)
        Y = ring.y + sqrt(ring.r ** 2 - (X - ring.x) ** 2+0.001)
        Y = np.append(Y, (2 * ring.y - Y))
        X = np.append(X, X[::-1])
        ax.plot(
            X, Y, Z, color="black",  linewidth = 2
        )


# Choosing set of parameters

#from Parameters import Rings, name

from Parameters_anisotropic import Rings, delta_x, delta_z, name

fig = plt.figure(figsize = (20, 20))
ax = fig.add_subplot(1, 1, 1, projection = '3d')


# Draw only one orientation (any for all)
# and show only bordered rings
# to decrease number of rings

Orientation = "any"
OnlyBorder = False


Border_z = max([ring.z for ring in Rings])
for ring in Rings[:len(Rings)]:
    id_d = ring.pos == Orientation or Orientation == "any"

    Layer_z = [RING for RING in Rings if RING.z == ring.z]
    Border_x = max([ring.x for ring in Layer_z])
    Border_y = min([ring.y for ring in Layer_z])


    id_x = ring.x == Border_x or not OnlyBorder
    id_y = ring.y == Border_y or not OnlyBorder
    id_z = ring.z == Border_z or not OnlyBorder

    if id_d and (id_x or id_y or id_z):
        plot_circle(ring)

# Other parameters of plot
plt.title(f"{name} structure", size = 36) #Size
#plt.ylim(-18, 18)
ax.set_xlabel("x", size = 36)
ax.set_ylabel("y", size = 36)
ax.set_zlabel("z", size = 36)
#ax.view_init(elev=2, azim=90)
plt.savefig(f"Plots/3D-structure/{name}-rings")
plt.show()
