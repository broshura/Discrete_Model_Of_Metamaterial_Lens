# Visualization of each part of modeling

import matplotlib.pyplot as plt
from Parameters import Rings
from matplotlib import rcParams
from Simple import Impedance_real, Impedance_imag, Omega
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

# Making subplots and figure of Imaginary part of Impedance on responding ring

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

surf = ax.scatter(
[ring.x for ring in Rings if ring.pos==Orientation],            # X-data of centers
[ring.y for ring in Rings if ring.pos==Orientation],            # Y-data of centers
[ring.z for ring in Rings if ring.pos==Orientation],            # Z-data of centers
color = "b")                                                    # Coloring and so one

# Other parameters of plot
plt.title(f"Centers of all {Orientation}oriented rings", size = 36) #Size
ax.set_ylim(-18, 18)
ax.set_xlabel("x", size = 36)
ax.set_ylabel("y", size = 36)
ax.set_zlabel("z", size = 36)

#plt.show()

