# Discrete_Model_Of_Metamaterial_Lens
This repository contains program written on Python 3 to calculate the properties of finite metamaterial samples, which takes into account the discrete structure of metamaterials based on capacitively loaded ring resonators.

## Code

There are function and parameter files used in modeling.
Parameter files contain variable `name` that usually used in
the end of txt-data files.

### `Solving.py`
Main file that combine all part of modeling, use different 
method different sets of parameters for analyzing.

### `Impedance_matrix.py`
Calculating mutual part of impedance, that based only at geomerty
of the problem. Result saved in Data folder as `Data-name.txt` to avoid unnecessary repeating.

### `Comparing_matrixes.py`
To debug code results from integrating compared with verified matrix from matlab code, stored at `matrica.txt`

Modeling based on MRI-Lenz from "Exact modelling method for discrete finite
metamaterial lens
M. Lapine, L. Jelinek, R. Marque´s,
 M.J. Freire" experiment.
### `Simple_Method.py`
Solving matrix equation with currents, voltage and Impedance using Gauss method. 
The result saved in `Responding-name.txt` file as a table for real and imaginary part.

There necessary to use real parameters of each ring
and coils, that described in `Parameters-name.py` files.

### `Fast_Method.py`
Fast Method for solving matrix equation which uses symmetry of Impedance Matrix, so
in certain numeration it turns to 3D Toeplitz matrix and solved numerically on Krylov subspaces.

### `Plotting.py`
Compare results of modeling, using tracked data of plots from "Exact modelling method for discrete finite
metamaterial lens
M. Lapine, L. Jelinek, R. Marque´s,
 M.J. Freire" experiment and calculated frequency-response characteristics
from Data folder.

Tracked experiment data named `RespondingRingIm.txt` and `RespondingRingReal` for real and imaginary part,
modeling from matlab code based on similar integrate method have `model.txt` at the end.


### `Rings_visualize.py/.ipynb`
Visualizing of mutual location of rings in structure.
You can configure visibility of each layer and orientation using `Only_Border` and 
`Orientation`.

_To rotate structure straightaway use Jupiter (`.ipynb`)_.

### `Permeablity_anisotropic.py`
Calculate effective megnetic permeablity for structure anisotropic structure.

**(Does not veryfied)**

### `Ring_class.py`
Object class to simplify working with rings orientation and coordinates.

Each ring contain has to get its position, orientation and radius. To see each one you can use `print(ring)`.
### `Parameters_MRI.py`
Contains parameters of split-ring, responding coil and structure
for MRI-Lenz from "Exact modelling method for discrete finite
metamaterial lens
M. Lapine, L. Jelinek, R. Marque´s,
 M.J. Freire" experiment.

### `Parameters_anisotropic.py`

Parameter of researched anisotropic structure of non-capacitive rings,
where all rings has similar orientation and little distance
between layers.

### Data 
There is files with impedance matrixes(`Data-name.txt`) and their sums(`SumM-name.txt`) in format
for each set of parameters to avoid extra calcultions.

### Plots

There is plots of all studied dependencies and mutual location of rings in structure.

## Reports

This folder contains links for progress reports on overleaf, presentation and supplementary materials.

## Old_linza

There is matlab version of similar code for only MRI set of parameters to compare result with verified results.
