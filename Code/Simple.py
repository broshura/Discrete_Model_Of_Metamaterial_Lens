from Parameters import Z_0, V
import numpy as np
from scipy import linalg
with open("DATA/Data.txt", "r") as res:
    RES = res.read()
    Z = [[k for k in x.split(" ")] for x in RES.split("\n")]
Z = [[complex(Z[i][k]) for k in range(len(Z[i]))] for i in range(len(Z)-1)]
Z = np.array(Z)
Z = Z + np.eye(len(Z))* Z_0
I = linalg.solve(Z, V)
