"""
Calculate the accumulated average of the exponential estimator <exp(-beta*spr)>
for the isotope fractionation.

accumulated.out:
first column : simulation time
second column : the instantanous exponential estimator
third column : the accumulated average of the exponential estimator
fourth column : 1000ln(<exp>)
5-7th column : for the other type of isotope substitution
"""

import numpy as np
from math import log

# the number of steps that allows the system to equilibrate
# for water liquid with pile_g thermostat, equilibration time ~ 0.1 ps
equilibriumstep = 4000

table = np.loadtxt("simulation.out")

(asize, bsize) = np.shape(table)
print(asize, bsize)

out = np.zeros((asize - equilibriumstep, 7))

count = 0

for a in range(asize):
    if a >= equilibriumstep:
        out[count][0] = table[a][1]
        out[count][1] = table[a][9]
        out[count][4] = table[a][12]
        if count == 0:
            out[count][2] = table[a][9]
            out[count][5] = table[a][12]
        else:
            out[count][2] = out[count - 1][2] * count / (count + 1) + table[a][9] / (
                count + 1
            )
            out[count][5] = out[count - 1][5] * count / (count + 1) + table[a][12] / (
                count + 1
            )
        out[count][3] = 1000.0 * log(out[count][2])
        out[count][6] = 1000.0 * log(out[count][5])
        count += 1

np.savetxt("accumulated.out", out)
