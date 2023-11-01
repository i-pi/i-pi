import numpy as np

table1 = np.loadtxt("accumulated.out.w")
table2 = np.loadtxt("accumulated.out.v")


(asize1, bsize) = np.shape(table1)
(asize2, bsize) = np.shape(table2)

if asize1 > asize2:
    asize = asize2
else:
    asize = asize1

out = np.zeros((asize, 3))

count = 0

for a in range(asize):
    out[a][0] = table1[a][0]
    out[a][1] = table1[a][3] - table2[a][3]
    out[a][2] = table1[a][6] - table2[a][6]


np.savetxt("1000loga.out", out)
