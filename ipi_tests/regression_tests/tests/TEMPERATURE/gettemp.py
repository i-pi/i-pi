import numpy as np
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser(description="Calculate the average temperature")
parser.add_argument('fin', help='output properties file')
args = parser.parse_args()

properties = ['Temperature', 'T (atom=O)', 'T (atom=H)', 'T (atom=0)', 'T (bead=1)', 'T (nm=0)', 'T (nm=4)', 'T (nm=7)']

data = np.genfromtxt(args.fin)
for i,pi in enumerate(properties):
    print(pi,'(K) : ', np.average(data[200:,i+3]))
    plt.plot(data[200:,1],data[200:,i+3],label=pi)
plt.show()
