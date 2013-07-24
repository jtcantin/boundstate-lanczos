import numpy as np
import sys
import string
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata
from mpl_toolkits.mplot3d import Axes3D

dataFilename = sys.argv[1]

dataFile = open(dataFilename, 'r')
junk = dataFile.readline()

na = []
nb = []
energy = []

for line in dataFile:
    pieces = string.split(line);
    na.append(int(pieces[0]))
    nb.append(int(pieces[1]))
    energy.append(float(pieces[2]))


    #for i in range(0,len(na)):
    #print na[i], " ", nb[i], " ", energy[i]

#Plot the values
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

naMesh, nbMesh = np.meshgrid(na, nb)
energyData = griddata(na, nb, energy, naMesh, nbMesh)

ax.plot_wireframe(naMesh, nbMesh, energyData)
ax.set_xlabel("Theta Points")
ax.set_ylabel("Phi Points")
ax.set_zlabel("Expectation Value (kJ/mol)")

plt.show()
