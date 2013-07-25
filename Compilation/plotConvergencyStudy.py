import numpy as np
import sys
import string
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata
from mpl_toolkits.mplot3d import Axes3D
import scipy as sp
import scipy.ndimage

def gradient2D(x_data, y_data, z_data):

    zMat = [[ 0.0 for i in range(0, max(x_data))] for j in range(0, max(y_data))]

    #print zMat

    #Make a z 2D array
    l=0
    for i in range(0, max(x_data)):
        for j in range(0, max(y_data)):
            zMat[i][j] = z_data[l]
            l += 1
            
    zArray = np.array(zMat, dtype=np.float)
    
    zLaplacian = np.zeros(shape=(max(x_data),max(y_data)))

    sp.ndimage.filters.laplace(zArray, output=zLaplacian, mode='constant', cval=0.0)

    #zGradients = np.gradient(zArray)

    #zGradientxMag2 = np.dot(zGradients[0], zGradients[0])
    # zGradientyMag2 = np.dot(zGradients[1], zGradients[1])

    #zGradientMag2 = zGradientxMag2 + zGradientyMag2
    #zGradientMag2 = 

    #print zGradients
    print np.unravel_index(np.abs(zLaplacian).argmin(), zLaplacian.shape)
    print zLaplacian.argmin()


################
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

gradient2D(na, nb, energy)




#Plot the values
#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')

#naMesh, nbMesh = np.meshgrid(na, nb)
#energyData = griddata(na, nb, energy, naMesh, nbMesh)

#ax.plot_wireframe(naMesh, nbMesh, energyData)
#ax.set_xlabel("Theta Points")
#ax.set_ylabel("Phi Points")
#ax.set_zlabel("Expectation Value (kJ/mol)")

#plt.show()

