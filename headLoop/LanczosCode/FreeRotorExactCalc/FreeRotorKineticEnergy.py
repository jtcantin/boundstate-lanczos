import numpy as np
import sys

numPoints = int(sys.argv[1])

mass_gmol = 2.015650064 #g/mol
x_max = 1e-9 #m
l_max = 0
RotConst = 0.709652722117881 #kJ/mol

hBar = 1.054572E-34 #kJ/mol
Nav = 6.02214E+23 #/mol
J_kJ = 1000. #J/kJ
g_kg = 1000. #g/kg

kJmol_J = 1. / J_kJ * Nav

print kJmol_J

mass_kg = mass_gmol / g_kg / Nav

dx = x_max / numPoints

maxRotKinEng = RotConst * l_max * (l_max + 1)

print "Mass: ", mass_kg, "kg"
print "x_max: ", x_max, "m"
print "dx: ", dx, "m"
print "l_max: ", l_max
print "Maximum Rotational Kinetic Energy: ", maxRotKinEng, "kJ/mol"

kinEngPreFactor = hBar**2 / 2. / mass_kg / dx**2

kinEngPreFactor_kJmol = kinEngPreFactor * kJmol_J
print "Kinetic Energy Pre-Factor: ", kinEngPreFactor_kJmol, "kJ/mol"

T = np.random.rand(numPoints,numPoints)
for i in range(0,numPoints):
    for ip in range(0,numPoints):
        if i==ip:
            T[i,ip] = kinEngPreFactor_kJmol * (-1.)**(i-ip) * np.pi * np.pi / 3
        elif i!=ip:
            T[i,ip] = kinEngPreFactor_kJmol * (-1.)**(i-ip) * 2 / (i-ip)**2

eigVal, eigVec = np.linalg.eig(T)

#T_X2D = np.kron(T,np.eye(numPoints))
#T_X3D = np.kron(T_X2D,np.eye(numPoints))

#T_Y2D = np.kron(np.eye(numPoints),T)
#T_Y3D = np.kron(T_Y2D,np.eye(numPoints))

#T_Z2D = np.kron(np.eye(numPoints),np.eye(numPoints))
#T_Z3D = np.kron(T_Z2D,T)

#T_3D = T_X3D + T_Y3D + T_Z3D

#eigVal3D, eigVec3D = np.linalg.eig(T_3D)

#print "Kinetic Energy Eigenvalues (1-D): ",  sorted(eigVal)[0:5]
#print "Distinct Kinetic Energy Eigenvalues (3-D): ", sorted(3*eigVal)[0:10]
#print "Kinetic Energy Eigenvalues (3-D): ", sorted(eigVal3D)[0:30]
#print eigVec

print "Kinetic Energy Eigenvalues (1-D): "
for i in range(0,numPoints):
    print i, " ", sorted(eigVal)[i]

print " "
print " "

#print "Kinetic Energy Eigenvalues (3-D): "
#for i in range(0,numPoints**3):
    #    print i, " ", sorted(eigVal3D)[i]
            
