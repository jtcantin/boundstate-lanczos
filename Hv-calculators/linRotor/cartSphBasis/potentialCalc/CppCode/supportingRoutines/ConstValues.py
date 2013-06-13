from MMTK import *
import numpy as np

print np.pi
print Units.eps0

epsNew = 8.854187817e-12 #F/m, NIST
Na = 6.02214e23 #/mol, NIST
eC = 1.602177e-19 #C, NIST
epsMMTK = epsNew * (1e12)**2 /((eC)**2) /Na *(1./1000.) * (1e-9) **3

print epsNew
print epsMMTK

hbar_MMTK = Units.hbar
h_SI = 6.62606957e-34 #J.s, NIST
hbar_SI = h_SI /(2.0*np.pi)
hbar_kJmolps = hbar_SI / 1000 * Na * 1e12 # hbar_SI / 1000[J/kJ] * Na[/mol] * 1e12[ps/s]

print hbar_MMTK
print hbar_SI
print hbar_kJmolps
