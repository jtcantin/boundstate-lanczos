from sys import argv
import string
import SphFort
import numpy as np

l_max = string.atoi(argv[1])
theta = string.atof(argv[2])
phi = string.atof(argv[3])

for m in range(-l_max, l_max+1):
    for l in range(abs(m), l_max+1):
        print l, " ", m, " ", SphFort.ylm(l,m,theta,phi)
