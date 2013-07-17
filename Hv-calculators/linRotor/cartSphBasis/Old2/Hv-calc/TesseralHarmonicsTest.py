import SphFort
import numpy as np

for m in range(-5, 6):
    for l in range(abs(m), 6):
           print l, " ", m, " ", SphFort.ylm(l, m, 2.070796326794897, 2.070796326794897)
