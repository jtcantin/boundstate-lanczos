################################################################################
# Test program for the Module that contains functions for Lanczos Diagonalization
################################################################################
# Started on 18 Apr 2013 by Joshua Cantin
################################################################################

import LanczosDiag as LD

lm_n_table = LD.lm_in_n(3)

print lm_n_table
print len(lm_n_table)

print LD.n_to_lm(5)
print LD.n_to_lm(13)
print LD.n_to_lm(21)
print LD.n_to_lm(3)
