
export OMP_NUM_THREADS=8
echo SmallBasisTest
time ./LanczosTOBY LanczosInputFileTobySmallBasisTest2.txt HvInputFileTobySmallBasisTest2.txt > tests

