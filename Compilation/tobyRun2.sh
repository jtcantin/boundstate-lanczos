
export OMP_NUM_THREADS=8
echo Alavi_SPCE_smallCage11_ortho_x04_dx002_l5_tp70_cp18
time ./LanczosTOBY LanczosInputFileToby.txt HvInputFileToby12.txt > /scratch/jtcantin/Alavi_SPCE_smallCage11_ortho_x04_dx002_l5_tp70_cp18

export OMP_NUM_THREADS=8
echo SmallBasisTest
time ./LanczosTOBY LanczosInputFileTobySmallBasisTest2.txt HvInputFileTobySmallBasisTest2.txt > /scratch/jtcantin/SmallBasisTest2

