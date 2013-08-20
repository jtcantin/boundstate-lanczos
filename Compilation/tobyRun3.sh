
export OMP_NUM_THREADS=8
echo SPCElanczoslogBacicSmallCage0
time ./LanczosTOBYnoAcosIndex LanczosInputFileToby5.txt HvInputFileToby0.txt > /scratch/jtcantin/SPCElanczoslogBacicSmallCage0

export OMP_NUM_THREADS=8
echo SPCElanczoslogBacicSmallCage1
time ./LanczosTOBYnoAcosIndex LanczosInputFileToby5.txt HvInputFileToby1.txt > /scratch/jtcantin/SPCElanczoslogBacicSmallCage1

export OMP_NUM_THREADS=8
echo SPCElanczoslogBacicSmallCage2
time ./LanczosTOBYnoAcosIndex LanczosInputFileToby5.txt HvInputFileToby2.txt > /scratch/jtcantin/SPCElanczoslogBacicSmallCage2

export OMP_NUM_THREADS=8
echo SPCElanczoslogBacicSmallCage3
time ./LanczosTOBYnoAcosIndex LanczosInputFileToby5.txt HvInputFileToby3.txt > /scratch/jtcantin/SPCElanczoslogBacicSmallCage3

export OMP_NUM_THREADS=8
echo SPCElanczoslogBacicSmallCage4
time ./LanczosTOBYnoAcosIndex LanczosInputFileToby5.txt HvInputFileToby4.txt > /scratch/jtcantin/SPCElanczoslogBacicSmallCage4

export OMP_NUM_THREADS=8
echo SPCElanczoslogBacicSmallCage5
time ./LanczosTOBYnoAcosIndex LanczosInputFileToby5.txt HvInputFileToby5.txt > /scratch/jtcantin/SPCElanczoslogBacicSmallCage5

export OMP_NUM_THREADS=8
echo SPCElanczoslogBacicSmallCage6
time ./LanczosTOBYnoAcosIndex LanczosInputFileToby5.txt HvInputFileToby6.txt > /scratch/jtcantin/SPCElanczoslogBacicSmallCage6

export OMP_NUM_THREADS=8
echo SPCElanczoslogBacicSmallCage7
time ./LanczosTOBYnoAcosIndex LanczosInputFileToby5.txt HvInputFileToby7.txt > /scratch/jtcantin/SPCElanczoslogBacicSmallCage7
