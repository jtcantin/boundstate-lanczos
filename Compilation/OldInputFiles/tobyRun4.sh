
export OMP_NUM_THREADS=8
echo SPCElanczoslogBacicSmallCage12
time ./LanczosTOBYgrid LanczosInputFileToby5.txt HvInputFileToby12.txt > /scratch/jtcantin/SPCElanczoslogBacicSmallCage12

export OMP_NUM_THREADS=8
echo SPCElanczoslogBacicSmallCage13
time ./LanczosTOBYgrid LanczosInputFileToby5.txt HvInputFileToby13.txt > /scratch/jtcantin/SPCElanczoslogBacicSmallCage13

export OMP_NUM_THREADS=8
echo SPCElanczoslogBacicSmallCage14
time ./LanczosTOBYgrid LanczosInputFileToby5.txt HvInputFileToby14.txt > /scratch/jtcantin/SPCElanczoslogBacicSmallCage14

export OMP_NUM_THREADS=8
echo SPCElanczoslogBacicSmallCage15
time ./LanczosTOBYgrid LanczosInputFileToby5.txt HvInputFileToby15.txt > /scratch/jtcantin/SPCElanczoslogBacicSmallCage15
