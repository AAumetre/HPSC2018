#!/bin/bash
#
module load openmpi/1.8.4/gcc-4.9.2
#
#SBATCH --job-name=test_mpi
#SBATCH --mail-user=cmoureau@student.uliege.be
#SBATCH --mail-type=ALL
#SBATCH --output=scalability_strong_implicit.txt
#
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=16
#      ####################  maximum number of cores per node on this cluster
#                            to be sure to be alone on the node

#####  #SBATCH --exclusive
#####  #      ############# OR --exclusive to get a full node for the job

#SBATCH --time=30:00
#SBATCH --mem-per-cpu=3000


#mpicc -g -o Octopus octopus.c explicit.c implicit.c COO_CSR_BSR.c fileIO.c algorithms.c -O0 -lm -std=c99 -fopenmp

echo "Implicit";
for ii in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16;
echo -e "\n\n**" $ii;
OMP_NUM_THREADS=$ii
do
  for jj in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16;
  do
  echo -e "\n\n**" $jj;
  t0=$(date +%s%3N)
  mpiexec -np $jj ./Octopus param_big.dat 1;
  t1=$(date +%s%3N)
  echo "Elapsed wall clock time = " $(( $t1-$t0 )) " milliseconds."
  done
done
