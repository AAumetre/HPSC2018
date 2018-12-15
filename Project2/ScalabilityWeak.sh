#!/bin/bash
#
module load openmpi/1.8.4/gcc-4.9.2
#
#SBATCH --job-name=test_mpi
#SBATCH --mail-user=cmoureau@student.uliege.be
#SBATCH --mail-type=ALL
#SBATCH --output=scalability_weak.txt
#
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=4
#      ####################  maximum number of cores per node on this cluster
#                            to be sure to be alone on the node

#####  #SBATCH --exclusive
#####  #      ############# OR --exclusive to get a full node for the job

#SBATCH --time=5:00
#SBATCH --mem-per-cpu=3000


#mpicc -g -o Octopus octopus.c explicit.c implicit.c COO_CSR_BSR.c fileIO.c algorithms.c -O0 -lm -std=c99 -fopenmp

echo "Explicit";

echo "Tiny";
t0=$(date +%s%3N)
mpiexec ./Octopus param_tiny.dat 0;
t1=$(date +%s%3N)
echo "Elapsed wall clock time = " $(( $t1-$t0 )) " milliseconds."

echo "Small";
t0=$(date +%s%3N)
mpiexec ./Octopus param_small.dat 0;
t1=$(date +%s%3N)
echo "Elapsed wall clock time = " $(( $t1-$t0 )) " milliseconds."

echo "Medium";
t0=$(date +%s%3N)
mpiexec ./Octopus param_medium.dat 0;
t1=$(date +%s%3N)
echo "Elapsed wall clock time = " $(( $t1-$t0 )) " milliseconds."

echo "Big";
t0=$(date +%s%3N)
mpiexec ./Octopus param_big.dat 0;
t1=$(date +%s%3N)
echo "Elapsed wall clock time = " $(( $t1-$t0 )) " milliseconds."
