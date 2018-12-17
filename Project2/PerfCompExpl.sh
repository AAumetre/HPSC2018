#!/bin/bash
#
module load openmpi/1.8.4/gcc-4.9.2
#
#SBATCH --job-name=test_mpi
#SBATCH --mail-user=cmoureau@student.uliege.be
#SBATCH --mail-type=ALL
#SBATCH --output=perfComp_expl.txt
#
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=4
#      ####################  maximum number of cores per node on this cluster
#                            to be sure to be alone on the node

#####  #SBATCH --exclusive
#####  #      ############# OR --exclusive to get a full node for the job

#SBATCH --time=20:00
#SBATCH --mem-per-cpu=3000


#mpicc -g -o Octopus octopus.c explicit.c implicit.c COO_CSR_BSR.c fileIO.c algorithms.c -O0 -lm -std=c99 -fopenmp

echo "Explicit";
OMP_NUM_THREADS=4

echo "Big";
t0=$(date +%s%3N)
mpiexec --bind-to none -np 4 ./Octopus param_big.dat 0;
t1=$(date +%s%3N)
echo "Elapsed wall clock time = " $(( $t1-$t0 )) " milliseconds."

echo "D greater";
t0=$(date +%s%3N)
mpiexec --bind-to none -np 4 ./Octopus param_big_Dgreater.dat 0;
t1=$(date +%s%3N)
echo "Elapsed wall clock time = " $(( $t1-$t0 )) " milliseconds."

echo "D smaller";
t0=$(date +%s%3N)
mpiexec --bind-to none -np 4 ./Octopus param_big_Dsmaller.dat 0;
t1=$(date +%s%3N)
echo "Elapsed wall clock time = " $(( $t1-$t0 )) " milliseconds."

echo "V greater";
t0=$(date +%s%3N)
mpiexec --bind-to none -np 4 ./Octopus param_big_Vgreater.dat 0;
t1=$(date +%s%3N)
echo "Elapsed wall clock time = " $(( $t1-$t0 )) " milliseconds."

echo "V really greater";
t0=$(date +%s%3N)
mpiexec --bind-to none -np 4 ./Octopus param_big_VreallyGreater.dat 0;
t1=$(date +%s%3N)
echo "Elapsed wall clock time = " $(( $t1-$t0 )) " milliseconds."
