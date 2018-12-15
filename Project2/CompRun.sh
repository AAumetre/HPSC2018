clear

#gcc CSR_BSR_example.c COO_CSR_BSR.c algorithms.c -o Main -lm -std=c99
#./Main

#module load openmpi/1.8.4/gcc-4.9.2
OMP_NUM_THREADS=2

t0=$(date +%s%3N)
mpicc -g -o Octopus octopus.c explicit.c implicit.c COO_CSR_BSR.c fileIO.c algorithms.c -lm -std=c99 -fopenmp
mpiexec -np 2 ./Octopus param_medium.dat 1
t1=$(date +%s%3N)
  echo "Elapsed wall clock time = " $(( $t1-$t0 )) " milliseconds."

#mpicc -g -o Octopus octopus.c explicit.c implicit.c COO_CSR_BSR.c fileIO.c algorithms.c -lm -std=c99 -fopenmp
#mpiexec -np 2 time ./Octopus param_medium.dat 0

#mpicc -g -o Octopus octopus.c explicit.c implicit.c COO_CSR_BSR.c fileIO.c algorithms.c -lm -std=c99 -fopenmp
#mpiexec -np 4 ./Octopus param_medium.dat 1
