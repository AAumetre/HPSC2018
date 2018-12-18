clear

#gcc CSR_BSR_example.c COO_CSR_BSR.c algorithms.c -o Main -lm -std=c99
#./Main

module load openmpi/1.8.4/gcc-4.9.2

mpicc -g -o Octopus octopus.c explicit.c implicit.c COO_CSR_BSR.c fileIO.c algorithms.c -lm -std=c99 -fopenmp
OMP_NUM_THREADS=1

#for ii in 1 2 3 4 5 6 7 8;
#do
#export OMP_NUM_THREADS=$ii
#t0=$(date +%s%3N)
#mpiexec -np 2 ./Octopus param_medium.dat 1
#t1=$(date +%s%3N)
  #echo "Elapsed wall clock time = " $(( $t1-$t0 )) " milliseconds."
#done


#mpicc -g -o Octopus octopus.c explicit.c implicit.c COO_CSR_BSR.c fileIO.c algorithms.c -lm -std=c99 -fopenmp
#mpiexec -np 2 time ./Octopus param_medium.dat 0

#mpicc -g -o Octopus octopus.c explicit.c implicit.c COO_CSR_BSR.c fileIO.c algorithms.c -lm -std=c99 -fopenmp
mpiexec -np 1 ./Octopus param_big.dat 1

#echo "Implicit";

#echo "Tiny";
#t0=$(date +%s%3N)
#mpiexec -np 4 ./Octopus param_tiny.dat 1;
#t1=$(date +%s%3N)
#echo "Elapsed wall clock time = " $(( $t1-$t0 )) " milliseconds."

#echo "Small";
#t0=$(date +%s%3N)
#mpiexec -np 4 ./Octopus param_small.dat 1;
#t1=$(date +%s%3N)
#echo "Elapsed wall clock time = " $(( $t1-$t0 )) " milliseconds."

#echo "Medium";
#t0=$(date +%s%3N)
#mpiexec -np 4 ./Octopus param_medium.dat 1;
#t1=$(date +%s%3N)
#echo "Elapsed wall clock time = " $(( $t1-$t0 )) " milliseconds."

#echo "Big";
#t0=$(date +%s%3N)
#mpiexec -np 4 ./Octopus param_big.dat 1;
#t1=$(date +%s%3N)
#echo "Elapsed wall clock time = " $(( $t1-$t0 )) " milliseconds."
