clear
mpicc -g -o Main explicit.c COO_CSR_BSR.c algorithms.c  -lm -std=c99 -fopenmp
mpiexec -np 10 ./Main param_big.dat 1

#mpicc -g -o Main explicitCSR.c COO_CSR_BSR.c algorithms.c -lm -std=c99 -fopenmp
#mpiexec -np 2 ./Main param_small.dat

#mpicc -g -o Main implicit.c COO_CSR_BSR.c algorithms.c -lm -std=c99 -fopenmp
#mpiexec -np 2 ./Main param_small.dat 1

#gcc CSR_BSR_example.c COO_CSR_BSR.c algorithms.c -o Main -lm -std=c99
#./Main
