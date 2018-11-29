clear
#mpicc -g -o Main explicitCSR.c -lm -std=c99 -fopenmp
#mpiexec -np 4 ./Main param_small.dat

#mpicc -g -o Main explicit.c CSR_BSR.c algorithms.c  -lm -std=c99 -fopenmp
#mpiexec -np 2 ./Main param_small.dat 1

#mpicc -g -o Main explicitCSR.c -lm -std=c99 -fopenmp
#mpiexec -np 2 ./Main param_small.dat

gcc CSR_BSR_example.c CSR_BSR.c algorithms.c -o Main -lm -std=c99
./Main
