clear
#mpicc -g -o Main explicitCSR.c -lm -std=c99 -fopenmp
#mpiexec -np 10 time ./Main param_big.dat

mpicc -g -o Main explicit.c -lm -std=c99 -fopenmp
mpiexec -np 10 time ./Main param_medium.dat

#gcc -o Main CSR_BSR_example.c -lm -std=c99
#./Main
