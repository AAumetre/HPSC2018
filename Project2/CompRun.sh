clear
#mpicc -g -o Main explicitCSR.c -lm -std=c99 -fopenmp
#mpiexec -np 4 time ./Main param_small.dat

mpicc -g -o Main explicit.c -lm -std=c99 -fopenmp
mpiexec -np 14 time ./Main param_medium.dat

#mpicc -g -o Main explicitCSR.c -lm -std=c99 -fopenmp
#mpiexec -np 1 time ./Main param_medium.dat

#gcc -o Main CSR_BSR_example.c -lm -std=c99
#./Main
