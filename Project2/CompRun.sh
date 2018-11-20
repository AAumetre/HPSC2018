clear
#mpicc -g -o Main explicit.c -lm -std=c99 -fopenmp
#mpiexec -np 12 time ./Main param_tiny.dat

gcc -o Main CSR_BSR_example.c -lm -std=c99
./Main
