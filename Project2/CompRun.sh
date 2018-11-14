clear
mpicc -o Main explicit.c -lm -std=c99 -fopenmp
mpiexec -np 1 ./Main param_tiny.dat
