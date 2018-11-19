clear
mpicc -g -o Main explicit.c -lm -std=c99 -fopenmp
mpiexec -np 3 ./Main param_tiny.dat
