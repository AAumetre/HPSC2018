clear
mpicc -o Main explicit.c -lm -std=c99 -fopenmp
#mpiexec -np 4 ./Main param_tiny.dat
./Main param_tiny.dat
