clear
#mpicc -g -o Main explicitCSR.c -lm -std=c99 -fopenmp
#mpiexec -np 4 ./Main param_small.dat

#mpicc -g -o Main explicit.c -lm -std=c99 -fopenmp
#mpiexec -np 10 ./Main param_medium.dat

#mpicc -g -o Main explicitCSR.c -lm -std=c99 -fopenmp
#mpiexec -np 2 ./Main param_small_bigger.dat

gcc -o Main CSR_BSR_example.c -lm -std=c99
./Main
