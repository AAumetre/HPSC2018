clear
#mpicc -g -o Main explicit.c COO_CSR_BSR.c algorithms.c  -lm -std=c99 -fopenmp
#mpiexec -np 4 ./Main param_medium.dat 1

#mpicc -g -o Main explicitCSR.c COO_CSR_BSR.c algorithms.c -lm -std=c99 -fopenmp
#mpiexec -np 2 ./Main param_small.dat

#mpicc -g -o Main implicit.c COO_CSR_BSR.c algorithms.c -lm -std=c99 -fopenmp
#mpiexec -np 2 ./Main param_small.dat 1

#gcc CSR_BSR_example.c COO_CSR_BSR.c algorithms.c -o Main -lm -std=c99
#./Main

mpicc -g -o Octopus octopus.c explicit.c implicit.c COO_CSR_BSR.c fileIO.c algorithms.c -lm -std=c99 -fopenmp
mpiexec -np 2 ./Octopus param_medium.dat 0

mpicc -g -o Octopus octopus.c explicit.c implicit.c COO_CSR_BSR.c fileIO.c algorithms.c -lm -std=c99 -fopenmp
mpiexec -np 4 ./Octopus param_medium.dat 0
