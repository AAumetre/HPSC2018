/*=======================================================================================
*	This code was written by:                                                          *
*								Antonin Aumètre - antonin.aumetre@gmail.com            *
*								Céline Moureau -  cemoureau@gmail.com          		   *
*	For: High Performance Scientific course at ULiège, 2018-19                         *
*	Project 2                                                                          *
*                             														   *
*	Originally uploaded to: https://github.com/Cobalt1911                              *
*	Under GNU General Public License 11/2018                                           *
=======================================================================================*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>

#include "COO_CSR_BSR.h"
#include "algorithms.h"
#include "fileIO.h"

Param readDat(char *filename){
	FILE *file = fopen(filename, "r");
	Param parameters;

	if (file != NULL)
	{
		fscanf(file, "%lf %lf %lf %lf %lf %lf %lf %lf %u %lf", &parameters.h, &parameters.m, &parameters.L, &parameters.Tmax, &parameters.vx, &parameters.vy, &parameters.vz, &parameters.D, &parameters.S, &parameters.rthreshold);
		//printf("The parameters are : h %lf  m %lf  L %lf Tmax %lf %lf %lf %lf %lf %u %lf", parameters.h, parameters.m, parameters.L, parameters.Tmax, parameters.vx, parameters.vy, parameters.vz, parameters.D, parameters.S, parameters.rthreshold);

		fclose(file);
	}
	else
	{
		puts("File ERR0R !");
		exit(1);
	}

	return parameters;
}

void writeSequentialFile(double values[], int iteration, int N){
	char file_name[20];
	sprintf(file_name, "results/c_%d.dat",iteration);

	FILE *pNewFile = fopen(file_name, "wb");
	//Writing the header
	int header[] = {N};
	fwrite(header, sizeof(int), 1, pNewFile);

	fwrite(values, sizeof(double), N, pNewFile);
	fclose(pNewFile);
}