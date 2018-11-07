/*=======================================================================================
*	This code was written by: 
*								Antonin Aumètre - antonin.aumetre@gmail.com
*								Céline Moureau -  cemoureau@gmail.com
*	For: High Performance Scientific course at ULiège, 2018-19
*	Project 2
*
*	Under GNU General Public License 11/2018
=======================================================================================*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <malloc.h>
#include <stdbool.h>

#include "matrices.h"

#define DEBUG 0

/* STRUCTURES AND FUNCTIONS TO BE IMPLEMENTED
* Sparse matrices structure cf. MPI class
*
* Conversion function, from natural matrix to sparse ?
* Scalar product
* Vector sparse matrix product
* Scalar vector product
* Vector sum
* Euclidian norm
*/



	int main(int argc, char **argv){

		bsr_matrix mat; // Creates an empty BSR matrix
		/*
		1 2 3 3
		1 0 1 2 
		0 0 2 1
		0 0 0 6
		*/
		bsr_init(&mat, 4,4,2,3); // Sets its size and block size
		// Do stuff with it
		mat.block_row_offsets[0]=0;
		mat.block_row_offsets[1]=2;
		mat.block_row_offsets[2]=3;
  		mat.block_columns[0] = 0;
  		mat.block_columns[1] = 1;
  		mat.block_columns[2] = 1;

  		mat.values[0] = 1;
  		mat.values[1] = 2;
  		mat.values[2] = 1;
  		mat.values[3] = 0;

  		mat.values[4] = 3;
  		mat.values[5] = 3;
  		mat.values[6] = 1;
  		mat.values[7] = 2;

  		mat.values[8] = 2;
  		mat.values[9] = 1;
  		mat.values[10] = 0;
  		mat.values[11] = 6;

  		for (int j =0 ; j<4 ; ++j){
  			for (int i=0 ; i<4 ; ++i){
  				printf("%f ", bsr_get(&mat, j,i));
  			}
  			printf("\n");
  			printf("\n");
  			printf("\n");
  		}

  		printf("Value: %f\n", bsr_get(&mat, 0,3));

		bsr_free(&mat);

		double natural[] = {1.0, 2.0, 3.0, 3.0,
						    1.0, 0.0, 1.0, 2.0,
						    0.0, 0.0, 2.0, 1.0,
						    0.0, 0.0, 0.0, 6.0};
		bsr_matrix mat2; 
		natural_to_bsr(natural, &mat2, 4, 2);

	//======================= PRE-PROCESSING ============================//

	//======================= ALGORITHM =================================//

	//======================= POST-PROCESSING ===========================//
	
	//======================= END OF PROGRAM ============================//

	return(0);
}
