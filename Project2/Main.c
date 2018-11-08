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

#include "CSR_BSR.h"

#define DEBUG 0

/* STRUCTURES AND FUNCTIONS TO BE IMPLEMENTED
* Vector sum
*/



	int main(int argc, char **argv){
		/*
		bsr_matrix mat; // Creates an empty BSR matrix
		/*
		1 2 3 3
		1 0 1 2 
		0 0 2 1
		0 0 0 6
		*/
		/*
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
  				printf("%.0f ", bsr_get(&mat, j,i));
  			}
  			printf("\n");
  		}
		bsr_free(&mat);

  		double natural2[] = {1.0, 2.0, 3.0,   0.0, 0.0, 0.0,
						     1.0, 0.0, 1.0,   0.0, 0.0, 0.0,
						     1.0, 0.0, 2.0,   0.0, 0.0, 0.0,

						     0.0, 0.0, 0.0,   6.0, 1.0, 2.0,
						     0.0, 0.0, 0.0,   1.0, 1.0, 2.0,
						     0.0, 0.0, 0.0,   6.0, 1.0, 2.0};
		bsr_matrix mat3; 
		natural_to_bsr(natural2, &mat3, sqrt(sizeof(natural2)/sizeof(natural2[0])), 3);
		for (int j =0 ; j<6 ; ++j){
  			for (int i=0 ; i<6 ; ++i){
  				printf("%.0f ", bsr_get(&mat3, j,i));
  			}
  			printf("\n");
  		}
  		bsr_free(&mat3);
  		*/

  		// Vectors
  		/*csr_vector vec;
  		csr_vector vec2;
  		double vec_nat[] =  {1,0,2,0,0,1};
  		double vec_nat2[] = {1,1,1,0,3,0};
  		csr_vector_init(&vec, vec_nat, 6);
  		csr_vector_init(&vec2, vec_nat2, 6);

  		csr_vector_free(&vec);
  		csr_vector_free(&vec2);*/

  		// Matrix vector product
  		/*csr_vector vec3;
  		csr_vector vec4;
  		double vec_nat3[] =  {0,1,4,0,0,0};
  		csr_vector_init(&vec3, vec_nat3, 6);

  		double natural3[] = {1.0, 2.0, 3.0,   0.0, 0.0, 0.0,
						     1.0, 0.0, 1.0,   0.0, 0.0, 0.0,
						     1.0, 0.0, 2.0,   0.0, 0.0, 0.0,

						     0.0, 0.0, 0.0,   6.0, 1.0, 2.0,
						     0.0, 0.0, 0.0,   1.0, 1.0, 2.0,
						     0.0, 0.0, 0.0,   6.0, 1.0, 2.0};
		bsr_matrix mat4; 
		natural_to_bsr(natural3, &mat4, sqrt(sizeof(natural3)/sizeof(natural3[0])), 3);

		bsr_matrix_vector(&mat4, &vec3, &vec4);
  		printf("Result:\n");
		for (int i = 0; i < 6; ++i){
  			printf("%.0f\n", csr_vector_get(&vec4, i));
  		}
  		double scale = 0.5;
  		printf("Scaling by a factor of %.2f:\n", scale);
  		csr_vector_scale(&vec4, scale);
		for (int i = 0; i < 6; ++i){
  			printf("%.0f\n", csr_vector_get(&vec4, i));
  		}

		bsr_free(&mat4);
  		csr_vector_free(&vec3);
  		csr_vector_free(&vec4);
  		*/

		

	//======================= PRE-PROCESSING ============================//

	//======================= ALGORITHM =================================//

	//======================= POST-PROCESSING ===========================//
	
	//======================= END OF PROGRAM ============================//

	return(0);
}
