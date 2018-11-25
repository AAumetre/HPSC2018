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
#include <malloc.h>
#include <stdbool.h>
#include <math.h>

#include "CSR_BSR.h"

/*=====================================================================================
* Provides examples for the use of the CSR/BSR library
=====================================================================================*/


int main(int argc, char **argv){

	/* Example codes for the CSR_BSR library */

	// Manually setting a BSR matrix
	/*
	1 2 3 3
	1 0 1 2 
	0 0 2 1
	0 0 0 6
	*/
	bsr_matrix mat; // Creates an empty BSR matrix
	bsr_init(&mat, 4,4,2,3); // Sets its size and block size
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

	// Prints out the BSR matrix
	printf("Prints out the BSR matrix, previously defined by hand\n");
	for (int j =0 ; j<4 ; ++j){
		for (int i=0 ; i<4 ; ++i){
			printf("%.0f ", bsr_get(&mat, j,i));
		}
		printf("\n");
	}
	bsr_free(&mat);

	// Let's define a matrix the usual way
	printf("Prints out the BSR matrix, previously defined from a natural matrix\n");
	double natural2[] = {1.0, 2.0, 3.0,   0.0, 0.0, 0.0,
						 1.0, 0.0, 1.0,   0.0, 0.0, 0.0,
						 1.0, 0.0, 2.0,   0.0, 0.0, 0.0,

						 0.0, 0.0, 0.0,   6.0, 1.0, 2.0,
						 0.0, 0.0, 0.0,   1.0, 1.0, 2.0,
					   	 0.0, 0.0, 0.0,   6.0, 1.0, 2.0};
	bsr_matrix mat3; 
	convert_natural_to_bsr(natural2, &mat3, sqrt(sizeof(natural2)/sizeof(natural2[0])), 3); // Convert it to BSR
	// And print it out to check
	for (int j =0 ; j<6 ; ++j){
		for (int i=0 ; i<6 ; ++i){
			printf("%.0f ", bsr_get(&mat3, j,i));
		}
		printf("\n");
	}
	printf("\n");


	// CSR vectors
	// Matrix vector product, vector norm
	csr_vector vec3;
	csr_vector vec4; // Used to store the result
	double vec_nat3[] =  {0,1,4,0,2,0};
	csr_vector_init(&vec3, vec_nat3, 6);

	bsr_matrix_vector(&mat3, &vec3, &vec4);
	printf("Does a BSR matrix, CSR vector product\n");
	for (int i = 0; i < 6; ++i){
		printf("%.0f\n", csr_vector_get(&vec4, i));
	}

	// Let's normalize that result
	double scale = 1/csr_vector_norm(&vec4);
	printf("\nScaling by a factor of %.4f (normalizing):\n", scale);
	csr_vector_scale(&vec4, scale);
	for (int i = 0; i < 6; ++i){
		printf("%.4f\n", csr_vector_get(&vec4, i));
	}
	printf("\n");
	bsr_free(&mat3);
	csr_vector_free(&vec3);
	csr_vector_free(&vec4);


	// Vector addition
	csr_vector vec5;
	csr_vector vec6;
	csr_vector vec7;
	double vec_nat5[] =  {0,1,4,3,1,0};
	double vec_nat6[] =  {0,1,4,0,0,0};
	csr_vector_init(&vec5, vec_nat5, 6);
	csr_vector_init(&vec6, vec_nat6, 6);

	printf("Adding two vectors");
	csr_vector_add(&vec5, &vec6, &vec7);
	printf("\n");
	for (int i = 0; i < 6; ++i){
		printf("%.0f\n", csr_vector_get(&vec7, i));
	}
	printf("\n");
	csr_vector_free(&vec5);
	csr_vector_free(&vec6);
	csr_vector_free(&vec7);


	// Setting values on a dynamically allocated vector
	printf("Dynamically allocates a CSR vector");
	csr_vector concentration;
	csr_vector_init_empty(&concentration, 500);
	// Setting a value
	for (int i = 0; i < 500; i+=1){ // Segmentation faults for different values of incrementation
		csr_vector_set(&concentration, i, i);
	}
	printf("\n");
	for (int j = 0; j < 20; ++j){
		for (int i = 0; i < 25; ++i){
			printf("%.0f ", csr_vector_get(&concentration, i+25*j));
		}
		printf("\n");
	}
	printf("\n");

	// findIndex function test
	printf("Testing the findIndex function:\n");
	int A[] = {0,1,2,3,7,8,9};
	for (int i = 0; i < 10; ++i){
		printf("They key %d can be found at index %d.\n", i, findIndex(A, i, 7));
	}
	printf("\n");

	// Test the function that gives each nodes its workload
	int *share = shareWorkload(41,16);
	printf("The best sharing is m=%d, n=%d\n", share[0], share[1]);

	// Key insertion algorithm test
	printf("\nKey insertion algorithm test\n");
	int rows[] = {0,1,2,5,8,9, 0};
	insertKey_int(rows, 7, 7);
	for (int i = 0; i < 7; ++i){
		printf("%d ", rows[i]);
	}

	// Resizing function test
	csr_vector test_vector;
	csr_vector_init(&test_vector, vec_nat5, 6);
	printf("\n\nSize before resizing: %d bytes\n", sizeof(test_vector.rows)); // doesn't work, return the size of the pointer
	//void csr_vector_resize(csr_vector *vector, int new_nnzb)


	return(0);
}
