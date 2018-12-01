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
#include <math.h>
#include <time.h>

#include "COO_CSR_BSR.h"
#include "algorithms.h"
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

	bsr_matrix_csr_vector(&mat3, &vec3, &vec4);
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

	// CSR equals function test
	csr_vector vec10;
	csr_vector_init(&vec10, vec_nat5, 6);
	csr_vector vec11;
	csr_vector_init_empty(&vec11, 6);
	printf("\n\nCSR equals function test\nFirst vector:      ");
	for (int i = 0; i < 6; ++i){
		printf("%.1f ", csr_vector_get(&vec10, i));
	}
	csr_vector_equals(&vec11, &vec10);
	printf("\nAfter affectation: ");
	for (int i = 0; i < 6; ++i){
		printf("%.1f ", csr_vector_get(&vec11, i));
	}
	printf("\n");

	/* =============================== */
	// Natural vectors
	nat_vector nat1;
	nat_vector_init(&nat1, 6);
	for (int i=0; i < 6; ++i){
		nat1.values[i] = 10*i;
	}
	printf("\nNatural vector: \n");
	for (int i = 0; i < 6; ++i){
		printf("%.1f ", nat1.values[i]);
	}
	printf("\nEuclidian norm: %.2f\n", nat_vector_norm(&nat1));
	nat_vector_scale(&nat1, 1/nat_vector_norm(&nat1));
	printf("Normalized natural vector: \n");
	for (int i = 0; i < 6; ++i){
		printf("%.1f ", nat1.values[i]);
	}
	printf("\n\n");

	// COO matrix
	printf("COO matrices test\n");
	coo_matrix coo1;
	coo_matrix_init_empty(&coo1, 6, 6);
	convert_natural_to_coo(natural2, &coo1, 36);
	printf("Rows: ");
	for (int i = 0; i < coo1.nnzb; ++i){
		printf("%d ", coo1.rows[i]);
	}
	printf("\nColumns: ");
	for (int i = 0; i < coo1.nnzb; ++i){
		printf("%d ", coo1.columns[i]);
	}
	printf("\nValues: ");
	for (int i = 0; i < coo1.nnzb; ++i){
		printf("%.0f ", coo1.values[i]);
	}
	printf("\n");
	for (int i = 0; i < 6; ++i){
		for (int j = 0; j < 6; ++j){
			printf("%.0f ", coo_matrix_get(&coo1, i, j));
		}
		printf("\n");
	}

	printf("\nCOO matrix and natural vector product\n");
	for (int i = 0; i < 6; ++i){
		nat1.values[i] = 1;
	}
	nat_vector nat2;
	nat_vector_init(&nat2, 6);
	coo_matrix_nat_vector(&coo1, &nat1, &nat2);
	for (int i = 0; i < 6; ++i){
		printf("%.1f ", nat2.values[i]);
	}
	printf("\n");

	//nat_vector_free(&nat1);
	nat_vector_free(&nat2);
	coo_matrix_free(&coo1);

	//Performance test for the COO/natural product.
	printf("\nPerformance test for the COO/natural product.\n");
	nat_vector nat3;
	nat_vector_init(&nat3, 12);
	nat_vector res3;
	nat_vector_init(&res3, 12);

	for (int i=0; i < 12; ++i){
		nat3.values[i] = i;
	}
	coo_matrix coo_mat;
	coo_matrix_init_empty(&coo_mat, 12, 12);
	double natural8[] = {1.0, 2.0, 3.0,   0.0, 0.0, 0.0,   1.0, 2.0, 3.0,   1.0, 0.0, 0.0,

						 1.0, 0.0, 1.0,   0.0, 3.0, 0.0,   1.0, 2.0, 3.0,   0.0, 6.0, 0.0,

						 1.0, 3.0, 2.0,   0.0, 0.0, 0.0,   6.0, 7.0, 3.0,   0.0, 0.0, 1.0,

						 0.0, 0.0, 0.0,   6.0, 1.0, 4.0,   1.0, 2.0, 3.0,   0.0, 0.5, 1.0,

						 0.0, 1.0, 0.0,   1.0, 1.0, 2.0,   0.0, 2.0, 3.0,   0.0, 1.0, 2.0,

					   	 0.0, 0.0, 0.0,   6.0, 1.0, 2.0,   1.0, 5.0, 3.0,   2.0, 0.0, 0.0,

					   	 1.0, 2.0, 3.0,   0.0, 0.0, 0.0,   2.0, 2.0, 3.0,   0.0, 0.8, 1.4,

						 1.0, 0.0, 1.0,   0.0, 1.0, 5.0,   1.0, 2.0, 3.0,   8.0, 4.0, 0.0,

						 1.0, 0.0, 2.0,   0.0, 0.0, 0.0,   1.0, 4.0, 3.0,   0.0, 0.0, 1.0,

						 0.0, 1.0, 0.0,   6.0, 0.0, 2.0,   1.0, 7.0, 3.0,   1.0, 0.1, 0.4,

						 0.0, 0.0, 0.0,   1.0, 1.0, 2.0,   1.0, 2.0, 3.0,   0.0, 4.0, 0.0,

					   	 0.0, 0.0, 0.0,   6.0, 1.0, 2.0,   1.0, 2.0, 3.0,   0.0, 0.0, 0.0};
	convert_natural_to_coo(natural8, &coo_mat, 12*12);
	coo_matrix_scale(&coo_mat, 5e-4);

	for (int i = 0; i < 12; ++i){
		for (int j = 0; j < 12; ++j){
			printf("%.5f ", coo_matrix_get(&coo_mat, i, j));
		}
		printf("\n");
	}

	clock_t begin = clock();
	clock_t end = clock();
	double time_spent;
	// Apply the matrix n times to the vector
	for (int i = 0; i < 1000; ++i){
		coo_matrix_nat_vector(&coo_mat, &nat3, &res3);
		nat_vector_equals(&nat3, &res3);
	}
	for (int i=0; i < 12; ++i){
		printf("%f\n", res3.values[i]);
	}


	end = clock();
	time_spent = (double)(end - begin) / (CLOCKS_PER_SEC);
	printf("\nJob done in %2.4lf s.\n", time_spent);
	
	return(0);
}
