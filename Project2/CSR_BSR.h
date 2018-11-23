/*=======================================================================================
*	This code was written by:                                                          *
*								Antonin Aumètre - antonin.aumetre@gmail.com            *
*								Céline Moureau -  cemoureau@gmail.com                  *
*	For: High Performance Scientific course at ULiège, 2018-19                         *
*	Project 2                                                                          *
*                             														   *
*	Originally uploaded to: https://github.com/Cobalt1911                              *
*	Under GNU General Public License 11/2018                                           *
=======================================================================================*/

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include "algorithms.h"
/*=====================================================================================
* Contains all the necessary functions to handle BSR matrices and CSR vectors
=====================================================================================*/

// Structure
// CSR & BSR reference: https://medium.com/@jmaxg3/101-ways-to-store-a-sparse-matrix-c7f2bf15a229
typedef struct bsr_matrix bsr_matrix;
struct bsr_matrix{
	int nrows;
	int ncolumns;
	int block_size; // evenly divides nrows and ncolumns
	int nnzb; // # of non-zero *blocks*
	int n_elements_per_block;
	unsigned int *block_row_offsets;
	unsigned int *block_columns;
	double *values;
};

typedef struct csr_vector csr_vector;
struct csr_vector{
	int nrows;
	int nnzb; // # of non-zero values
	unsigned int *rows;
	double *values;
};

/*=====================================================================================*/

// Prototypes
void bsr_init(bsr_matrix *matrix, int nrows, int ncolumns, int block_size, int nnzb);
double bsr_get(bsr_matrix *matrix, int i, int j);
int natural_to_bsr(double *natural, bsr_matrix *matrix, int size, int block_size);
void bsr_free(bsr_matrix *matrix);
/*==============*/
void csr_vector_init(csr_vector *vector, double *natural, int nrows);
void csr_vector_init_empty(csr_vector *vector, int nrows);
double csr_vector_get(csr_vector *vector, int index);
void csr_vector_scale(csr_vector *vector, double scale);
int csr_vector_add(csr_vector *P, csr_vector *Q, csr_vector *R);
int csr_vector_set(csr_vector *vector, double value, int index);
double csr_vector_scalar(csr_vector *P, csr_vector *Q);
double csr_vector_norm(csr_vector *P);
void csr_vector_free(csr_vector *vector);
/*==============*/
int bsr_matrix_vector(bsr_matrix *matrix, csr_vector *vector, csr_vector *csr_result_vector);


/*=====================================================================================*/


/*============== BSR Matrix functions ===================*/
// Initialization function, sets all the variables and arrays from the structure
void bsr_init(bsr_matrix *matrix, int nrows, int ncolumns, int block_size, int nnzb){
	matrix->nrows = nrows;
	matrix->ncolumns = ncolumns;
	matrix->block_size = block_size;
	matrix->nnzb = nnzb;
	matrix->n_elements_per_block = block_size*block_size;

	matrix->block_row_offsets = malloc(sizeof(int) * (nrows / block_size + 1));
	matrix->block_columns = malloc(sizeof(int) * nnzb);
	matrix->values = malloc(sizeof(double) * nnzb*block_size*block_size);
}

// Fetches a value from the BSR matrix
// TODO : Simplify this once it's 100% working
double bsr_get(bsr_matrix *matrix, int i_, int j_) {
	// Find the corresponding block
	int i = floor(i_/matrix->block_size); // Block-row index
	int j = floor(j_/matrix->block_size); // Block-column index

	// Use CSR look-up algorithm to locate the block
	bool found = false;
	int block_row_offset = matrix->block_row_offsets[i];
	int step = 0;
	while (step < matrix->block_row_offsets[i+1] - block_row_offset){
		if (matrix->block_columns[block_row_offset+step] != j) ++step;
		else{
			found = true;
			break;
		}
	}
	if (!found) return 0; // If the block is not present in the CSR representation, it's a 0
	// At this point, the block has been found
	int offset = (block_row_offset+step)*matrix->n_elements_per_block;
	int i_block_start = matrix->block_size*i;
	int j_block_start = matrix->block_size*j;
	i_ -= i_block_start;
	j_ -= j_block_start;

offset += i_*matrix->block_size + j_;
return matrix->values[offset];
}

// Takes a matrix written as a 1D array and stores it as a bsr matrix!
int natural_to_bsr(double *natural, bsr_matrix *matrix, int size, int block_size) {

	if (size%block_size != 0){ // The block size is incompatible with the matrix size
		printf("!!! The block size is incompatible with the matrix size.\n");
		return -1;
	}
	unsigned int b_size = size/block_size;
	unsigned int *temp_block_row_offsets = malloc(sizeof(int) * (b_size + 1));
	unsigned int *temp_block_columns = malloc(sizeof(int) * (int)(b_size*b_size));
	double *temp_values = malloc(sizeof(double) * size*size);
	long value_index = 0, block_count = 0;
	char *block_matrix = malloc(sizeof(char)*b_size*b_size);
	int temp_row_index = 0;
	int temp_col_index = 0;

	/*======== Creation of the CSR matrix for the block matrix =========*/
	temp_block_row_offsets[0] = 0;
	// Browse the natural matrix, block by block
	for (int j=0 ; j<size ; j+=block_size){
		for (int i=0 ; i<size ; i+=block_size){
			// Looking at a sub-matrix of size block_size x block_size,
			// check if there is at least one non-zero value
			bool empty = true;
			for (int l=0 ; l<block_size ; ++l){
				for (int k=0 ; k<block_size ; ++k){
					if (natural[i+k + size*(j+l)] != 0){
						empty = false;
						break;
					}
				}
				if (!empty)break;
			}

			if (!empty){
				block_matrix[i/block_size + (j/block_size)*b_size] = 1;
				// Non-empty case, append temp_values
				for (int l=0 ; l<block_size ; ++l){
					for (int k=0 ; k<block_size ; ++k){
						temp_values[value_index] = natural[i+k + size*(j+l)];
						++value_index;
					}
				}
				temp_block_columns[temp_col_index] = i/block_size;
				++temp_col_index;
				++block_count;
			}
			else{
				// Empty case
				block_matrix[i/block_size + (j/block_size)*b_size] = 0;
			}
		}
		++temp_row_index;
		temp_block_row_offsets[temp_row_index] = block_count;
	}
	temp_block_row_offsets[temp_row_index+1] = block_count;

	/*======== Allocation of the BSR matrix =========*/
	// Give the values to the receiving bsr matrix
	bsr_init(matrix, size, size, block_size, block_count);
	for (int i=0 ; i<block_count*block_size*block_size ; ++i){
		matrix->values[i] = temp_values[i];
	}
	for (int i = 0; i < b_size+1; ++i){
		matrix->block_row_offsets[i] = temp_block_row_offsets[i];
	}
	printf("\n");
	for (int i = 0; i < block_count; ++i){
		matrix->block_columns[i] = temp_block_columns[i];
	}
	
	free(temp_block_row_offsets);
	free(temp_block_columns);
	free(temp_values);
	free(block_matrix);
}

// Frees the memory, fly away !
void bsr_free(bsr_matrix *matrix) {
	free(matrix->block_row_offsets);
	free(matrix->block_columns);
	free(matrix->values);
}


/*============== CSR Vector functions ===================*/

// Initialization function, sets all the variables and arrays from the structure
void csr_vector_init(csr_vector *vector, double *natural, int nrows){
	vector->nrows = nrows;
	int temp_nnzb = 0;
	for (int i = 0; i < nrows; ++i){
		if (natural[i] != 0)++temp_nnzb;
	}
	vector->nnzb = temp_nnzb;
	vector->rows = malloc(sizeof(int) * temp_nnzb);
	vector->values = malloc(sizeof(double) * temp_nnzb);
	int index = 0;
	for (int i = 0; i < nrows; ++i){
		if (natural[i] != 0){
			vector->rows[index] = i;
			vector->values[index] = natural[i];
			++index;
		}
	}
}

// Initialization function, sets all the variables and arrays from the structure
void csr_vector_init_empty(csr_vector *vector, int nrows){
	vector->nrows = nrows;
	vector->nnzb = 0;
	vector->rows = malloc(sizeof(int)*0);
	vector->values = malloc(sizeof(double)*0);
}

// Returns the value of a vector at a given index
double csr_vector_get(csr_vector *vector, int index){
	if (index >= vector->nrows){
		printf("!!! Index out of bounds.\n");
		return -1;
	}

	int looking_index = findIndex(vector->rows, index, vector->nnzb);
	if (looking_index > 0) return vector->values[looking_index];
	return 0;
}

// Sets a value at a given index of a CSR vector
// TODO : do it efficiently !
int csr_vector_set(csr_vector *vector, double value, int index){
	if (index >= vector->nrows){
		printf("!!! Index out of bounds.\n");
		return -1;
	}

	// Does the key insertion in rows and gives the index at which the value was inserted
	int *new_rows = malloc(sizeof(int)*(vector->nnzb+1));
	bool isPresent = false;
	int target_index = vector->nnzb;
	int new_nnzb;
	for (int i = 0; i < vector->nnzb; ++i){ // This part needs to be optimized
		if (vector->rows[i] == index){
			isPresent = true;
			target_index = i;
			break;
		}
		if (vector->rows[i] > index){
			target_index = i;
			break;
		}
	}

	// Case where the key is not already in the list
	if (!isPresent && value != 0){
		new_nnzb = vector->nnzb+1;
		double *new_values = malloc(sizeof(double)*new_nnzb);
		// Create new lists with the correct values
		for (int i = 0; i < new_nnzb; ++i){
			if (i == target_index){
				new_rows[i] = index; // Insertion
				new_values[i] = value;
			}
			else if(i > target_index){
				new_rows[i] = vector->rows[i-1];
				new_values[i] = vector->values[i-1];
			}
			else {
				new_rows[i] = vector->rows[i];
				new_values[i] = vector->values[i];
			}
		}

		// Setting the new values
		int *error_rows = realloc(vector->rows, new_nnzb*sizeof(int));
		int *error_values = realloc(vector->values, new_nnzb*sizeof(double));
		if (error_rows == NULL){
			printf("Error in rows allocation.\n");
			exit(0);
		}
		if (error_values == NULL){
			printf("Error in values allocation.\n");
			exit(0);
		}

		vector->nnzb = new_nnzb;
		vector->rows = error_rows;
		vector->values = error_values;
		for (int i = 0; i < new_nnzb; ++i){
			vector->rows[i] = new_rows[i];
			vector->values[i] = new_values[i];
		}
		free(new_values);
	}
	free(new_rows);

	// Case where the key is present, only updating the value
	if (isPresent){
		vector->values[target_index] = value;
	}
}

// Scales a CSR vector
void csr_vector_scale(csr_vector *vector, double scaling_factor){
	for (int i = 0; i < vector->nnzb; ++i){
		vector->values[i] *= scaling_factor;
	}
}

// Sums two CSR vectors and stores the result in a third vector
int csr_vector_add(csr_vector *P, csr_vector *Q, csr_vector *R){
	if (P->nrows != Q->nrows){
		printf("!!! Vector dimensions mismatch.\n");
		return -1;
	}

	double *result = malloc(sizeof(double)*P->nrows);
	int *adding_list = malloc(sizeof(int)*P->nrows);
	int max_nnzb = 0;
	if (P->nnzb > Q->nnzb)max_nnzb = P->nnzb;
	else max_nnzb = Q->nnzb;

	// Build the list of rows with non-zero values
	int new_nnzb = merge_sorted_lists(P->rows, Q->rows, adding_list, P->nnzb, Q->nnzb);

	// Compute the sum
	for (int i = 0; i < new_nnzb; ++i){
		result[i] = csr_vector_get(P, adding_list[i]) + csr_vector_get(Q, adding_list[i]);
	}

	// Do a half manual initilization, as the rows and values are already computed
	double empty[] = {0};
	csr_vector_init(R, empty, P->nrows);
	R->nnzb = new_nnzb;
	for (int i = 0; i < new_nnzb; ++i){
		R->rows[i] = adding_list[i];
		R->values[i] = result[i];
	}
	free(result);
	free(adding_list);
}

// Does a scalar product between two CSR vectors
double csr_vector_scalar(csr_vector *P, csr_vector *Q){
	double result = 0;
	int index_Q = 0;

	// Finding the intersection between the rows (non-zero values)
	for (int i = 0; i < P->nnzb; ++i){
		while (Q->rows[index_Q] < P->rows[i]){
			++index_Q;
		}
		if (Q->rows[index_Q] == P->rows[i]){
			result += P->values[i]*Q->values[index_Q];
		}
	}
	return result;
}

// Computes the Euclidian norm of a CSR vector
double csr_vector_norm(csr_vector *P){
	return sqrt(csr_vector_scalar(P,P));
}

// Frees the memory
void csr_vector_free(csr_vector *vector){
	free(vector->rows);
	free(vector->values);
}


/*============== BSR Matrix & CSR Vector functions ===================*/
// Does a BSR matrix/vector product
int bsr_matrix_vector(bsr_matrix *matrix, csr_vector *vector, csr_vector *csr_result_vector){
	if (vector->nrows != matrix->nrows){
		printf("!!! Matrix and vector dimensions mismatch.\n");
		return -1;
	}

	double *row_vector = malloc(sizeof(double)*matrix->nrows);
	double *result_vector = malloc(sizeof(double)*matrix->nrows);

	for (int i = 0; i < matrix->nrows; ++i){
		// Extract vectors by row from the matrix
		for (int j = 0; j<matrix->nrows; ++j){
			row_vector[j] = bsr_get(matrix, i, j);
		}
		// Convert it to CSR vector
		csr_vector csr_row_vector;
		csr_vector_init(&csr_row_vector, row_vector, matrix->nrows);
		// Do the scalar product
		result_vector[i] = csr_vector_scalar(vector, &csr_row_vector);
	}
	csr_vector_init(csr_result_vector, result_vector, matrix->nrows);
	free(row_vector);
	free(result_vector);
}
