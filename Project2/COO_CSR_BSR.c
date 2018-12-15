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
#include "COO_CSR_BSR.h"

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
void convert_natural_to_bsr(double *natural, bsr_matrix *matrix, int size, int block_size) {

	if (size%block_size != 0){ // The block size is incompatible with the matrix size
		printf("!!! The block size is incompatible with the matrix size.\n");
		exit(1);
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


/*====================== COO Matrix ======================================================*/

// Initialization function, sets all the variables and arrays from the structure
void coo_matrix_init_empty(coo_matrix *matrix, int nrows, int ncolumns){
	matrix->nrows = nrows;
	matrix->ncolumns = ncolumns;
	matrix->nnzb = 0;
	matrix->rows = malloc(sizeof(int));
	matrix->columns = malloc(sizeof(int));
	matrix->values = malloc(sizeof(double));
	matrix->mem_size = 1;
}

// Fetches a value from the COO matrix
double coo_matrix_get(coo_matrix *matrix, int i, int j){
	int row_index = 0;
	for (int k = 0; k < matrix->nnzb; ++k){
		if (matrix->rows[k] == i && matrix->columns[k] == j){
			//printf("k: %d\n", k);
			return matrix->values[k];
		}
	}
	return 0;
}

// Takes a matrix written as a 1D array and stores it as a COO matrix!
void convert_natural_to_coo(double *natural, coo_matrix *matrix, int size){
	if (size != matrix->nrows*matrix->ncolumns){
		printf("!!! Matrix dimensions mismatch.\n");
		exit(1);
	}
	matrix->nnzb = -1;
	for (int i = 0; i < size; ++i){
		if (natural[i] != 0){
			++ matrix->nnzb;
			coo_matrix_resize(matrix, matrix->nnzb);
			matrix->rows[matrix->nnzb] = floor(i/matrix->nrows);
			matrix->columns[matrix->nnzb] = i - matrix->nrows*floor(i/matrix->nrows);
			matrix->values[matrix->nnzb] = natural[i];
		}
	}
}

// Resizes the memory allocated to a csr_vector type
void coo_matrix_resize(coo_matrix *matrix, int new_nnzb){
	// Check if new space is needed
	if (matrix->mem_size < new_nnzb){
		matrix->mem_size *= 2;
		unsigned int *ptr_rows = realloc(matrix->rows, matrix->mem_size*sizeof(int));
		unsigned int *ptr_columns = realloc(matrix->columns, matrix->mem_size*sizeof(int));
		double *ptr_values = realloc(matrix->values, matrix->mem_size*sizeof(double));
		if (ptr_rows == NULL || ptr_columns == NULL || ptr_values == NULL){
			printf("\nMem ERR0R !\n");
			exit(1);
		}
		matrix->rows = ptr_rows;
		matrix->columns = ptr_columns;
		matrix->values = ptr_values;
	}
}

// Scales a COO matrix
void coo_matrix_scale(coo_matrix *matrix, double scale){
	for (int i = 0; i < matrix->nnzb; ++i){
		matrix->values[i] *= scale;
	}
}

// Frees the memory, fly away !
void coo_matrix_free(coo_matrix *matrix){
	free(matrix->rows);
	free(matrix->columns);
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
	vector->mem_size = temp_nnzb;
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
	vector->rows = malloc(sizeof(int));
	vector->values = malloc(sizeof(double));
	vector->mem_size = 1;
}

// Returns the value of a vector at a given index
double csr_vector_get(csr_vector *vector, int index){
	if (index >= vector->nrows){
		printf("!!! Index out of bounds.\n");
		return -1;
	}
	int looking_index = findIndex(vector->rows, index, vector->nnzb);
	if (looking_index >= 0) return vector->values[looking_index];
	return 0;
}

// Sets a value at a given index of a CSR vector
void csr_vector_set(csr_vector *vector, double value, int index){
	if (index >= vector->nrows){
		printf("!!! Index out of bounds.\n");
		exit(1);
	}

	// Does the key insertion in rows and gives the index at which the value was inserted
	int target_index;
	// Checks if a value already exists at this index
	target_index = findIndex(vector->rows, index, vector->nnzb);
	if (target_index != -1){ // Case where the key is present, only updating the value
		vector->values[target_index] = value;
	}
	else if (value != 0){
		int new_nnzb = vector->nnzb+1;
		// Increase size of rows and values
		csr_vector_resize(vector, new_nnzb);
		vector->nnzb = new_nnzb;
		// Insert the key in rows
		target_index = insertKey_int(vector->rows, index, vector->nnzb);
		// Insert value in values at the same index
		insertKeyAt_double(vector->values, value, new_nnzb, target_index);
	}
}


// Scales a CSR vector
void csr_vector_scale(csr_vector *vector, double scaling_factor){
	for (int i = 0; i < vector->nnzb; ++i){
		vector->values[i] *= scaling_factor;
	}
}

// Sums two CSR vectors and stores the result in a third vector
void csr_vector_add(csr_vector *P, csr_vector *Q, csr_vector *R){
	if (P->nrows != Q->nrows){
		printf("!!! Vector dimensions mismatch.\n");
		exit(1);
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
	if (P->nrows != Q->nrows){
		printf("!!! Vector dimensions mismatch.\n");
		exit(1);
	}
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

// Resizes the memory allocated to a csr_vector type
void csr_vector_resize(csr_vector *vector, int new_nnzb){
	// Check if new space is needed
	if (vector->mem_size < new_nnzb){
		if (vector->mem_size == 0) vector->mem_size = 1;
		vector->mem_size *= 2;
		unsigned int *ptr_rows = realloc(vector->rows, vector->mem_size*sizeof(int));
		double *ptr_values = realloc(vector->values, vector->mem_size*sizeof(double));
		// Should check if !NULL
		if (ptr_rows == NULL || ptr_values == NULL){
			printf("\nMem ERR0R !\n");
			exit(1);
		}
		vector->rows = ptr_rows;
		vector->values = ptr_values;
	}
}

// Copies the data from Q into P ( P = Q )
void csr_vector_equals(csr_vector *P, csr_vector *Q){
	if (P->nrows != Q->nrows){
		printf("!!! Vector dimensions mismatch.\n");
		exit(1);
	}
	csr_vector_resize(P, Q->nnzb);
	P->nnzb = Q->nnzb;
	for (int i = 0; i < P->nnzb; ++i){
		P->rows[i] = Q->rows[i];
		P->values[i] = Q->values[i];
	}
}

// Frees the memory
void csr_vector_free(csr_vector *vector){
	free(vector->rows);
	free(vector->values);
}


/*============== Natural vector functions ===================*/

// Initialization function
void nat_vector_init(nat_vector *vector, int nrows){
	vector->nrows = nrows;
	vector->values = calloc(nrows, sizeof(double));
}

// Scales a natural vector
void nat_vector_scale(nat_vector *vector, double scaling_factor){
	for (int i = 0; i < vector->nrows; ++i){
		vector->values[i] *= scaling_factor;
	}
}

// Sums two natural vectors and stores the result in a third vector
void nat_vector_add(nat_vector *P, nat_vector *Q, nat_vector *R){
	if (P->nrows != Q->nrows || P->nrows != R->nrows){
		printf("!!! Vector dimensions mismatch.\n");
		exit(1);
	}
	for (int i = 0; i < P->nrows; ++i){
		R->values[i] = P->values[i] + Q->values[i];
	}
}

// Does a scalar product between two natural vectors
double nat_vector_scalar(nat_vector *P, nat_vector *Q){
	if (P->nrows != Q->nrows){
		printf("!!! Vector dimensions mismatch.\n");
		return -1;
	}
	double result = 0;
	for (int i = 0; i < P->nrows; ++i){
		result += P->values[i]*Q->values[i];
	}
	return result;
}

// Computes the Euclidian norm of a natural vector
double nat_vector_norm(nat_vector *P){
	return sqrt(nat_vector_scalar(P,P));
}

// Copies the data from Q into P ( P = Q )
void nat_vector_equals(nat_vector *P, nat_vector *Q){
	if (P->nrows != Q->nrows){
		printf("!!! Vector dimensions mismatch.\n");
		exit(1);
	}
	for (int i = 0; i < P->nrows; ++i){
		P->values[i] = Q->values[i];
	}
}

// Frees the memory
void nat_vector_free(nat_vector *vector){
	free(vector->values);
}

/*============== Matrix & Vector functions ===================*/
// Does a BSR matrix/vector product
void bsr_matrix_csr_vector(bsr_matrix *matrix, csr_vector *vector, csr_vector *csr_result_vector){
	if (vector->nrows != matrix->nrows){
		printf("!!! Matrix and vector dimensions mismatch.\n");
		exit(1);
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

// Does a COO matrix/ natural vector product
void coo_matrix_nat_vector(coo_matrix *matrix, nat_vector *vector, nat_vector *result){
	if (vector->nrows != matrix->ncolumns){
		printf("!!! Matrix and vector dimensions mismatch.\n");
		exit(1);
	}
	for (int i = 0; i < matrix->nnzb; ++i){
		result->values[matrix->rows[i]] += matrix->values[i]*vector->values[matrix->columns[i]];
	}
}

/*============== Direct implementation ===================*/
// Computes the square of the Euclidian norm of a natural vector
double SquaredNorm(double* vect, size_t sizeVect)
{
  double normVal=0;
	#pragma omp for
  for(size_t i = 0; i< sizeVect; i++)
	normVal+=vect[i]*vect[i];
  return normVal;
}

// Does a scalar product between two natural vectors
double VectorProduct(double* vect1, double* vect2, size_t sizeVect)
{
  double prod=0;
  for(size_t i = 0; i< sizeVect; i++)
	prod+=vect1[i]*vect2[i];
  return prod;
}

// Sums two natural vectors and stores the result in a third vector
void SumVect(double* vectToStore, double* vect1, double* vect2, double multVal, size_t sizeVect)
{
  for(size_t i = 0; i< sizeVect; i++)
	vectToStore[i] = vect1[i] + vect2[i]*multVal;
}

// Pre-computes the application of A to a vector p
void Ap(double* p, double* Apresult, size_t nodeX, size_t nodeY, size_t thicknessMPI, double h, double m, double vx, double vy, double vz, double D, int rank, int world_size)
{

  for(size_t index = 0; index< nodeX*nodeY*thicknessMPI; index++)
  {
	bool onZBoundary = false;
	int k = floor(index/(nodeX*nodeY));
	int j = floor((index-k*nodeX*nodeY)/nodeX);
	int i = index - k * nodeX * nodeY - j * nodeX;
	double Ddivh2 = D/(h*h);
	if (rank == 0) onZBoundary = (index<nodeX*nodeY);
	if (rank == world_size-1) onZBoundary = (index>=(thicknessMPI-1)*nodeX*nodeX);
	if(!(i==0 || i == nodeX-1 || j==0 || j == nodeY-1 || onZBoundary))
	{
	  k+=1;
	  Apresult[index] = p[i+j*nodeX+k*nodeX*nodeY]*(1/m + 6*Ddivh2) + p[i-1+j*nodeX+k*nodeX*nodeY]*(-vx/(2*h) - Ddivh2) + p[i+1+j*nodeX+k*nodeX*nodeY]*(vx/(2*h) - Ddivh2) + p[i+(j-1)*nodeX+k*nodeX*nodeY]*(-vy/(2*h) - Ddivh2) + p[i+(j+1)*nodeX+k*nodeX*nodeY]*(vy/(2*h) - Ddivh2) + p[i+j*nodeX+(k-1)*nodeX*nodeY]*(-vz/(2*h) -Ddivh2) + p[i+j*nodeX+(k+1)*nodeX*nodeY]*(vz/(2*h) - Ddivh2);
	  Apresult[index] *= m;
	}
  }

}