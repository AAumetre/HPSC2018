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
/*=====================================================================================
* Contains all the necessary functions to handle BSR matrices and CSR vectors
=====================================================================================*/

// CSR & BSR reference: https://medium.com/@jmaxg3/101-ways-to-store-a-sparse-matrix-c7f2bf15a229

// Structure
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

typedef struct coo_matrix coo_matrix;
struct coo_matrix{
	int nrows;
	int ncolumns;
	int nnzb;
	unsigned int *rows;
	unsigned int *columns;
	double *values;
	int mem_size;
};

typedef struct csr_vector csr_vector;
struct csr_vector{
	int nrows;
	int nnzb; // # of non-zero values
	unsigned int *rows;
	double *values;
	int mem_size; // in # of elements, not bytes
};

typedef struct nat_vector nat_vector;
struct nat_vector{
	int nrows;
	double *values;
};

/*====================== BSR Matrix ======================================================*/

// Initialization function, sets all the variables and arrays from the structure
void bsr_init(bsr_matrix *matrix, int nrows, int ncolumns, int block_size, int nnzb);

// Fetches a value from the BSR matrix
double bsr_get(bsr_matrix *matrix, int i, int j);

// Takes a matrix written as a 1D array and stores it as a BSR matrix!
void convert_natural_to_bsr(double *natural, bsr_matrix *matrix, int size, int block_size);

// Frees the memory, fly away !
void bsr_free(bsr_matrix *matrix);

/*====================== COO Matrix ======================================================*/

// Initialization function, sets all the variables and arrays from the structure
void coo_matrix_init_empty(coo_matrix *matrix, int nrows, int ncolumns);

// Fetches a value from the COO matrix
double coo_matrix_get(coo_matrix *matrix, int i, int j);

// Takes a matrix written as a 1D array and stores it as a COO matrix!
void convert_natural_to_coo(double *natural, coo_matrix *matrix, int size);

// Resizes the memory allocated to a csr_vector type
void coo_matrix_resize(coo_matrix *matrix, int new_nnzb);

// Frees the memory, fly away !
void coo_matrix_free(coo_matrix *matrix);

/*====================== CSR Vector ======================================================*/
// Initialization function, sets all the variables and arrays from the structure
void csr_vector_init(csr_vector *vector, double *natural, int nrows);

// Initialization function, sets all the variables and arrays from the structure
void csr_vector_init_empty(csr_vector *vector, int nrows);

// Returns the value of a vector at a given index
double csr_vector_get(csr_vector *vector, int index);

// Sets a value at a given index of a CSR vector
void csr_vector_set(csr_vector *vector, double value, int index);

// Scales a CSR vector
void csr_vector_scale(csr_vector *vector, double scale);

// Sums two CSR vectors and stores the result in a third vector
void csr_vector_add(csr_vector *P, csr_vector *Q, csr_vector *R);

// Does a scalar product between two CSR vectors
double csr_vector_scalar(csr_vector *P, csr_vector *Q);

// Computes the Euclidian norm of a CSR vector
double csr_vector_norm(csr_vector *P);

// Resizes the memory allocated to a csr_vector type
void csr_vector_resize(csr_vector *vector, int new_nnzb);

// Copies the data from Q into P ( P = Q )
void csr_vector_equals(csr_vector *P, csr_vector *Q);

// Scales a COO matrix
void coo_matrix_scale(coo_matrix *matrix, double scale);

// Frees the memory
void csr_vector_free(csr_vector *vector);

/*====================== Natural Vector ==================================================*/
// Initialization function
void nat_vector_init(nat_vector *vector, int nrows);

// Scales a natural vector
void nat_vector_scale(nat_vector *vector, double scale);

// Sums two natural vectors and stores the result in a third vector
void nat_vector_add(nat_vector *P, nat_vector *Q, nat_vector *R);

// Does a scalar product between two natural vectors
double nat_vector_scalar(nat_vector *P, nat_vector *Q);

// Computes the Euclidian norm of a natural vector
double nat_vector_norm(nat_vector *P);

// Copies the data from Q into P ( P = Q )
void nat_vector_equals(nat_vector *P, nat_vector *Q);

// Frees the memory
void nat_vector_free(nat_vector *vector);

/*====================== BSR CSR ========================================================*/
// Does a BSR matrix/ CSR vector product
void bsr_matrix_csr_vector(bsr_matrix *matrix, csr_vector *vector, csr_vector *csr_result_vector);

// Does a BSR matrix/ CSR vector product
void coo_matrix_nat_vector(coo_matrix *matrix, nat_vector *vector, nat_vector *result);
