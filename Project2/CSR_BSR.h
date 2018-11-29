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
	int mem_size; // in # of elements, not bytes
};

typedef struct nat_vector nat_vector;
struct nat_vector{
	int nrows;
	double *values;
};

/*=====================================================================================*/
// Prototypes

void bsr_init(bsr_matrix *matrix, int nrows, int ncolumns, int block_size, int nnzb);
double bsr_get(bsr_matrix *matrix, int i, int j);
int convert_natural_to_bsr(double *natural, bsr_matrix *matrix, int size, int block_size);
void bsr_free(bsr_matrix *matrix);
/*==============*/
void csr_vector_init(csr_vector *vector, double *natural, int nrows);
void csr_vector_init_empty(csr_vector *vector, int nrows);
double csr_vector_get(csr_vector *vector, int index);
void csr_vector_scale(csr_vector *vector, double scale);
int csr_vector_add(csr_vector *P, csr_vector *Q, csr_vector *R);
int csr_vector_set(csr_vector *vector, double value, int index);
void csr_vector_resize(csr_vector *vector, int new_nnzb);
double csr_vector_scalar(csr_vector *P, csr_vector *Q);
double csr_vector_norm(csr_vector *P);
int csr_vector_equals(csr_vector *P, csr_vector *Q);
void csr_vector_free(csr_vector *vector);
/*==============*/
void nat_vector_init(nat_vector *vector, int nrows);
void nat_vector_scale(nat_vector *vector, double scale);
int nat_vector_add(nat_vector *P, nat_vector *Q, nat_vector *R);
double nat_vector_scalar(nat_vector *P, nat_vector *Q);
double nat_vector_norm(nat_vector *P);
int nat_vector_equals(nat_vector *P, nat_vector *Q);
void nat_vector_free(nat_vector *vector);
/*==============*/
int bsr_matrix_vector(bsr_matrix *matrix, csr_vector *vector, csr_vector *csr_result_vector);
