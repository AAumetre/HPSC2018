
/*=======================================================================================
*	This code was written by: 
*								Antonin Aumètre - antonin.aumetre@gmail.com
*								Céline Moureau -  cemoureau@gmail.com
*	For: High Performance Scientific course at ULiège, 2018-19
*	Project 2
*
*	Under GNU General Public License 11/2018
=======================================================================================*/

// Structure
// BSR https://medium.com/@jmaxg3/101-ways-to-store-a-sparse-matrix-c7f2bf15a229
// Implementation
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

// Prototypes
void bsr_init(bsr_matrix *matrix, int nrows, int ncolumns, int block_size, int nnzb);
double bsr_get(bsr_matrix *matrix, int i, int j);
void bsr_set(bsr_matrix *matrix, int i, int j, double value);
void bsr_free(bsr_matrix *matrix);

// Implementation
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
// TODO : Find a way to tell wether a block exists or not (is empty or not)
double bsr_get(bsr_matrix *matrix, int i_, int j_) {
   // Find the corresponding block
  int i = floor(i_/matrix->block_size); // Block-row index
  int j = floor(j_/matrix->block_size); // Block-column index

  //printf("(i,j) = (%d,%d)\n", i, j);

  // Use CSR look-up algorithm to locate the block
  bool found = false;
  int block_row_offset = matrix->block_row_offsets[i];
  //printf("Block row offset = %d\n", block_row_offset);
  int step = 0;
  while (step < matrix->block_row_offsets[i+1] - block_row_offset){
  	if (matrix->block_columns[block_row_offset+step] != j) ++step;
  	else{
		found = true;
		break;
	}
  }
  if (!found) return 0; // If the block is not present in the CSR representation, it's a 0
  //printf("Step value : %d\n", step);

  // At this point, the block has been found
  int offset = (block_row_offset+step)*matrix->n_elements_per_block;
  //printf("Offset: %d\n", offset);
  int i_block_start = matrix->block_size*i;
  int j_block_start = matrix->block_size*j;
  //printf("(i,j)_start = (%d,%d)\n", i_block_start, j_block_start);
  i_ -= i_block_start;
  j_ -= j_block_start;

  offset += i_*matrix->block_size + j_;
  //printf("Offset: %d\n", offset);
  return matrix->values[offset];
}

// Takes a matrix written as a 1D array and stores it as a bsr matrix
void natural_to_bsr(double *natural, bsr_matrix *matrix, int size, int block_size) {

	unsigned int *temp_block_row_offsets = malloc(sizeof(int) * (size / block_size + 1));
	unsigned int *temp_block_columns = malloc(sizeof(int) * (int)(size*size/block_size/block_size));
	double *temp_values = malloc(sizeof(double) * size*size);
	long value_index = 0, block_count = 0;

	for (int i=0 ; i<size ; i+=block_size){
		for (int j=0 ; j<size ; j+=block_size){
			// Looking at a sub-matrix of size block_size x block_size

			// Check if there is at least one non-zero value
			bool empty = true;
			for (int k=0 ; k<block_size ; ++k){
				for (int l=0 ; l<block_size ; ++l){
					if (natural[i*block_size+k + size*(j*block_size+l)] != 0){
						empty = false;
						break;
					} 
				}
				if (!empty) break;
			}

			if (!empty){
				// Non-empty case, append temp_values
				for (int k=0 ; k<block_size ; ++k){
					for (int l=0 ; l<block_size ; ++l){
						temp_values[value_index] = natural[i*block_size+k + size*(j*block_size+l)];
						++value_index;
					}
				}
				++block_count;
			}
			else{// Empty case

			}
		}
	}
	bsr_init(matrix, size, size, block_size, block_count);
	for (int i=0 ; i<block_count*block_size*block_size ; ++i){
		matrix->values[i] = temp_values[i];
		//printf("%f ", temp_values[i]);
	}
}

// Frees the memory, fly away !
void bsr_free(bsr_matrix *matrix) {
	free(matrix->block_row_offsets);
	free(matrix->block_columns);
	free(matrix->values);
}