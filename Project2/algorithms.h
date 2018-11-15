/*=======================================================================================
*	This code was written by:
*								Antonin Aumètre - antonin.aumetre@gmail.com
*								Céline Moureau -  cemoureau@gmail.com
*	For: High Performance Scientific course at ULiège, 2018-19
*	Project 2
*
*	Under GNU General Public License 11/2018
=======================================================================================*/

/*=====================================================================================
* Contains general purpose algorithms, can be used anywhere
=====================================================================================*/

// Prototypes
int merge_sorted_lists(int list_A[], int list_B[], int list_C[], int size_A, int size_B);
int * getCommListSlices(unsigned int size);

/*=====================================================================================*/
// Implementation

// Merges two sorted lists and stores the result in a third one. Removes multiple occurences.
// Return the length of the merged list.
int merge_sorted_lists(int list_A[], int list_B[], int list_C[], int size_A, int size_B) {
	int i=0, j=0, k=0;

      // Do the first iteration out of the loop to manage the case k=0
	if (list_A[i] <= list_B[j]) {
		list_C[k] = list_A[i];
		++i;
	}
	else {
		list_C[k] = list_B[j];
		++j;
	}
	++k;
    // Main loop
	while (i < size_A && j < size_B) {
		if (list_A[i] <= list_B[j]) {
			list_C[k] = list_A[i];
			++i;
			++k;
		}
		else if (list_B[j] != list_C[k-1]) {
			list_C[k] = list_B[j];
			++j;
			++k;
		}
		else {
			++j;
		}
	}
    // If one list is exhausted, fill C with the remains of the other
	if (i < size_A) {
		for (int p = i; p < size_A; ++p) {
			if (list_A[p] != list_C[k-1]){
				list_C[k] = list_A[p];
				++k;
			}
		}
	}
	else {
		for (int p = j; p < size_B; ++p) {
			if (list_B[p] != list_C[k-1]){
				list_C[k] = list_B[p];
				++k;
			}
		}
	}
	return k;
}

// Creates a list of the directionnal communications to be done
int * getCommListSlices(unsigned int size){
	int *list = malloc(sizeof(int)*4*(size-1));

	int max_value_begin;
	if (size%2 == 0) max_value_begin = size;
	else max_value_begin = size - 1;

	for (int i = 0; i < max_value_begin; ++i){
		list[i] = i;
	}
	for (int i = 0; i < max_value_begin; ++i){ // Permutations
		if (i%2 == 0) list[max_value_begin+i] = list[i+1];
		else list[max_value_begin+i] = list[i-1];
	}

	int max_value;
	if (size%2 == 0) max_value = size - 1;
	else max_value = size;
	for (int i = 1; i < max_value; ++i){
		list[2*max_value_begin+i-1] = i;
	}
	for (int i = 0; i < max_value-1; ++i){ // Permutations
		if (i%2 ==0) list[2*max_value_begin+max_value-1+i] = list[2*max_value_begin+1+i];
		else list[2*max_value_begin+max_value-1+i] = list[2*max_value_begin-1+i];
	}

	for (int i = 0; i < 4*(size-1); ++i)
		print("%d ", list[i]);
	return list;
}
