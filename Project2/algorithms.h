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