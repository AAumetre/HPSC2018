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
int sorted_list_insertion_int(int target[], int key, int target_size);

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

// Inserts a key into a sorted list, returns the index at which the key was inserted
// !!!! If the key is already present, this function does not add an other occurence
// TODO : do it efficiently !
int sorted_list_insertion_int(int *target, int key, int target_size){
	int *resulting_list = malloc(sizeof(int)*(target_size+1));
	bool isPresent = false;
	int index;
	for (int i = 0; i < target_size; ++i){
		if (target[i] == key){
			isPresent = true;
			index = i;
			break;
		}
		if (target[i] > key){
			index = i;
			break;
		}
	}

	// Case where the key is not already in the list	
	if (!isPresent){
		for (int i = 0; i < target_size+1; ++i){
			if (i == index) resulting_list[i] = key; // Insertion
			if (i > index) resulting_list[i+1] = target[i];
			else resulting_list[i] = target[i]; 
		}
		realloc(target, sizeof(int)*(target_size+1)); // Resizing
		memcpy(target, resulting_list, target_size+1);
	}

	return index;
}

// Inserts a key at a given location, in a sorted list
void sorted_list_insertion_index_int(int target[], int key, int index, int target_size){

}