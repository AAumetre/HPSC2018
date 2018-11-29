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
* Contains general purpose algorithms, can be used anywhere
=====================================================================================*/

// Merges two sorted lists and stores the result in a third one. Removes multiple occurences.
// Return the length of the merged list.
int merge_sorted_lists(int list_A[], int list_B[], int list_C[], int size_A, int size_B);

// Creates a list of the directionnal communications to be done
int * getCommListSlices(unsigned int size);

// Natural round-off
double myRound(double value);

// Finds the index at which a key is present in a sorted list and returns it. Returns -1 if the key was not found.
int findIndex(int list[], int key, int size);

// Function that shares the workload
int * shareWorkload(int problemSize, int nNodes);

// Inserts a key in a list, returns the index, does not prevent redundance.
int insertKey_int(int list[], int key, int size);

// Inserts a key in a list, returns the index, does not prevent redundance.
void insertKeyAt_double(double list[], double key, int size, int target_index);