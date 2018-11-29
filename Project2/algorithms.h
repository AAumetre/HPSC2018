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

// Prototypes
int merge_sorted_lists(int list_A[], int list_B[], int list_C[], int size_A, int size_B);
int * getCommListSlices(unsigned int size);
double myRound(double value);
int findIndex(int list[], int key, int size);
int * shareWorkload(int problemSize, int nNodes);
int insertKey_int(int list[], int key, int size);
void insertKeyAt_double(double list[], double key, int size, int target_index);