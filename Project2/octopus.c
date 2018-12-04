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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int explicit_solver(int argc, char *argv[]);
int implicit_solver(int argc, char *argv[]);

int main(int argc, char *argv[]){
	if (argc != 3) {
		printf("Wrong arguments\n");
		printf("Please use the function with ./exe param.dat\n");
		return 1;
	}

	if (!strcmp(argv[2], "0")){
		explicit_solver(argc, argv);
	}
	if (!strcmp(argv[2], "1")){
		implicit_solver(argc, argv);		
	}
	return 0;
}