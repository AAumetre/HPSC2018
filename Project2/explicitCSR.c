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
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

#include "CSR_BSR.h"
#include "fileIO.h"

#define initConcentration 1 //[g/m3]


int main(int argc, char *argv[])
{
	if (argc != 2) {
		printf("Wrong arguments\n");
		printf("Please use the function with ./exe param.dat\n");
		return 1;
	}

	// Initalizes MPI
	MPI_Init(NULL, NULL);
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	//Check out the time
	clock_t begin = clock();
	clock_t end = clock();
	double time_spent;

	// Declare variables for the main loop
	size_t iteration = 0;
	bool onBoundary = false;
	bool onZBoundary = false;
	bool valueOnBoundary = false;
	int *stopFlags = calloc(world_size, sizeof(int));
	int *stopFlagsFromOthers = calloc(world_size, sizeof(int));
	bool stopFlag = false;
	bool isIdle = false;

	// Retrieving data from the .dat file
	Param parameters = readDat(argv[1]);
	size_t nodeX = (int)(parameters.L/parameters.h) + 1;
	size_t nodeY = nodeX, nodeZ =nodeX;
	size_t stopTime = parameters.Tmax/parameters.m;

	if (rank == 0) {
		printf("We have %d nodes\n", world_size);
		printf("Number of slices: %zu\n", nodeZ);
		printf("=====================================\n");
	}
	MPI_Barrier(MPI_COMM_WORLD);

	// If there are more nodes than possible slices, the surplus nodes are set to idle
	int *ranks = calloc(world_size, sizeof(int));
	if (world_size > nodeZ){ // Change the value of world_size to that further calculations are still valid
		world_size = nodeZ;
		for (int i = 0; i < nodeZ; ++i){
			ranks[i] = i;
		}
		if (rank >= world_size){ // Idle this node
			printf("Node %d was set to idle.\n", rank);
			isIdle = true;
			stopFlag = true;
		}
	}
	else{
		for (int i = 0; i < world_size; ++i){
			ranks[i] = i;
		}
	}

	// If there are some idling nodes, exclude them from the new communicator
	MPI_Group world_group;
	MPI_Comm_group(MPI_COMM_WORLD, &world_group);
	MPI_Group sub_world_group;
	MPI_Group_incl(world_group, world_size, ranks, &sub_world_group);
	// Define the new communiator associated with the sub-group
	MPI_Comm SUB_COMM;
	MPI_Comm_create_group(MPI_COMM_WORLD, sub_world_group, 0, &SUB_COMM);

	// Assign each node its nomber of slices (thicknessMPI)
	int nbAdditionalSlices = 0;
	size_t thicknessMPI = 0;
	int *share = shareWorkload(nodeZ, world_size);
	if (rank == world_size-1 && world_size != 1){
		thicknessMPI = share[0];
		nbAdditionalSlices = share[0]-share[1];
	}
	else thicknessMPI = share[1];

	// Some prints
	if (rank == 0) printf("Vx, Vy, Vz: %.2f %.2f %.2f\n", parameters.vx, parameters.vy, parameters.vz);
	MPI_Barrier(MPI_COMM_WORLD);
	if (!isIdle) printf("Thickness = %zu, for rank %d\n", thicknessMPI, rank);
	MPI_Barrier(MPI_COMM_WORLD);
	if (rank == 0) printf("=====================================\n");

	// Get memory for the concentration values
	csr_vector concentration;
	csr_vector_init_empty(&concentration, nodeX*nodeY*thicknessMPI);
	csr_vector c_;
	csr_vector_init_empty(&c_, nodeX*nodeY*(thicknessMPI+2));
	
	// Give initial concentration value
	int middleSliceIndex = floor(nodeZ/2);
	int rankMiddle = floor(middleSliceIndex/(thicknessMPI-nbAdditionalSlices));
	if (rank == rankMiddle){
		int initValueIndex = nodeX*floor(nodeY/2) + floor(nodeX/2) + nodeX*nodeY*(middleSliceIndex-rank*(thicknessMPI-nbAdditionalSlices));
		csr_vector_set(&c_, initConcentration, initValueIndex); 
		printf("Initial value set on node #%d, at index %d\n", rank, initValueIndex);
	}
	MPI_Barrier(MPI_COMM_WORLD);

	/*================================================================================================
	#	Main loop
	================================================================================================*/
	if (rank == 0) printf("Iteration: ");
	while (iteration <= stopTime && !stopFlag)
	{
		if (rank == 0) printf("%ld ", iteration);
		// Search for boundaries
		size_t isXbound = 0;
		size_t index = 0;
		if (rank == 0) index = nodeX*nodeY;
		size_t stopIndex = nodeX*nodeY*thicknessMPI;
		if (rank == world_size-1) stopIndex -= nodeY*nodeX;

		// ============================== Compute internal values
		for(index=0; index<stopIndex; index++)
		{
			int stage = floor(index/(nodeX*nodeY));
			int inStage0 = index-stage*nodeX*nodeY;

			if (!(isXbound==0 || isXbound==nodeX-1 || inStage0<nodeX || inStage0>=nodeX*nodeY-nodeX))
			{//I am in the domain
				int k = floor(index/(nodeX*nodeY));
				int j = floor((index-k*nodeX*nodeY)/nodeX);
				int i = index - k * nodeX * nodeY - j * nodeX;

				int kbis = k+1;
				int jbis = floor((index+nodeX*nodeY-kbis*nodeX*nodeY)/nodeX);
				int ibis = index+nodeX*nodeY - kbis * nodeX * nodeY - jbis * nodeX;

				// Somehow, the csr_vector_get always returns 0.
				double cValue = 
					csr_vector_get(&c_, ibis+jbis*nodeX+kbis*nodeX*nodeY) +
					parameters.m * parameters.D * (csr_vector_get(&c_, ibis+1+jbis*nodeX+kbis*nodeX*nodeY)+csr_vector_get(&c_, ibis+(jbis+1)*nodeX+kbis*nodeX*nodeY)+
					csr_vector_get(&c_, ibis+jbis*nodeX+(kbis+1)*nodeX*nodeY)-6*csr_vector_get(&c_, ibis+jbis*nodeX+kbis*nodeX*nodeY)+
					csr_vector_get(&c_, ibis-1+jbis*nodeX+kbis*nodeX*nodeY)+csr_vector_get(&c_, ibis+(jbis-1)*nodeX+kbis*nodeX*nodeY)+
					csr_vector_get(&c_, ibis+jbis*nodeX+(kbis-1)*nodeX*nodeY))/pow(parameters.h,2) -
					parameters.m * parameters.vx * (csr_vector_get(&c_, ibis+1+jbis*nodeX+kbis*nodeX*nodeY)-csr_vector_get(&c_, ibis-1+jbis*nodeX+kbis*nodeX*nodeY))/(2*parameters.h) -
					parameters.m * parameters.vy * (csr_vector_get(&c_, ibis+(jbis+1)*nodeX+kbis*nodeX*nodeY)-csr_vector_get(&c_, ibis+(jbis-1)*nodeX+kbis*nodeX*nodeY))/(2*parameters.h) -
					parameters.m * parameters.vz * (csr_vector_get(&c_, ibis+jbis*nodeX+(kbis+1)*nodeX*nodeY)-csr_vector_get(&c_, ibis+jbis*nodeX+(kbis-1)*nodeX*nodeY))/(2*parameters.h);
				csr_vector_set(&concentration, cValue, i+j*nodeX+k*nodeX*nodeY); 

				// Check the stopping conditions on the boundaries of the domain
				if (rank == 0) onZBoundary = (index<2*nodeX*nodeY);
				if (rank == world_size-1) onZBoundary = (index>=(thicknessMPI-2)*nodeX*nodeX);
				onBoundary = ((isXbound == nodeX-2) || (isXbound == 1) || (inStage0 >= nodeX && inStage0 <= 2*nodeX-1) || (inStage0 >= nodeY*nodeY-2*nodeX));
				if ((onBoundary || onZBoundary)  && (csr_vector_get(&concentration, i+j*nodeX+k*nodeX*nodeY) > 1e-12)){
				  valueOnBoundary=true;
				}
			}

			isXbound++;
			if(isXbound==nodeX) isXbound = 0;
		}

	// Send your status to all the other nodes
	for (int i = 0; i < world_size; ++i){
		if(valueOnBoundary)stopFlags[i] = 1;
	}
	MPI_Allgather(stopFlags, 1, MPI_INT, stopFlagsFromOthers, 1, MPI_INT, SUB_COMM);
	for (int i = 0; i < world_size; ++i){
		if(stopFlagsFromOthers[i] == 1){
		stopFlag = true;
		break;
		}
	}

		// ============================== Send and receive neighboring values
		int *commList = getCommListSlices(world_size);
		for (int commIndex=0 ; commIndex<4*(world_size-1) ; commIndex += 2) {
			// Get sender (commList[commIndex]) & receiver (commList[commIndex+1])
			bool isSender = false;
			bool isReceiver = false;
			if (rank == commList[commIndex]) isSender = true;
			if (rank == commList[commIndex+1]) isReceiver = true;
			int klocal = 0;
			int i = 0;
			int j = 0;
			double *CSRbuf = malloc(2*sizeof(double));
			if (isSender)
			{
				if (commList[commIndex] > commList[commIndex+1]){ // If sender ID is greater than receiver ID
					// Send the upper boundary values
					klocal = 0;
					for (size_t index = 0; index < nodeX*nodeY; index ++){
						j = floor(index/nodeX);
						i = index - j*nodeX;
						CSRbuf[0] = csr_vector_get(&concentration,i+j*nodeX+klocal*nodeX*nodeY);
						MPI_Send(&CSRbuf, 1, MPI_DOUBLE, commList[commIndex+1], 0, SUB_COMM);
					}
				}
				else
				{
					// Send the lower boundary values
					klocal = thicknessMPI-1;
					for (size_t index = 0; index < nodeX*nodeY; index ++){
						j = floor(index/nodeX);
						i = index - j*nodeX;
						CSRbuf[0] = csr_vector_get(&concentration,i+j*nodeX+klocal*nodeX*nodeY);
						MPI_Send(&CSRbuf, 1, MPI_DOUBLE, commList[commIndex+1], 0, SUB_COMM);
					}
				}
			}
			else if (isReceiver){
				if (commList[commIndex] > commList[commIndex+1]){ // If sender ID is greater than receiver ID
					// Get the upper boundary values
					klocal = thicknessMPI-1;
					for (size_t index = 0; index < nodeX*nodeY; index ++){
						j = floor(index/nodeX);
						i = index - j*nodeX;
						MPI_Recv(&CSRbuf[1], 1, MPI_DOUBLE, commList[commIndex], 0, SUB_COMM, MPI_STATUS_IGNORE);
						csr_vector_set(&c_, CSRbuf[1], i+j*nodeX+klocal*nodeX*nodeY);
					}
				}
				else
				{
					// Get the lower boundary values
					klocal = 0;
					for (size_t index = 0; index < nodeX*nodeY; index ++){
						j = floor(index/nodeX);
						i = index - j*nodeX;
						MPI_Recv(&CSRbuf[1], 1, MPI_DOUBLE, commList[commIndex], 0, SUB_COMM, MPI_STATUS_IGNORE);
						csr_vector_set(&c_, CSRbuf[1], i+j*nodeX+klocal*nodeX*nodeY);
					}
				}
			}
		}


		for(int index=0; index<stopIndex ; index++)
		{
			int k = floor(index/(nodeX*nodeY));
			int j = floor((index-k*nodeX*nodeY)/nodeX);
			int i = index - k * nodeX * nodeY - j * nodeX;

			int kbis = k+1;
			int jbis = floor((index+nodeX*nodeY-kbis*nodeX*nodeY)/nodeX);
			int ibis = index+nodeX*nodeY - kbis * nodeX * nodeY - jbis * nodeX;
			csr_vector_set(&c_, csr_vector_get(&concentration,i+j*nodeX+k*nodeX*nodeY), ibis+jbis*nodeX+kbis*nodeX*nodeY);
		}

		// Check if files should be savec
		if (iteration%parameters.S == 0){
			// ============================== Writing the output file
			// Use of the MPI file IO functions
			int data_size = thicknessMPI*nodeX*nodeY; // doubles, data every node has (this is a number, not bytes!)
			MPI_File output_file;
			char file_name[20];
			sprintf(file_name, "results/c_%ld.dat",iteration);
			unsigned int N[] = {nodeX};

			MPI_File_open(SUB_COMM, file_name, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &output_file);
			if (rank == 0) MPI_File_write(output_file, N, 1, MPI_UNSIGNED, MPI_STATUS_IGNORE);

			MPI_Offset disp;
			disp = rank*(thicknessMPI-nbAdditionalSlices)*nodeX*nodeY*sizeof(double) + sizeof(unsigned int); // Displacement in bytes
			// Set the view the current node has
			MPI_File_seek(output_file, disp, MPI_SEEK_SET);

			int chunk_size = nodeX*nodeY*thicknessMPI;
			double *buffer = malloc(chunk_size*sizeof(double));
			for(int i=0; i < thicknessMPI ; ++i){ // For all slices
				for(int j=0; j < nodeY ; ++j){ // For all y values in the current slice
					for(int k=0; k < nodeX ; ++k){ // For all x values in the current line
						buffer[k + nodeX*j + nodeX*nodeY*i] = csr_vector_get(&concentration,k + nodeX*j + nodeY*nodeX*i); // Put all slices of data in the buffer
					}
				}
			}
			MPI_File_write(output_file, buffer, chunk_size, MPI_DOUBLE, MPI_STATUS_IGNORE);
			MPI_File_close(&output_file);
			free(buffer);
		}
		++iteration;
	}

	// Cleaning
	free(stopFlagsFromOthers);
	free(stopFlags);
	free(share);
	free(ranks);
	MPI_Group_free(&world_group);
	MPI_Group_free(&sub_world_group);
	if (!isIdle) MPI_Comm_free(&SUB_COMM);
	if (rank == 0){
		end = clock();
		time_spent = (double)(end - begin) / (CLOCKS_PER_SEC);
		printf("\nJob done in %2.4lf s, using %d nodes.\n", time_spent, world_size);
	}
	csr_vector_free(&concentration);
	csr_vector_free(&c_);
	MPI_Finalize();
	return 0;
}