#include <stdio.h>
#include <math.h>
#include <stdlib.h>
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
	if (rank == 0)MPI_File_delete("out.dat", MPI_INFO_NULL);

	// Retrieving data from the .dat file
	Param parameters = readDat(argv[1]);
	size_t nodeX = (int)(parameters.L/parameters.h) + 1;
	size_t nodeY = nodeX, nodeZ =nodeX;
	size_t centerIndex = nodeX*nodeY*floor(nodeZ/2)+ floor(nodeY/2)*nodeX + floor(nodeX/2);
	size_t stopTime = parameters.Tmax/parameters.m;

	// Assign each node its nomber of slices (thicknessMPI)
	size_t thicknessMPI = myRound(value);
	if (rank == world_size-1)
	{
		if (thicknessMPI>nodeZ/world_size)
			thicknessMPI--;
		else
			thicknessMPI++;
	}

	// Some prints
	if (rank == 0) printf("We have %d processes\n", world_size);
	if (rank == 0) printf("Number of slices: %zu\n", nodeZ);
	printf("Thickness = %zu, for rank %d\n", thicknessMPI, rank);


	// Get memory for the concentration values
	double *concentration = calloc(nodeX*nodeY*thicknessMPI, sizeof(double));
	double *c_ = calloc(nodeX*nodeY*(thicknessMPI+2), sizeof(double)); // +2 to store the values coming from the previous and next slices
	if (concentration == NULL || c_ == NULL) {
		puts("Mem ERR0R !");
		exit(1);
	}

	// Give initial concentration value
	size_t kCenter = floor(centerIndex/(nodeX*nodeY));
	size_t klocalCenter = kCenter - floor(world_size/2)*thicknessMPI;
	if (rank == floor(world_size/2))
		c_[nodeX*nodeY + nodeX/2+ nodeZ/2 * nodeX + klocalCenter*nodeX*nodeZ] = initConcentration; // !!! check with k

	// Declare variables for the main loop
	size_t iteration = 0;
	bool onBoundary = false;
	bool onZBoundary = false;
	bool valueOnBoundary = false;
	int *stopFlags = calloc(world_size, sizeof(int));
	int *stopFlagsFromOthers = calloc(world_size, sizeof(int));
	bool stopFlag = false;

	/*================================================================================================
	#	Main loop
	================================================================================================*/
	while (iteration <= stopTime && !stopFlag) // ! valueOnBoundary, un seul process s'arrete !
	{
		//if (rank == 1 && iteration>0) printf("enter the iteration loop from process %d\n", rank);
		//printf("iteration %zu from process %d\n", iteration, rank);
		//research for boundaries
		size_t isXbound = 0;
		size_t index = 0;
		if (rank == 0) index = nodeX*nodeY;
		size_t stopIndex = nodeX*nodeY*thicknessMPI;
		if (rank == world_size-1) stopIndex -= nodeY*nodeX;

		// ============================== Compute internal values
		for(index=0; index<stopIndex; index++)
		{
			//if (rank == 1 && iteration>0) printf("enter the for index %zu from process %d\n", index, rank);
			int stage = floor(index/(nodeX*nodeY)); // !!! check with k
			//printf("Stage %d from process %zu\n", stage, rank);
			int inStage0 = index-stage*nodeX*nodeY; // !!! check with k

			if (!(isXbound==0 || isXbound==nodeX-1 || inStage0<nodeX || inStage0>=nodeX*nodeY-nodeX)) // !!! -1 enlevé ici
			{//I am in the domain
				int k = floor(index/(nodeX*nodeY)); // !!! check with k
				int j = floor((index-k*nodeX*nodeY)/nodeX);
				int i = index - k * nodeX * nodeY - j * nodeX;
				/*if (rank == 1 && iteration >0)
				{
					printf("i j k %ld %ld %ld\n", i, j, k);
					printf("vector %ld\n", i+j*nodeX+k*nodeX*nodeY);
					printf("vector size %ld\n", nodeX*nodeY*thicknessMPI);
				}*/

				int kbis = k+1;
				int jbis = floor((index+nodeX*nodeY-kbis*nodeX*nodeY)/nodeX);
				int ibis = index+nodeX*nodeY - kbis * nodeX * nodeY - jbis * nodeX;

				/*if (rank == 1 && iteration >0)
				{
					printf("Cprev i j k %ld %ld %ld\n", ibis, jbis, kbis);
					printf("Cprev vector %ld\n", ibis+jbis*nodeX+kbis*nodeX*nodeY);
					printf("Cprev vector size %ld\n", nodeX*nodeY*(thicknessMPI+2));
				}*/

				//printf("index, i , j, k : %d %d %d %d\n", index, i, j, k);

				concentration[i+j*nodeX+k*nodeX*nodeY] = 
					c_[ibis+jbis*nodeX+kbis*nodeX*nodeY] +
					parameters.m * parameters.D * (c_[ibis+1+jbis*nodeX+kbis*nodeX*nodeY]+c_[ibis+(jbis+1)*nodeX+kbis*nodeX*nodeY]+
					c_[ibis+jbis*nodeX+(kbis+1)*nodeX*nodeY]-6*c_[ibis+jbis*nodeX+kbis*nodeX*nodeY]+
					c_[ibis-1+jbis*nodeX+kbis*nodeX*nodeY]+c_[ibis+(jbis-1)*nodeX+kbis*nodeX*nodeY]+
					c_[ibis+jbis*nodeX+(kbis-1)*nodeX*nodeY])/pow(parameters.h,2) -
					parameters.m * parameters.vx * (c_[ibis+1+jbis*nodeX+kbis*nodeX*nodeY]-c_[ibis-1+jbis*nodeX+kbis*nodeX*nodeY])/(2*parameters.h) -
					parameters.m * parameters.vy * (c_[ibis+(jbis+1)*nodeX+kbis*nodeX*nodeY]-c_[ibis+(jbis-1)*nodeX+kbis*nodeX*nodeY])/(2*parameters.h) -
					parameters.m * parameters.vz * (c_[ibis+jbis*nodeX+(kbis+1)*nodeX*nodeY]-c_[ibis+jbis*nodeX+(kbis-1)*nodeX*nodeY])/(2*parameters.h);

				// Check the stopping conditions on the boundaries of the domain
				if (rank == 0) onZBoundary = (index<2*nodeX*nodeY);
				if (rank == world_size-1) onZBoundary = (index>=(thicknessMPI-2)*nodeX*nodeX);
				//{printf("onZ comparison from process %d\n", rank); onZBoundary = (index<=2*nodeX*nodeY || index >(thicknessMPI-2)*nodeX*nodeX); printf("onZ comparison works from process %d\n", rank);}
				onBoundary = ((isXbound == nodeX-2) || (isXbound == 1) || (inStage0 >= nodeX && inStage0 <= 2*nodeX-1) || (inStage0 >= nodeY*nodeY-2*nodeX));
				//printf("onZ %d, onbound %d from process %d\n", onZBoundary, onBoundary, rank);
				if ((onBoundary || onZBoundary)  && (concentration[i+j*nodeX+k*nodeX*nodeY] != 0)){
					valueOnBoundary=true;
					//printf("onbound comparison works at index %zu for process %d at iteration %zu and C is %f\n", index, rank, iteration, concentration[index]);
				}

				// Send your status to all the other nodes
				for (int i = 0; i < world_size; ++i){
					if(valueOnBoundary)stopFlags[i] = 1;
				}
				MPI_Allgather(stopFlags, 1, MPI_INT, stopFlagsFromOthers, 1, MPI_INT, MPI_COMM_WORLD);
				printf("On node #%d, at iteration %d and index %d stop flags from others: \n", rank, iteration, index);
				for (int i = 0; i < world_size; ++i){
					printf("%d ", stopFlagsFromOthers[i]);
					if(stopFlagsFromOthers[i] == 1){
						stopFlag = true;
						break;
					}
				}
			}

			isXbound++;
			if(isXbound==nodeX) isXbound = 0;
		}

		// ============================== Send and receive neighboring values
		int *commList = getCommListSlices(world_size);
		for (int commIndex=0 ; commIndex<4*(world_size-1) ; commIndex += 2) {
			// Get sender (commList[commIndex]) & receiver (commList[commIndex+1])
			//printf("Node #%d has to send to node #%d.\n", commList[commIndex], commList[commIndex+1]);
			bool isSender = false;
			bool isReceiver = false;
			if (rank == commList[commIndex]) isSender = true;
			if (rank == commList[commIndex+1]) isReceiver = true;
			//if (isSender) printf("Line %d in the commList: %d is sending to %d\n", commIndex, commList[commIndex], commList[commIndex+1]);
			int klocal = 0;
			int i = 0;
			int j = 0;
			if (isSender)
			{
				//printf("Node #%d tries to send to %d.\n", rank, commList[commIndex+1]);
				if (commList[commIndex] > commList[commIndex+1])
					{ // If sender ID is greater than receiver ID
						// Send the upper boundary values
						klocal = 0;
						for (size_t index = 0; index < nodeX*nodeY; index ++){
							j = floor(index/nodeX);
							i = index - j*nodeX;
							//printf("index, i , j : %d %d %d\n", index, i, j);
							MPI_Send(&concentration[i+j*nodeX+klocal*nodeX*nodeY], 1, MPI_DOUBLE, commList[commIndex+1], 0, MPI_COMM_WORLD);
						}
						//printf("Node #%d has successfully sent its data to node #%d.\n", rank, commList[commIndex+1]);
					}
					else
					{
						// Send the lower boundary values
						klocal = thicknessMPI-1;
						for (size_t index = 0; index < nodeX*nodeY; index ++){
							j = floor(index/nodeX);
							i = index - j*nodeX;
							//printf("index, i , j : %d %d %d\n", index, i, j);
							MPI_Send(&concentration[i+j*nodeX+klocal*nodeX*nodeY], 1, MPI_DOUBLE, commList[commIndex+1], 0, MPI_COMM_WORLD);
						}
						//printf("Node #%d has successfully sent its data to node #%d.\n", rank, commList[commIndex+1]);
					}
			}
			else if (isReceiver){
				//printf("Node #%d tries to receive from node #%d.\n", rank, commList[commIndex]);
				if (commList[commIndex] > commList[commIndex+1]){ // If sender ID is greater than receiver ID
					// Get the upper boundary values
					klocal = thicknessMPI-1;
					for (size_t index = 0; index < nodeX*nodeY; index ++){
							j = floor(index/nodeX);
							i = index - j*nodeX;
						MPI_Recv(&c_[i+j*nodeX+klocal*nodeX*nodeY], 1, MPI_DOUBLE, commList[commIndex], 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					}
				}
				else
				{
					// Get the lower boundary values
					klocal = 0;
					for (size_t index = 0; index < nodeX*nodeY; index ++){
							j = floor(index/nodeX);
							i = index - j*nodeX;
						MPI_Recv(&c_[i+j*nodeX+klocal*nodeX*nodeY], 1, MPI_DOUBLE, commList[commIndex], 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
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
			c_[ibis+jbis*nodeX+kbis*nodeX*nodeY] = concentration[i+j*nodeX+k*nodeX*nodeY];
			//printf("%f ", concentration[i]);
		}

		//printf("iteration %zu ended\n", iteration);
		++iteration;
	}


	/*================================================================================================
	#	Writing the output file
	================================================================================================*/
	// Use of the MPI file IO functions
	int data_size = thicknessMPI*nodeX*nodeY; // doubles, data every node has (this is a number, not bytes!)
	// Prepare the data types and dispalcements for MPI
	MPI_File output_file;
	MPI_Datatype data_type;
	MPI_Type_contiguous(data_size, MPI_DOUBLE, &data_type);
	MPI_Type_commit(&data_type);
	MPI_Offset disp;
	disp = data_size*rank*sizeof(double) + sizeof(char); // Displacement in bytes
	// Actually open/create the file and set the view the current node has
	MPI_File_open(MPI_COMM_WORLD, "out.dat", MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &output_file);
	MPI_File_seek(output_file, disp, MPI_SEEK_SET);

	int isXbound = 0;
	int chunk_size = nodeX;
	int *buffer = malloc(chunk_size*sizeof(double));
	char *carrier_return[] = {"\n"};
	for(int i=0; i<thicknessMPI ; ++i){ // For all slices
		for(int j=0; j<nodeY ; ++j){ // For all y values in the current slice
			for(int k=0; k<nodeX ; ++k){ // For all x values in the current line
				buffer[i] = concentration[i]; // Put a line of data in the buffer
			}
			// Once we have a line, write it in the file
			MPI_File_write(output_file, buffer, chunk_size, MPI_DOUBLE, MPI_STATUS_IGNORE);
			MPI_File_write(output_file, carrier_return, 1, MPI_CHAR, MPI_STATUS_IGNORE);
		}
	}

	MPI_File_close(&output_file);
	printf("Node #%d has finished writing %ld bytes in the output file.\n", rank, data_size*sizeof(double)+thicknessMPI*nodeY*sizeof(char));

	// Cleaning
	free(stopFlagsFromOthers);
	free(stopFlags);
	free(buffer);
	free(concentration);
	free(c_);
	MPI_Finalize();
	return 0;
}
