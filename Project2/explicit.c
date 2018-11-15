#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <mpi.h>


#include "CSR_BSR.h"
#include "readData.h"

#define initConcentration 1 //[g/m3]



int main(int argc, char *argv[])
{
	if (argc != 2) {
		printf("Wrong arguments\n");
		printf("Please use the function with ./exe param.dat\n");
		return 1;
	}

	printf("Reading file %s ...\n", argv[1]);
	Param parameters = readDat(argv[1]);

	size_t nodeX = (int)(parameters.L/parameters.h) + 1;
	size_t nodeY = nodeX, nodeZ =nodeX;
	printf("\nNumber of nodes: %zu\n", nodeX);

	size_t centerIndex = nodeX*nodeY*floor(nodeZ/2)+ floor(nodeY/2)*nodeX + floor(nodeX/2);
	printf("Index of center: %zu\n", centerIndex);


	MPI_Init(NULL, NULL);
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);


	/*
	// Example for the getCommListSlices() function
	int *commList0 = getCommListSlices(20);
	for (int i = 0; i < 4*19; i+=2)
	{
		printf("%d %d\n", commList0[i], commList0[i+1]);
	}*/

	char bonusSlice=0;

	//if (nodeZ%rank) bonusSlice = 1;
	size_t thicknessMPI = (int)(nodeZ/world_size);
	if (rank == world_size-1) thicknessMPI++;

	size_t kCenter = floor(centerIndex/(nodeX*nodeY));
	size_t klocalCenter = kCenter - floor(world_size/2)*thicknessMPI;
	double *concentration = calloc(nodeX*nodeY*thicknessMPI, sizeof(double));
	double *c_ = calloc(nodeX*nodeY*(thicknessMPI+2), sizeof(double)); // +2 to store the values coming from the previous and next slices

	if (concentration == NULL || c_ == NULL) {
		puts("Mem ERR0R !");
		exit(1);
	}

	if (rank == floor(world_size/2))
		c_[nodeX/2+ nodeZ/2 * nodeX + klocalCenter*nodeX*nodeZ] = initConcentration;

	size_t stopTime = parameters.Tmax/parameters.m;
	printf("Stop time: %zu\n", stopTime);
	size_t iteration = 0;
	bool onBoundary = false;
	bool onZBoundary = false;
	bool valueOnBoundary= false;

	while (iteration <= stopTime && !valueOnBoundary)
	{
		//printf("iteration %zu\n", iteration);
		//research for boundaries
		size_t isXbound = 0;
		size_t index = 0;
		if (rank == 0)
			index=nodeX*nodeY;
		size_t stopIndex = nodeX*nodeY*thicknessMPI;
		if (rank == world_size-1)
			stopIndex -=nodeY*nodeX;

		// Compute internal values
		for(index; index<stopIndex; index++)
		{
			int stage = floor(index/(nodeX*nodeY));
			int inStage0 = index-stage*nodeX*nodeY;

			if (!(isXbound==0 || isXbound==nodeX-1 || inStage0<nodeX || inStage0>=nodeX*nodeY-nodeX-1))
			{//I am in the domain
				int k = floor(index/(nodeX*nodeY));
				int j = floor((index-k*nodeX*nodeY)/nodeX);
				int i = index - k * nodeX * nodeY - j * nodeX;
				concentration[i+j*nodeX+k*nodeX*nodeY] = c_[i+j*nodeX+k*nodeX*nodeY] +
				parameters.m * parameters.D * (c_[i+1+j*nodeX+k*nodeX*nodeY]+c_[i+(j+1)*nodeX+k*nodeX*nodeY]+
					c_[i+j*nodeX+(k+1)*nodeX*nodeY]-6*c_[i+j*nodeX+k*nodeX*nodeY]+
					c_[i-1+j*nodeX+k*nodeX*nodeY]+c_[i+(j-1)*nodeX+k*nodeX*nodeY]+
					c_[i+j*nodeX+(k-1)*nodeX*nodeY])/pow(parameters.h,2) -
				parameters.m * parameters.vx * (c_[i+1+j*nodeX+k*nodeX*nodeY]-c_[i-1+j*nodeX+k*nodeX*nodeY])/(2*parameters.h) -
				parameters.m * parameters.vy * (c_[i+(j+1)*nodeX+k*nodeX*nodeY]-c_[i+(j-1)*nodeX+k*nodeX*nodeY])/(2*parameters.h) -
				parameters.m * parameters.vz * (c_[i+j*nodeX+(k+1)*nodeX*nodeY]-c_[i+j*nodeX+(k-1)*nodeX*nodeY])/(2*parameters.h);

				if (rank==0 || rank ==world_size-1) onZBoundary = (index<=2*nodeX*nodeY || index >(thicknessMPI-2)*nodeX*nodeX);
				onBoundary = ((isXbound == nodeX-2) || (isXbound == 1) || (inStage0 >= nodeX && inStage0 <= 2*nodeX-2) || (inStage0 >= nodeY*nodeY-2*nodeX-1));
				if ((onBoundary||onZBoundary)  && concentration[i+j*nodeX+k*nodeX*nodeY] != 0) valueOnBoundary=true;
			}

			isXbound++;
			if(isXbound==nodeX) isXbound = 0;
		}

		// Send and receive neighboring values
		int *commList = getCommListSlices(world_size);
		for (int commIndex=0 ; commIndex<4*(world_size-1) ; commIndex += 2) {
			// Get sender (commList[commIndex]) & receiver (commList[commIndex+1])
			//printf("%d to %d\n", commList[commIndex], commList[commIndex+1]);
			bool isSender = false;
			bool isReceiver = false;
			if (rank == commList[commIndex]) isSender = true;
			if (rank == commList[commIndex+1]) isReceiver = true;
			//if (isSender) printf("Line %d in the commList: %d is sending to %d\n", commIndex, commList[commIndex], commList[commIndex+1]);
			// Need a if here, depending on the rank so that everyone knows what to do asynchronously
			if (isSender || isReceiver)
			{
				for (size_t index = 0; index < nodeX*nodeY; index ++)
				{
					int klocal = 0;
					int j = floor((index-klocal*nodeX*nodeY)/nodeX);
					int i = index - klocal * nodeX * nodeY - j * nodeX;
					if (isSender)
					{
						// Send the lower boundary values
						if (rank!=0)
							MPI_Send(&concentration[i+j*nodeX+klocal*nodeX*nodeY], 1, MPI_DOUBLE, commList[commIndex+1], 0, MPI_COMM_WORLD);
						// Send the uppwer boundary values
						klocal = thicknessMPI-1;//max
						if (rank!=world_size-1)
							MPI_Send(&concentration[i+j*nodeX+klocal*nodeX*nodeY], 1, MPI_DOUBLE, commList[commIndex+1], 0, MPI_COMM_WORLD);
					}
					else
					{
						// Get the upper boundary values
						klocal = thicknessMPI-1;//max
						if (rank!=world_size-1)
							MPI_Recv(&c_[i+j*nodeX+klocal*nodeX*nodeY], 1, MPI_DOUBLE, commList[commIndex], 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
						// Get the lower boundary values
						klocal = 0;//max
						if (rank!=0)
							MPI_Recv(&c_[i+j*nodeX+klocal*nodeX*nodeY], 1, MPI_DOUBLE, commList[commIndex], 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					}
				}
			}
		}
		/*if (!(iteration%500))
		{
			for(int i=0; i<nodeX*nodeY*nodeZ ; i++)
			{
				c_[i] = concentration[i];
				printf("%f ", concentration[i]);
			}

			printf("\n \n \n");
		} */
		++iteration;
	}

	MPI_Finalize();
	return 0;
}
