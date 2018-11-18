#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <mpi.h>


#include "CSR_BSR.h"
#include "readData.h"

#define initConcentration 1 //[g/m3]

double round(double value, bool* isSup)
{
	if (floor(value)-value > value- ceil(value))
	{
		isSup = true;
		return floor(value)
	}
	else
	{
		isSup = false;
		return ceil(value)
	}
}


int main(int argc, char *argv[])
{
	if (argc != 2) {
		printf("Wrong arguments\n");
		printf("Please use the function with ./exe param.dat\n");
		return 1;
	}

	MPI_Init(NULL, NULL);
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	if (rank == 3) printf("We have %d processes\n", world_size);
	//if (rank == 0) printf("Reading file %s ...\n", argv[1]);
	Param parameters = readDat(argv[1]);

	size_t nodeX = (int)(parameters.L/parameters.h) + 1;
	size_t nodeY = nodeX, nodeZ =nodeX;
	if (rank == 3) printf("Number of nodes: %zu\n", nodeX);

	size_t centerIndex = nodeX*nodeY*floor(nodeZ/2)+ floor(nodeY/2)*nodeX + floor(nodeX/2);
	//printf("Index of center: %zu\n", centerIndex);

	bool isSup= true;
	size_t thicknessMPI = round(nodeZ/world_size, &isSup);
	if (rank == world_size-1)
	{
		if (isSup)
			thicknessMPI--;
		else
			thicknessMPI++;
	}
	printf("Thickness = %zu, for rank %d", thicknessMPI, rank);

	size_t kCenter = floor(centerIndex/(nodeX*nodeY));
	size_t klocalCenter = kCenter - floor(world_size/2)*thicknessMPI;

	double *concentration = calloc(nodeX*nodeY*thicknessMPI, sizeof(double));
	double *c_ = calloc(nodeX*nodeY*(thicknessMPI+2), sizeof(double)); // +2 to store the values coming from the previous and next slices
	//printf("size of c_prev %zu from process %d\n", nodeX*nodeY*(thicknessMPI+2), rank);

	if (concentration == NULL || c_ == NULL) {
		puts("Mem ERR0R !");
		exit(1);
	}

	// Give initial concentration value
	if (rank == floor(world_size/2))
		c_[nodeX*nodeY + nodeX/2+ nodeZ/2 * nodeX + klocalCenter*nodeX*nodeZ] = initConcentration; // !!! check with k

	size_t stopTime = parameters.Tmax/parameters.m;
	//printf("Stop time: %zu\n", stopTime);
	size_t iteration = 0;
	bool onBoundary = false;
	bool onZBoundary = false;
	bool valueOnBoundary= false;

	while (iteration <= stopTime && !valueOnBoundary) // ! valueOnBoundary, un seul process s'arrete !
	{
		//printf("iteration %zu from process %d\n", iteration, rank);
		//research for boundaries
		size_t isXbound = 0;
		size_t index = 0;
		if (rank == 0)
			index= nodeX*nodeY;
		size_t stopIndex = nodeX*nodeY*thicknessMPI;
		if (rank == world_size-1)
			stopIndex -=nodeY*nodeX;
		//if (rank == 3) printf("stop index %zu from process %d\n", stopIndex, rank);

		// Compute internal values
		for(index=113; index<stopIndex; index++)
		{
			//if (rank == 3) printf("enter the for index %zu from process %d\n", index, rank);
			int stage = floor(index/(nodeX*nodeY)); // !!! check with k
			//printf("Stage %d from process %zu\n", stage, rank);
			int inStage0 = index-stage*nodeX*nodeY; // !!! check with k

			if (!(isXbound==0 || isXbound==nodeX-1 || inStage0<nodeX || inStage0>=nodeX*nodeY-nodeX-1))
			{//I am in the domain
				int k = floor(index/(nodeX*nodeY)); // !!! check with k
				int j = floor((index-k*nodeX*nodeY)/nodeX);
				int i = index - k * nodeX * nodeY - j * nodeX;
				/*if (rank == 3)
				{
					printf("i j k %d %d %d\n", i, j, k);
					printf("vector %d\n", i+j*nodeX+k*nodeX*nodeY);
					printf("vector size %d\n", nodeX*nodeY*thicknessMPI);
				}*/
				int kbis = k+1;
				int jbis = floor((index+nodeX*nodeY-kbis*nodeX*nodeY)/nodeX);
				int ibis = index+nodeX*nodeY - kbis * nodeX * nodeY - jbis * nodeX;

				/*if (rank == 3)
				{
					printf("Cprev i j k %d %d %d\n", ibis, jbis, kbis);
					printf("Cprev vector %d\n", ibis+jbis*nodeX+kbis*nodeX*nodeY);
					printf("Cprev vector size %d\n", nodeX*nodeY*(thicknessMPI+2));
				}*/

				concentration[i+j*nodeX+k*nodeX*nodeY] = c_[ibis+jbis*nodeX+kbis*nodeX*nodeY] + // !!! check with k
					parameters.m * parameters.D * (c_[ibis+1+jbis*nodeX+kbis*nodeX*nodeY]+c_[ibis+(jbis+1)*nodeX+kbis*nodeX*nodeY]+
					c_[ibis+jbis*nodeX+(kbis+1)*nodeX*nodeY]-6*c_[ibis+jbis*nodeX+kbis*nodeX*nodeY]+
					c_[ibis-1+jbis*nodeX+kbis*nodeX*nodeY]+c_[ibis+(jbis-1)*nodeX+kbis*nodeX*nodeY]+
					c_[ibis+jbis*nodeX+(kbis-1)*nodeX*nodeY])/pow(parameters.h,2) -
					parameters.m * parameters.vx * (c_[ibis+1+jbis*nodeX+kbis*nodeX*nodeY]-c_[ibis-1+jbis*nodeX+kbis*nodeX*nodeY])/(2*parameters.h) -
					parameters.m * parameters.vy * (c_[ibis+(jbis+1)*nodeX+kbis*nodeX*nodeY]-c_[ibis+(jbis-1)*nodeX+kbis*nodeX*nodeY])/(2*parameters.h) -
					parameters.m * parameters.vz * (c_[ibis+jbis*nodeX+(kbis+1)*nodeX*nodeY]-c_[ibis+jbis*nodeX+(kbis-1)*nodeX*nodeY])/(2*parameters.h);


				if (rank==0 || rank ==world_size-1) onZBoundary = (index<=2*nodeX*nodeY || index >(thicknessMPI-2)*nodeX*nodeX);
				//{printf("onZ comparison from process %d\n", rank); onZBoundary = (index<=2*nodeX*nodeY || index >(thicknessMPI-2)*nodeX*nodeX); printf("onZ comparison works from process %d\n", rank);}
				onBoundary = ((isXbound == nodeX-2) || (isXbound == 1) || (inStage0 >= nodeX && inStage0 <= 2*nodeX-2) || (inStage0 >= nodeY*nodeY-2*nodeX-1));
				//printf("onZ %d, onbound %d from process %d\n", onZBoundary, onBoundary, rank);
				if ((onBoundary||onZBoundary)  && concentration[i+j*nodeX+k*nodeX*nodeY] != 0) valueOnBoundary=true;
				//{printf("onbound comparison from process %d\n", rank); valueOnBoundary=true; printf("on comparison works from process %d\n", rank);}
			}

			isXbound++;
			if(isXbound==nodeX) isXbound = 0;
			//if (rank == 3) printf("index reached %zu from process %d\n", index, rank);
		}

		//printf("\n\n for loop works :D from process %d at iteration %zu\n\n", rank, iteration);
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
				int klocal = 0;
				int j = floor((index-klocal*nodeX*nodeY)/nodeX);
				int i = index - klocal * nodeX * nodeY - j * nodeX;
				if (isSender)
				{
					if (commList[commIndex] > commList[commIndex+1])
						{ // If sender ID is greater than receiver ID
							// Send the upper boundary values
							klocal = 0;
							for (size_t index = 0; index < nodeX*nodeY; index ++)
								MPI_Send(&concentration[i+j*nodeX+klocal*nodeX*nodeY], 1, MPI_DOUBLE, commList[commIndex+1], 0, MPI_COMM_WORLD);
						}
						else
						{
							// Send the lower boundary values
							klocal = thicknessMPI-1;
							for (size_t index = 0; index < nodeX*nodeY; index ++)
								MPI_Send(&concentration[i+j*nodeX+klocal*nodeX*nodeY], 1, MPI_DOUBLE, commList[commIndex+1], 0, MPI_COMM_WORLD);
						}
					}
					else
					{
						if (commList[commIndex] > commList[commIndex+1]){ // If sender ID is greater than receiver ID
							// Get the upper boundary values
							klocal = thicknessMPI-1;
							for (size_t index = 0; index < nodeX*nodeY; index ++)
								MPI_Recv(&c_[i+j*nodeX+klocal*nodeX*nodeY], 1, MPI_DOUBLE, commList[commIndex], 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
						}
						else
						{
							// Get the lower boundary values
							klocal = 0;
							for (size_t index = 0; index < nodeX*nodeY; index ++)
								MPI_Recv(&c_[i+j*nodeX+klocal*nodeX*nodeY], 1, MPI_DOUBLE, commList[commIndex], 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
						}
					}
				}
			}
		/*if (!(iteration%500))
		{*/
			for(int index=0; index<stopIndex ; index++)
			{
				int k = floor(index/(nodeX*nodeY)); // !!! check with k
				int j = floor((index-k*nodeX*nodeY)/nodeX);
				int i = index - k * nodeX * nodeY - j * nodeX;

				int kbis = k+1;
				int jbis = floor((index+nodeX*nodeY-kbis*nodeX*nodeY)/nodeX);
				int ibis = index+nodeX*nodeY - kbis * nodeX * nodeY - jbis * nodeX;
				c_[ibis+jbis*nodeX+kbis*nodeX*nodeY] = concentration[i+j*nodeX+k*nodeX*nodeY];
				//printf("%f ", concentration[i]);
			}

			/*printf("\n \n \n");
		} */
			//printf("iteration %zu ended\n", iteration);
			++iteration;
		}

		MPI_Finalize();
		return 0;
	}
