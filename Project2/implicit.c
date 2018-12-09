/*=======================================================================================
* This code was written by:                                                          *
*               Antonin Aumètre - antonin.aumetre@gmail.com            *
*               Céline Moureau -  cemoureau@gmail.com                *
* For: High Performance Scientific course at ULiège, 2018-19                         *
* Project 2                                                                          *
*                                                            *
* Originally uploaded to: https://github.com/Cobalt1911                              *
* Under GNU General Public License 11/2018                                           *
=======================================================================================*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>

#include "COO_CSR_BSR.h"
#include "algorithms.h"
#include "fileIO.h"

#define initConcentration 1 //[g/m3]

/* ----------------------------------------------------------------------------*
* function: SquaredNorm
* -----------------------------------------------------------------------------*
* @ param:
* @ return:
* -----------------------------------------------------------------------------*
*
* ----------------------------------------------------------------------------*/
double SquaredNorm(double* vect, size_t sizeVect)
{
  double normVal=0;
  for(size_t i = 0; i< sizeVect; i++)
    normVal+=vect[i]*vect[i];
  return normVal;
}

/* ----------------------------------------------------------------------------*
* function: VectorProduct
* -----------------------------------------------------------------------------*
* @ param:
* @ return:
* -----------------------------------------------------------------------------*
*
* ----------------------------------------------------------------------------*/
double VectorProduct(double* vect1, double* vect2, size_t sizeVect)
{
  double prod=0;
  for(size_t i = 0; i< sizeVect; i++)
    prod+=vect1[i]*vect2[i];
  return prod;
}

/* ----------------------------------------------------------------------------*
* function: SumVect
* -----------------------------------------------------------------------------*
* @ param:
* @ return:
* -----------------------------------------------------------------------------*
*
* ----------------------------------------------------------------------------*/
void SumVect(double* vectToStore, double* vect1, double* vect2, double multVal, size_t sizeVect)
{
  for(size_t i = 0; i< sizeVect; i++)
    vectToStore[i] = vect1[i] + vect2[i]*multVal;
}

/* ----------------------------------------------------------------------------*
* function: Ap
* -----------------------------------------------------------------------------*
* @ param:
* @ return:
* -----------------------------------------------------------------------------*
*
* ----------------------------------------------------------------------------*/
void Ap(double* p, double* Apresult, size_t nodeX, size_t nodeY, size_t thicknessMPI, double h, double m, double vx, double vy, double vz, double D, int rank, int world_size)
{
  for(size_t index = 0; index< nodeX*nodeY*thicknessMPI; index++)
  {
    bool onZBoundary = false;
    int k = floor(index/(nodeX*nodeY));
    int j = floor((index-k*nodeX*nodeY)/nodeX);
    int i = index - k * nodeX * nodeY - j * nodeX;
    double Ddivh2 = D/(h*h);
    if (rank == 0) onZBoundary = (index<nodeX*nodeY);
    if (rank == world_size-1) onZBoundary = (index>=(thicknessMPI-1)*nodeX*nodeX);
    if(!(i==0 || i == nodeX-1 || j==0 || j == nodeY-1 || onZBoundary))
    {
      k+=1;
      Apresult[index] = p[i+j*nodeX+k*nodeX*nodeY]*(1/m + 6*Ddivh2) + p[i-1+j*nodeX+k*nodeX*nodeY]*(-vx/(2*h) - Ddivh2) + p[i+1+j*nodeX+k*nodeX*nodeY]*(vx/(2*h) - Ddivh2) + p[i+(j-1)*nodeX+k*nodeX*nodeY]*(-vy/(2*h) - Ddivh2) + p[i+(j+1)*nodeX+k*nodeX*nodeY]*(vy/(2*h) - Ddivh2) + p[i+j*nodeX+(k-1)*nodeX*nodeY]*(-vz/(2*h) -Ddivh2) + p[i+j*nodeX+(k+1)*nodeX*nodeY]*(vz/(2*h) - Ddivh2);
      Apresult[index] *= m;
    }
  }
}

/* ----------------------------------------------------------------------------*
* MAIN
* -----------------------------------------------------------------------------*
* @ usage
* -----------------------------------------------------------------------------*
*
* ----------------------------------------------------------------------------*/
int implicit_solver(int argc, char *argv[])
{
  // Initalizes MPI
	MPI_Init(NULL, NULL);
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  //pour après
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
	if (rank == 0) printf("threshold: %f\n", parameters.rthreshold);
	MPI_Barrier(MPI_COMM_WORLD);
	if (!isIdle) printf("Thickness = %zu, for rank %d\n", thicknessMPI, rank);
	MPI_Barrier(MPI_COMM_WORLD);
	if (rank == 0) printf("=====================================\n");

  double *concentration = calloc(nodeX*nodeY*thicknessMPI, sizeof(double));

  if (concentration == NULL)
  {
    puts("Mem ERR0R !");
    exit(1);
  }

  // Give initial concentration value
  int middleSliceIndex = floor(nodeZ/2);
  int rankMiddle = floor(middleSliceIndex/(thicknessMPI-nbAdditionalSlices));
  if (rank == rankMiddle)
  {
    int initValueIndex = nodeX*floor(nodeY/2) + floor(nodeX/2) + nodeX*nodeY*(middleSliceIndex-rank*(thicknessMPI-nbAdditionalSlices));
    concentration[initValueIndex] = initConcentration;
    printf("Initial value set on node #%d, at index %d\n", rank, initValueIndex);
  }
  MPI_Barrier(MPI_COMM_WORLD);

  size_t iteration = 0;
  while (iteration <= stopTime && !stopFlag)
  {
    if (rank == 0) printf("%ld ", iteration);
    //--------------------------------------------------------------------------
    //              Conjugate Gradient Method
    //--------------------------------------------------------------------------
    // x0 = 0
    double *concentrationSuiv = calloc(nodeX*nodeY*thicknessMPI, sizeof(double));

    // r0 = b - Ax0 = b
    double *r = calloc(nodeX*nodeY*thicknessMPI, sizeof(double));
    double *rsuiv = calloc(nodeX*nodeY*thicknessMPI, sizeof(double));

    if (r == NULL || rsuiv == NULL || concentrationSuiv == NULL)
    {
      puts("Mem ERR0R concentrationSuiv or r or rsuiv!");
      exit(1);
    }

    for (size_t copyIndex = 0; copyIndex < nodeX*nodeY*thicknessMPI; copyIndex++)
    {
        r[copyIndex] = concentration[copyIndex];
    }

    // p0 = r0
    double *p = calloc(nodeX*nodeY*thicknessMPI, sizeof(double));

    if (p == NULL)
    {
      puts("Mem ERR0R p!");
      exit(1);
    }

    for (size_t copyIndex = 0; copyIndex < nodeX*nodeY*thicknessMPI; copyIndex++)
    {
        p[copyIndex] = r[copyIndex];
    }

    double local_sum = SquaredNorm(r, nodeX*nodeY*thicknessMPI);
    double global_sum, savedR0;
    MPI_Allreduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, SUB_COMM);
    savedR0=global_sum;
    while(sqrt(global_sum)/sqrt(savedR0)>=parameters.rthreshold)
    {
      double *Apvect = calloc(nodeX*nodeY*thicknessMPI, sizeof(double));

      if (Apvect == NULL)
      {
        puts("Mem ERR0R Apvect!");
        exit(1);
      }
      //Get neighbours
      double *p2 = calloc(nodeX*nodeY*(thicknessMPI+2), sizeof(double));
      if (p2 == NULL)
      {
        puts("Mem ERR0R p!");
        exit(1);
      }

      for (size_t copyIndex = 0; copyIndex< nodeX*nodeY*thicknessMPI; copyIndex++)
      {
        p2[copyIndex+nodeX*nodeY] = p[copyIndex];
      }

      // ============================== Send and receive neighboring values
  		int *commList = getCommListSlices(world_size);
  		for (int commIndex=0 ; commIndex<4*(world_size-1) ; commIndex += 2)
  		{
  			bool isSender = false;
  			bool isReceiver = false;
  			if (rank == commList[commIndex]) isSender = true;
  			if (rank == commList[commIndex+1]) isReceiver = true;
  			int klocal = 0;
  			int i = 0;
  			int j = 0;
  			if (isSender)
  			{
  				if (commList[commIndex] > commList[commIndex+1])
  					{ // If sender ID is greater than receiver ID
  						// Send the upper boundary values
  						klocal = 0;
  						for (size_t index = 0; index < nodeX*nodeY; index ++){
  							j = floor(index/nodeX);
  							i = index - j*nodeX;
  							MPI_Send(&p[i+j*nodeX+klocal*nodeX*nodeY], 1, MPI_DOUBLE, commList[commIndex+1], 0, SUB_COMM);
  						}
  					}
  					else
  					{
  						// Send the lower boundary values
  						klocal = thicknessMPI-1;
  						for (size_t index = 0; index < nodeX*nodeY; index ++){
  							j = floor(index/nodeX);
  							i = index - j*nodeX;
  							MPI_Send(&p[i+j*nodeX+klocal*nodeX*nodeY], 1, MPI_DOUBLE, commList[commIndex+1], 0, SUB_COMM);
  						}
  					}
  			}
  			else if (isReceiver){
  				if (commList[commIndex] > commList[commIndex+1]){ // If sender ID is greater than receiver ID
  					// Get the upper boundary values
  					klocal = thicknessMPI+1;
  					for (size_t index = 0; index < nodeX*nodeY; index ++){
  							j = floor(index/nodeX);
  							i = index - j*nodeX;
  						MPI_Recv(&p2[i+j*nodeX+klocal*nodeX*nodeY], 1, MPI_DOUBLE, commList[commIndex], 0, SUB_COMM, MPI_STATUS_IGNORE);
  					}
  				}
  				else
  				{
  					// Get the lower boundary values
  					klocal = 0;
  					for (size_t index = 0; index < nodeX*nodeY; index ++){
  							j = floor(index/nodeX);
  							i = index - j*nodeX;
  						MPI_Recv(&p2[i+j*nodeX+klocal*nodeX*nodeY], 1, MPI_DOUBLE, commList[commIndex], 0, SUB_COMM, MPI_STATUS_IGNORE);
  					}
  				}
  			}
  		}
      Ap(p2, Apvect, nodeX, nodeY, thicknessMPI, parameters.h, parameters.m, parameters.vx, parameters.vy, parameters.vz, parameters.D, rank, world_size);
      double rtr = SquaredNorm(r, nodeX*nodeY*thicknessMPI);
      double global_rtr;
      MPI_Allreduce(&rtr, &global_rtr, 1, MPI_DOUBLE, MPI_SUM, SUB_COMM);
      double ptAp = VectorProduct(p, Apvect, nodeX*nodeY*thicknessMPI);
      double global_ptAp;
      MPI_Allreduce(&ptAp, &global_ptAp, 1, MPI_DOUBLE, MPI_SUM, SUB_COMM);
      double alpha = global_rtr / global_ptAp;

      // xi+1 = xi + alpha*p
      SumVect(concentrationSuiv, concentrationSuiv, p, alpha, nodeX*nodeY*thicknessMPI);

      // ri+1 = ri - alpha*Ap
      SumVect(rsuiv, r, Apvect, -alpha, nodeX*nodeY*thicknessMPI);

      // beta = ri+1^T*ri+1 / ri^T*ri
      double rSuitrSui = SquaredNorm(rsuiv, nodeX*nodeY*thicknessMPI);
      double global_rSuitrSui;
      MPI_Allreduce(&rSuitrSui, &global_rSuitrSui, 1, MPI_DOUBLE, MPI_SUM, SUB_COMM);
      double beta = global_rSuitrSui / global_rtr;

      // pi+1 = ri+1 + beta*p
      SumVect(p, rsuiv, p, beta, nodeX*nodeY*thicknessMPI);
      for (size_t copyIndex = 0; copyIndex < nodeX*nodeY*thicknessMPI; copyIndex++)
      {
          r[copyIndex] = rsuiv[copyIndex];
      }
      local_sum = SquaredNorm(r, nodeX*nodeY*thicknessMPI);
      MPI_Allreduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, SUB_COMM);

      free(Apvect);
      free(p2);
    }

    for(size_t index = 0; index< nodeX*nodeY*thicknessMPI; index++)
    {
      bool onZBoundary = false;
      int k = floor(index/(nodeX*nodeY));
      int j = floor((index-k*nodeX*nodeY)/nodeX);
      int i = index - k * nodeX * nodeY - j * nodeX;
      if (rank == 0) onZBoundary = (index<2*nodeX*nodeY);
      if (rank == world_size-1) onZBoundary = (index>=(thicknessMPI-2)*nodeX*nodeX);
      if((i<=1 || i >= nodeX-2 || j<=1 || j >= nodeY-2 || onZBoundary) && concentrationSuiv[index]> 5e-8)
      {
          printf("conventration on boundary = %e at index %zu\n", concentrationSuiv[index], index);
          printf("STOP\n");
          valueOnBoundary = true;
          break;
      }
    }

    // Send your status to all the other nodes
		for (int i = 0; i < world_size; ++i)
		{
			if(valueOnBoundary)stopFlags[i] = 1;
		}
		MPI_Allgather(stopFlags, 1, MPI_INT, stopFlagsFromOthers, 1, MPI_INT, SUB_COMM);
		for (int i = 0; i < world_size; ++i)
		{
			if(stopFlagsFromOthers[i] == 1)
			{
			stopFlag = true;
			break;
			}
		}

    for (size_t copyIndex = 0; copyIndex < nodeX*nodeY*thicknessMPI; copyIndex++)
    {
        concentration[copyIndex] = concentrationSuiv[copyIndex];
    }

    // Check if files should be saved
    if (iteration%parameters.S == 0)
    {
      // ============================== Writing the output file
      // Use of the MPI file IO functions
      int data_size = thicknessMPI*nodeX*nodeY; // doubles, data every node has (this is a number, not bytes!)
      MPI_File output_file;
      char file_name[30];
      sprintf(file_name, "results/c_%ld.dat",iteration);
      unsigned int N[] = {nodeX};

      MPI_File_open(SUB_COMM, file_name, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &output_file);
      if (rank == 0) MPI_File_write(output_file, N, 1, MPI_UNSIGNED, MPI_STATUS_IGNORE);
      MPI_Offset displacement;
      displacement = rank*(thicknessMPI-nbAdditionalSlices)*nodeX*nodeY*sizeof(double) + sizeof(unsigned int); // Displacement in bytes
      // Set the view the current node has
      MPI_File_seek(output_file, displacement, MPI_SEEK_SET);
      MPI_File_write(output_file, concentration, nodeX*nodeY*thicknessMPI, MPI_DOUBLE, MPI_STATUS_IGNORE);
      MPI_File_close(&output_file);
    }

    free(concentrationSuiv);
    free(p);
    free(r);
    free(rsuiv);

    iteration+=1;
  }
  MPI_Finalize();
	return 0;
}
