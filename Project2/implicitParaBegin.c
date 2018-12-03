#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define initConcentration 1 //[g/m3]

/* ----------------------------------------------------------------------------*
* Structure: Param
* -----------------------------------------------------------------------------*
* @ param:  h             double          spatial step [m]
*           m             double          time step [s]
*           L             double          space domain size [m]
*           Tmax          double          time domain size [s]
*           vx            double          horizontal velocity [m/s]
*           vy            double          vertical velocity [m/s]
*           vz            double          z-axis velocity [m/s]
*           D             double          diffusion coefficient [m^2/s]
*           S             unsigned int
*           rthreshold    double
* -----------------------------------------------------------------------------*
* Structure the parameters defining the domain
* ----------------------------------------------------------------------------*/
typedef struct Param{
  double h,m, L, Tmax, vx, vy, vz, D, rthreshold;
  unsigned int S;
} Param;

/* ----------------------------------------------------------------------------*
* function: readDat
* -----------------------------------------------------------------------------*
* @ param:  filename   char*            pointer on the dat file path to open
* @ return: Param      parameters       structure containing the parameters
* -----------------------------------------------------------------------------*
* Function reading a given dat file and returning a structure containing its
* data
* ----------------------------------------------------------------------------*/
Param readDat(char *filename){
	FILE *file = fopen(filename, "r");
	Param parameters;

	if (file != NULL)
	{
		fscanf(file, "%lf %lf %lf %lf %lf %lf %lf %lf %u %lf", &parameters.h, &parameters.m, &parameters.L, &parameters.Tmax, &parameters.vx, &parameters.vy, &parameters.vz, &parameters.D, &parameters.S, &parameters.rthreshold);
		//printf("The parameters are : h %lf  m %lf  L %lf Tmax %lf %lf %lf %lf %lf %u %lf", parameters.h, parameters.m, parameters.L, parameters.Tmax, parameters.vx, parameters.vy, parameters.vz, parameters.D, parameters.S, parameters.rthreshold);

		fclose(file);
	}
	else
	{
		puts("File ERR0R !");
		exit(1);
	}

	return parameters;
}

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
void Ap(double* p, double* Apresult, size_t nodeX, size_t nodeY, size_t nodeZ, double h, double m, double vx, double vy, double vz, double D)
{
  for(size_t index = 0; index< nodeX*nodeY*nodeZ; index++)
  {
    int k = floor(index/(nodeX*nodeY));
    int j = floor((index-k*nodeX*nodeY)/nodeX);
    int i = index - k * nodeX * nodeY - j * nodeX;
    double Ddivh2 = D/(h*h);
    if(!(i==0 || i == nodeX-1 || j==0 || j == nodeY-1 || k==0 || k==nodeZ-1))
    {
      //Apresult[i+j*nodeX+k*nodeX*nodeY] = p[i+j*nodeX+k*nodeX*nodeY]*(1/m + 6*Ddivh2) + p[i-1+j*nodeX+k*nodeX*nodeY]*(-vx/(2*h) - Ddivh2) + p[i+1+j*nodeX+k*nodeX*nodeY]*(vx/(2*h) - Ddivh2) + p[i+(j-1)*nodeX+k*nodeX*nodeY]*(-vy/(2*h) - Ddivh2) + p[i+(j+1)*nodeX+k*nodeX*nodeY]*(vy/(2*h) - Ddivh2) + p[i+j*nodeX+(k-1)*nodeX*nodeY]*(-vz/(2*h) -Ddivh2) + p[i+j*nodeX+(k+1)*nodeX*nodeY]*(vz/(2*h) - Ddivh2);
      Apresult[i+j*nodeX+k*nodeX*nodeY] = p[i+j*nodeX+k*nodeX*nodeY]*(1/m + 6*Ddivh2) + p[i-1+j*nodeX+k*nodeX*nodeY]*(-vx/(2*h) - Ddivh2) + p[i+1+j*nodeX+k*nodeX*nodeY]*(vx/(2*h) - Ddivh2) + p[i+(j-1)*nodeX+k*nodeX*nodeY]*(-vy/(2*h) - Ddivh2) + p[i+(j+1)*nodeX+k*nodeX*nodeY]*(vy/(2*h) - Ddivh2) + p[i+j*nodeX+(k-1)*nodeX*nodeY]*(-vz/(2*h) -Ddivh2) + p[i+j*nodeX+(k+1)*nodeX*nodeY]*(vz/(2*h) - Ddivh2);
      Apresult[i+j*nodeX+k*nodeX*nodeY] *= m;
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

  if (concentration == NULL) {
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

    while(sqrt(SquaredNorm(r, nodeX*nodeY*thicknessMPI))>=parameters.rthreshold)
    {
      // alpha = r^T*r / p^T*A*p
      double *Apvect = calloc(nodeX*nodeY*thicknessMPI, sizeof(double));

      if (Apvect == NULL)
      {
        puts("Mem ERR0R Apvect!");
        exit(1);
      }
      Ap(p, Apvect, nodeX, nodeY, thicknessMPI, parameters.h, parameters.m, parameters.vx, parameters.vy, parameters.vz, parameters.D);
      double alpha = SquaredNorm(r, nodeX*nodeY*thicknessMPI) / VectorProduct(p, Apvect, nodeX*nodeY*thicknessMPI);
      // xi+1 = xi + alpha*p
      SumVect(concentrationSuiv, concentrationSuiv, p, alpha, nodeX*nodeY*thicknessMPI);
      // ri+1 = ri - alpha*Ap
      SumVect(rsuiv, r, Apvect, -alpha, nodeX*nodeY*thicknessMPI);
      // beta = ri+1^T*ri+1 / ri^T*ri
      double beta = SquaredNorm(rsuiv, nodeX*nodeY*thicknessMPI)/SquaredNorm(r,  nodeX*nodeY*thicknessMPI);
      // pi+1 = ri+1 + beta*p
      SumVect(p, rsuiv, p, beta, nodeX*nodeY*thicknessMPI);
      for (size_t copyIndex = 0; copyIndex < nodeX*nodeY*thicknessMPI; copyIndex++)
      {
          r[copyIndex] = rsuiv[copyIndex];
      }
      free(Apvect);
    }

    for(size_t index = 0; index< nodeX*nodeY*thicknessMPI; index++)
    {
      int k = floor(index/(nodeX*nodeY));
      int j = floor((index-k*nodeX*nodeY)/nodeX);
      int i = index - k * nodeX * nodeY - j * nodeX;
      if((i<=1 || i >= nodeX-2 || j<=1 || j >= nodeY-2 || k<=1 || k<=nodeZ-2) && concentrationSuiv[index]> 100*parameters.rthreshold)
      {
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

    for(size_t i=0; i<nodeX*nodeY*thicknessMPI; i++)
    {
      if(concentrationSuiv[i] > 1e-12)
        printf("%zu %lf ,", i, concentrationSuiv[i]);
    }
    printf("\n");

    free(concentrationSuiv);
    free(p);
    free(r);
    free(rsuiv);

    iteration+=1;
  }

}
