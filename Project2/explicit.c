#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "CSR_BSR.h"
#include <mpi.h>


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
Param readDat(char *filename)
{
    FILE *file = fopen(filename, "r");
    Param parameters;

    if (file != NULL)
    {
        fscanf(file, "%lf %lf %lf %lf %lf %lf %lf %lf %u %lf", &parameters.h, &parameters.m, &parameters.L, &parameters.Tmax, &parameters.vx, &parameters.vy, &parameters.vz, &parameters.D, &parameters.S, &parameters.rthreshold);
        printf("The parameters are : h %lf  m %lf  L %lf Tmax %lf %lf %lf %lf %lf %u %lf", parameters.h, parameters.m, parameters.L, parameters.Tmax, parameters.vx, parameters.vy, parameters.vz, parameters.D, parameters.S, parameters.rthreshold);

        fclose(file);
    }
    else
    {
        puts("File ERR0R !");
        exit(1);
    }

    return parameters;
}


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

  char bonusSlice=0;
  
  //if (nodeZ%rank) bonusSlice = 1;
  size_t thicknessMPI = (int)(nodeZ/world_size);
  if (rank == world_size-1) thicknessMPI++;

  size_t kCenter = floor(centerIndex/(nodeX*nodeY));
  size_t klocalCenter = kCenter - floor(world_size/2)*thicknessMPI;
  double *concentration = calloc(nodeX*nodeY*thicknessMPI, sizeof(double));
  double *c_ = calloc(nodeX*nodeY*thicknessMPI, sizeof(double));

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
    printf("interation %zu\n", iteration);
    //research for boundaries
    size_t isXbound = 0;
    size_t index = 0;
    if (rank == 0)
      index=nodeX*nodeY;
    size_t stopIndex = nodeX*nodeY*thicknessMPI;
    if (rank == world_size-1)
      stopIndex -=nodeY*nodeX;

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
      if(isXbound==nodeX)
          isXbound = 0;

    }

    //passer aux voisins
    //find klocal = 0 et = max
    int *commList = getCommListSlices(world_size);
    for (int commIndex=0 ; commIndex<4*(world_size-1) ; commIndex += 2) {
      // Get sender commList[commIndex] & receiver commList[commIndex+1]
      printf("Sender: %d to %d\n", commList[commIndex], commList[commIndex+1]);
        for (size_t index = 0; index < nodeX*nodeY; index ++)
        {
          int klocal = 0;
          int j = floor((index-klocal*nodeX*nodeY)/nodeX);
          int i = index - klocal * nodeX * nodeY - j * nodeX;
          MPI_Send(&concentration[i+j*nodeX+klocal*nodeX*nodeY], 1, MPI_DOUBLE, commList[commIndex+1], 0, MPI_COMM_WORLD);
          MPI_Recv(&number, 1, MPI_DOUBLE, commList[commIndex], 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

          /*klocal = thicknessMPI-1;//max
          MPI_Send(&concentration[i+j*nodeX+klocal*nodeX*nodeY], 1, MPI_DOUBLE, commList[commIndex+1], 0, MPI_COMM_WORLD);
          MPI_Recv(&number, 1, MPI_DOUBLE, commList[commIndex], 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
          */
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
    }*/
    ++iteration;
  }



MPI_Finalize();
  return 0;
}
