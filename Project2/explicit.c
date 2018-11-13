#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "CSR_BSR.h"

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

  if (argc != 2)
  {
    printf("Wrong arguments\n");
    printf("Please use the function with ./exe param.dat\n");
    return 1;
  }

  printf("Reading file %s ...\n", argv[1]);
  Param parameters = readDat(argv[1]);

  size_t nodeX = (int)(parameters.L/parameters.h) + 1;
  ///// !!!! vérifier pour les indices pairs le centre!
  size_t nodeY = nodeX, nodeZ =nodeX;
  printf("Number of nodes: %zu\n", nodeX);

  double initConcentration = 1; //[g/m3]

  double *concentration = calloc(nodeX*nodeY*nodeZ, sizeof(double));

  if (concentration == NULL) {
    puts("Mem ERR0R !");
    exit(1);
  }

  double *concentrationPrev = calloc(nodeX*nodeY*nodeZ, sizeof(double));

  if (concentrationPrev == NULL) {
    puts("Mem ERR0R !");
    exit(1);
  }

  size_t centerIndex = nodeX*nodeY*floor(nodeZ/2)+ floor(nodeY/2)*nodeX + floor(nodeX/2);
  printf("Index of center: %zu\n", centerIndex);

  concentrationPrev[centerIndex] = initConcentration;
  //for(int i=0; i<nodeX*nodeY*nodeZ ; i++)
    //printf("%f ", concentration[i]);

  size_t stopTime = parameters.Tmax/parameters.m;
  printf("Stop time: %zu\n", stopTime);
  size_t iteration = 0;
  bool onBoundary = false;
  bool valueOnBoundary= false;

  while (iteration <= stopTime && !valueOnBoundary)
  {
    printf("interation %zu\n", iteration);
    //research for boundaries
    size_t isXbound = 0;
    for(size_t index = 0; index<nodeX*nodeY*nodeZ; index++)
    {
      int stage = floor(index/(nodeX*nodeY));
      int inStage0 = index-stage*nodeX*nodeY;

      if(isXbound==0)
        {}
      else if(isXbound==nodeX-1)
        {}
      else if(inStage0<nodeX)
        {}
      else if(inStage0>=nodeX*nodeY-nodeX-1)//vérifier les = !!!!!
        {}
      else if(index<=nodeX*nodeY-1)
        {}
      else if(index>=nodeX*nodeY*nodeZ-nodeX*nodeY)
        {}
      else
        {//I am in the domain
          int k = floor(index/(nodeX*nodeY));
          int j = floor((index-k*nodeX*nodeY)/nodeX);
          int i = index - k * nodeX * nodeY - j * nodeX;
          concentration[i+j*nodeX+k*nodeX*nodeY] = concentrationPrev[i+j*nodeX+k*nodeX*nodeY] +
           parameters.m * parameters.D * (concentrationPrev[i+1+j*nodeX+k*nodeX*nodeY]+concentrationPrev[i+(j+1)*nodeX+k*nodeX*nodeY]+
            concentrationPrev[i+j*nodeX+(k+1)*nodeX*nodeY]-6*concentrationPrev[i+j*nodeX+k*nodeX*nodeY]+
            concentrationPrev[i-1+j*nodeX+k*nodeX*nodeY]+concentrationPrev[i+(j-1)*nodeX+k*nodeX*nodeY]+
            concentrationPrev[i+j*nodeX+(k-1)*nodeX*nodeY])/pow(parameters.h,2) -
           parameters.m * parameters.vx * (concentrationPrev[i+1+j*nodeX+k*nodeX*nodeY]-concentrationPrev[i-1+j*nodeX+k*nodeX*nodeY])/(2*parameters.h) -
           parameters.m * parameters.vy * (concentrationPrev[i+(j+1)*nodeX+k*nodeX*nodeY]-concentrationPrev[i+(j-1)*nodeX+k*nodeX*nodeY])/(2*parameters.h) -
           parameters.m * parameters.vz * (concentrationPrev[i+j*nodeX+(k+1)*nodeX*nodeY]-concentrationPrev[i+j*nodeX+(k-1)*nodeX*nodeY])/(2*parameters.h);

          onBoundary = ((isXbound == nodeX-2) || (isXbound == 1) || (inStage0 >= nodeX && inStage0 <= 2*nodeX-2) || (inStage0 >= nodeY*nodeY-2*nodeX-1) || (index<=2*nodeX*nodeY-2) || (index >= nodeX*nodeY*nodeZ - 2*nodeX*nodeY));
          if (onBoundary  && concentration[i+j*nodeX+k*nodeX*nodeY] != 0){
            printf("coucou");
            valueOnBoundary=true;
          }

        }

      isXbound++;
      if(isXbound==nodeX)
          isXbound = 0;
    }

    if (!(iteration%500))
    {
      for(int i=0; i<nodeX*nodeY*nodeZ ; i++)
      {
        concentrationPrev[i] = concentration[i];
        printf("%f ", concentration[i]);
      }

      printf("\n \n \n");
    }
    ++iteration;
  }




  return 0;
}
