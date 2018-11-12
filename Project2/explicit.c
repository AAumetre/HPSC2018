#include <stdio.h>
#include <math.h>
#include <stdlib.h>

int main(int argc, char const *argv[])
{
  double h,m, L, Tmax, vx, vy, vz, D;
  h = 0.1;
  m = 0.001;
  L = 1;
  Tmax = 1;
  vx = 0;
  vy = -1;
  vz = 0;
  D = 0.5;

  size_t nodeX = (int)(L/h) + 1;
  ///// !!!! vérifier pour les indices pairs le centre!
  size_t nodeY = nodeX, nodeZ =nodeX;
  printf("Number of nodes: %zu\n", nodeX);

  double initConcentration = 1; //g/m3

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

  size_t stopTime = Tmax/m;
  printf("Stop time: %zu\n", stopTime);

  for(size_t iteration = 0; iteration<= stopTime; iteration++)
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
           m * D * (concentrationPrev[i+1+j*nodeX+k*nodeX*nodeY]+concentrationPrev[i+(j+1)*nodeX+k*nodeX*nodeY]+
            concentrationPrev[i+j*nodeX+(k+1)*nodeX*nodeY]-6*concentrationPrev[i+j*nodeX+k*nodeX*nodeY]+
            concentrationPrev[i-1+j*nodeX+k*nodeX*nodeY]+concentrationPrev[i+(j-1)*nodeX+k*nodeX*nodeY]+
            concentrationPrev[i+j*nodeX+(k-1)*nodeX*nodeY])/pow(h,2) -
           m * vx * (concentrationPrev[i+1+j*nodeX+k*nodeX*nodeY]-concentrationPrev[i-1+j*nodeX+k*nodeX*nodeY])/(2*h) - 
           m * vy * (concentrationPrev[i+(j+1)*nodeX+k*nodeX*nodeY]-concentrationPrev[i+(j-1)*nodeX+k*nodeX*nodeY])/(2*h) - 
           m * vz * (concentrationPrev[i+j*nodeX+(k+1)*nodeX*nodeY]-concentrationPrev[i+j*nodeX+(k-1)*nodeX*nodeY])/(2*h);
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
  }




  return 0;
}