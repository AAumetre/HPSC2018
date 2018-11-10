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
  printf("%zu\n", nodeX);

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
  printf("%zu\n", centerIndex);

  concentration[centerIndex] = initConcentration;
  //for(int i=0; i<nodeX*nodeY*nodeZ ; i++)
    //printf("%f ", concentration[i]);

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
        concentration[i+j*nodeX+k*nodeX*nodeY] = concentrationPrev[i+j*nodeX+k*nodeX*nodeY] + m * D * (concentrationPrev[i+1+j*nodeX+k*nodeX*nodeY]+concentrationPrev[i+(j+1)*nodeX+k*nodeX*nodeY]+concentrationPrev[i+j*nodeX+(k+1)*nodeX*nodeY]-6*concentrationPrev[i+j*nodeX+k*nodeX*nodeY]+concentrationPrev[i-1+j*nodeX+k*nodeX*nodeY]+concentrationPrev[i+(j-1)*nodeX+k*nodeX*nodeY]+concentrationPrev[i+j*nodeX+(k-1)*nodeX*nodeY])/pow(h,2) - m * vx * (concentrationPrev[i+1+j*nodeX+k*nodeX*nodeY]-concentrationPrev[i-1+j*nodeX+k*nodeX*nodeY])/(2*h) - m * vy * (concentrationPrev[i+(j+1)*nodeX+k*nodeX*nodeY]-concentrationPrev[i+(j-1)*nodeX+k*nodeX*nodeY])/(2*h) - m * vz * (concentrationPrev[i+j*nodeX+(k+1)*nodeX*nodeY]-concentrationPrev[i+j*nodeX+(k-1)*nodeX*nodeY])/(2*h);
      }

    isXbound++;
    if(isXbound==nodeX)
        isXbound = 0;
  }

  for(int i=0; i<nodeX*nodeY*nodeZ ; i++)
    printf("%f ", concentration[i]);


  return 0;
}
