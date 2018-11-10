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
  for(size_t i = 0; i<nodeX*nodeY*nodeZ; i++)
  {
    int stage = floor(i/(nodeX*nodeY));
    int inStage0 = i-stage*nodeX*nodeY;

    if(isXbound==0)
      {}
    else if(isXbound==nodeX-1)
      {}
    else if(inStage0<nodeX)
      {}
    else if(inStage0>=nodeX*nodeY-nodeX-1)//vérifier les = !!!!!
      {}
    else if(i<=nodeX*nodeY-1)
      {}
    else if(i>=nodeX*nodeY*nodeZ-nodeX*nodeY)
      {}
    else
      {//I am in the domain
        concentration[] = concentrationPrev[] + m * D * (concentrationPrev[i+1]+concentrationPrev[j+1]+concentrationPrev[k+1]-6*concentrationPrev[]+concentrationPrev[i-1]+concentrationPrev[j-1]+concentrationPrev[k-1])/pow(h,2) - m * vx * (concentrationPrev[i+1]-concentrationPrev[i-1])/(2*h) - m * vy * (concentrationPrev[j+1]-concentrationPrev[j-1])/(2*h) - m * vz * (concentrationPrev[k+1]-concentrationPrev[k-1])/(2*h);
      }

    isXbound++;
    if(isXbound==nodeX)
        isXbound = 0;
  }

  for(int i=0; i<nodeX*nodeY*nodeZ ; i++)
    printf("%.0f ", concentration[i]);


  return 0;
}
