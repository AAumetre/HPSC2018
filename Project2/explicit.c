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

  size_t nodeX = 3;//(int)(L/h) + 1;
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
      {concentration[i]=5;}//juste pour tester
    else if(isXbound==nodeX-1)
      {concentration[i]=5;}
    else if(inStage0<nodeX)
      {concentration[i]=5;}
    else if(inStage0>=nodeX*nodeY-nodeX-1)//vérifier les = !!!!!
      {concentration[i]=5;}
    else if(i<=nodeX*nodeY-1)
      {concentration[i]=5;}
    else if(i>=nodeX*nodeY*nodeZ-nodeX*nodeY)
      {concentration[i]=5;}
    else
      {//I am in the domain
        //mettre la moche formule
      }

    isXbound++;
    if(isXbound==nodeX)
        isXbound = 0;
  }

  for(int i=0; i<nodeX*nodeY*nodeZ ; i++)
    printf("%.0f ", concentration[i]);


  return 0;
}
