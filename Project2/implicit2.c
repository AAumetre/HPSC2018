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
      Apresult[i+j*nodeX+k*nodeX*nodeY] = p[i+j*nodeX+k*nodeX*nodeY]*(1/m + 6*Ddivh2) + p[i-1+j*nodeX+k*nodeX*nodeY]*(-vx/(2*h) - Ddivh2) + p[i+1+j*nodeX+k*nodeX*nodeY]*(vx/(2*h) - Ddivh2) + p[i+(j-1)*nodeX+k*nodeX*nodeY]*(-vy/(2*h) - Ddivh2) + p[i+(j+1)*nodeX+k*nodeX*nodeY]*(vy/(2*h) - Ddivh2) + p[i+j*nodeX+(k-1)*nodeX*nodeY]*(-vz/(2*h) -Ddivh2) + p[i+j*nodeX+(k+1)*nodeX*nodeY]*(vz/(2*h) - Ddivh2);
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

  // Retrieving data from the .dat file
	Param parameters = readDat(argv[1]);
	size_t nodeX = (int)(parameters.L/parameters.h) + 1;
	size_t nodeY = nodeX, nodeZ =nodeX;
	size_t stopTime = 2;//parameters.Tmax/parameters.m;

  double *concentration = calloc(nodeX*nodeY*nodeZ, sizeof(double));

  if (concentration == NULL) {
    puts("Mem ERR0R !");
    exit(1);
  }

  size_t centerIndex = nodeX*nodeY*floor(nodeZ/2)+ floor(nodeY/2)*nodeX + floor(nodeX/2);
  printf("Index of center: %zu\n", centerIndex);

  concentration[centerIndex] = initConcentration;

  printf("Avant iter \n");
  for(size_t i=0; i<nodeX*nodeY*nodeZ; i++)
  {
    if(concentration[i] != 0)
      printf("%zu %lf ,", i, concentration[i]);
  }
  printf("\n");

  size_t iteration = 0;
  while (iteration <= stopTime)
  {
    //--------------------------------------------------------------------------
    //              Conjugate Gradient Method
    //--------------------------------------------------------------------------
    // x0 = 0
    double *concentrationSuiv = calloc(nodeX*nodeY*nodeZ, sizeof(double));

    // r0 = b - Ax0 = b
    double *r = calloc(nodeX*nodeY*nodeZ, sizeof(double));
    double *rsuiv = calloc(nodeX*nodeY*nodeZ, sizeof(double));

    if (r == NULL || rsuiv == NULL || concentrationSuiv == NULL)
    {
      puts("Mem ERR0R concentrationSuiv or r or rsuiv!");
      exit(1);
    }

    for (size_t copyIndex = 0; copyIndex < nodeX*nodeY*nodeZ; copyIndex++)
    {
        r[copyIndex] = concentration[copyIndex];
    }

    // p0 = r0
    double *p = calloc(nodeX*nodeY*nodeZ, sizeof(double));

    if (p == NULL)
    {
      puts("Mem ERR0R p!");
      exit(1);
    }
    for (size_t copyIndex = 0; copyIndex < nodeX*nodeY*nodeZ; copyIndex++)
    {
        p[copyIndex] = r[copyIndex];
    }

    while(sqrt(SquaredNorm(r, nodeX*nodeY*nodeZ))>parameters.rthreshold)
    {
      // alpha = r^T*r / p^T*A*p
      double *Apvect = calloc(nodeX*nodeY*nodeZ, sizeof(double));

      if (Apvect == NULL)
      {
        puts("Mem ERR0R Apvect!");
        exit(1);
      }
      Ap(p, Apvect, nodeX, nodeY, nodeZ, parameters.h, parameters.m, parameters.vx, parameters.vy, parameters.vz, parameters.D);
      double alpha = SquaredNorm(r, nodeX*nodeY*nodeZ) / VectorProduct(p, Apvect, nodeX*nodeY*nodeZ);
      // xi+1 = xi + alpha*p
      SumVect(concentrationSuiv, concentrationSuiv, p, alpha, nodeX*nodeY*nodeZ);
      // ri+1 = ri - alpha*Ap
      SumVect(rsuiv, r, Apvect, -alpha, nodeX*nodeY*nodeZ);
      // beta = ri+1^T*ri+1 / ri^T*ri
      double beta = SquaredNorm(rsuiv, nodeX*nodeY*nodeZ)/SquaredNorm(r,  nodeX*nodeY*nodeZ);
      // pi+1 = ri+1 + beta*p
      SumVect(p, rsuiv, p, beta, nodeX*nodeY*nodeZ);
      for (size_t copyIndex = 0; copyIndex < nodeX*nodeY*nodeZ; copyIndex++)
      {
          r[copyIndex] = rsuiv[copyIndex];
      }

      for(size_t i=0; i<nodeX*nodeY*nodeZ; i++)
      {
        if(concentrationSuiv[i] > 1e-12)
          printf("%zu %lf ,", i, concentrationSuiv[i]);
      }
      printf("\n\n");
      free(Apvect);
    }

    for(size_t index = 0; index< nodeX*nodeY*nodeZ; index++)
    {
      int k = floor(index/(nodeX*nodeY));
      int j = floor((index-k*nodeX*nodeY)/nodeX);
      int i = index - k * nodeX * nodeY - j * nodeX;
      if(!(i<=1 || i >= nodeX-2 || j<=1 || j >= nodeY-2 || k<=1 || k<=nodeZ-2) && concentrationSuiv[index]> 1e-12)
      {
          printf("STOP\n");
          break;
      }
    }

    for (size_t copyIndex = 0; copyIndex < nodeX*nodeY*nodeZ; copyIndex++)
    {
        concentration[copyIndex] = concentrationSuiv[copyIndex];
    }

    printf("Iteration %zu \n", iteration);
    for(size_t i=0; i<nodeX*nodeY*nodeZ; i++)
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
