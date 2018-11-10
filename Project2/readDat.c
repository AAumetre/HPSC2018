#include <stdlib.h>
#include <stdio.h>

#define BUFFER_SIZE 1000

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
    char* read = malloc(BUFFER_SIZE*sizeof(char));
    if (read == NULL)
    {
        puts("Mem ERR0R !");
        exit(1);
    }
    Param parameters;

    if (file != NULL)
    {
        fscanf(file, "%lf %lf %lf %lf %lf %lf %lf %lf %u %lf", &parameters.h, &parameters.m, &parameters.L, &parameters.Tmax, &parameters.vx, &parameters.vy, &parameters.vz, &parameters.D, &parameters.S, &parameters.rthreshold);
        printf("The parameters are : %lf %lf %lf %lf %lf %lf %lf %lf %u %lf", parameters.h, parameters.m, parameters.L, parameters.Tmax, parameters.vx, parameters.vy, parameters.vz, parameters.D, parameters.S, parameters.rthreshold);

        fclose(file);
    }
    else
    {
        puts("File ERR0R !");
        exit(1);
    }

    free(read);

    return parameters;
}


/* ----------------------------------------------------------------------------*
*  MAIN
* -----------------------------------------------------------------------------*
* @ usage:    ./exe param.ppm
* -----------------------------------------------------------------------------*
* Main function reading a dat file
* ----------------------------------------------------------------------------*/

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

    return 0;
}
