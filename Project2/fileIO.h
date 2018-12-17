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
Param readDat(char *filename);


void writeSequentialFile(double values[], int iteration, int N);
