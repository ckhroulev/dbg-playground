// GSL stuff
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>

// MPI and C++ stuff
#include <mpi.h>
#include <vector>
using namespace std;

#include "Array2D.hh"

enum MASK_VALUES {NO_VALUE = -2, ICE_FREE = -1};

int function(double t, const double y[], double f[], void* params);

int read_dem(MPI_Comm com, int rank,
             const char *filename,
             vector<double> &X,
             vector<double> &Y,
             vector<double> &Z,
             vector<double> &thk);

int write_mask(MPI_Comm com, int rank,
               const char *filename,
               vector<double> &X,
               vector<double> &Y,
               Array2D<double> &Z);

void init_mask(double *thickness,
               Array2D<double> &mask,
               Array2D<double> &tmp);

int streamline(gsl_odeiv_system system,
               gsl_odeiv_step *step,
               int i, int j,
               int steps_per_cell,
               int path_length,
               double min_elevation,
               double max_elevation,
               Array2D<double> &mask,
               Array2D<double> &output);
