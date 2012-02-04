// GSL stuff
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>

#include <mpi.h>
#include <vector>

using namespace std;

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
               double *Z);

void init_mask(int Mx, int My,
               double *thickness,
               double *mask,
               double *tmp);

int streamline(gsl_odeiv_system system,
               gsl_odeiv_step *step,
               int i, int j,
               int steps_per_cell,
               int path_length,
               double min_elevation,
               double max_elevation,
               double *mask,
               double *output);
