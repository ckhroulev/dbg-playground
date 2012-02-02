
// GSL stuff
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_math.h>       // GSL_NAN

#include <mpi.h>
#include <vector>

using namespace std;

class DEM {
public:
  DEM(double *x, int Mx, double *y, int My, double *z, double *thickness);
  ~DEM();

  void evaluate(const double *position, double *elevation, double *thickness,
                double *f, double *jac);

  int find_cell(const double *position,
                int &i, int &j);

  double x_min();
  double x_max();

  double y_min();
  double y_max();

  double dx();
  double dy();

  double get_x(int i);
  double get_y(int j);

  int get_Mx();
  int get_My();

protected:
  double *x, *y, *z, *thk;
  int Mx, My;
  gsl_interp_accel *x_accel, *y_accel;
  double x_spacing, y_spacing, one_over_dx, one_over_dy;
  int i_last, j_last;

  void get_corner_values(int i, int j, double *data,
                         double &A, double &B, double &C, double &D);
};

int function(double t, const double y[], double f[], void* params);
int jacobian(double t, const double y[], double *dfdy, double dfdt[], void *params);

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
               double *elevation,
               double *thickness,
               double *mask,
               double *tmp);

int streamline(gsl_odeiv2_system system,
               gsl_odeiv2_step *step,
               int i, int j,
               double *mask,
               double *output);
