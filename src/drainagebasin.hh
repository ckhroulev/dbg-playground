
// GSL stuff
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_math.h>       // GSL_NAN

#include <mpi.h>
#include <vector>

using namespace std;

struct Cell {
  int i, j;
};

// Purely virtual DEM interface class.
class DEM {
public:
  DEM(double *x, int Mx, double *y, int My, double *z);
  virtual ~DEM();

  virtual void eval(double x, double y, double f[], double jac[]) = 0;
  virtual double elevation(double x, double y) = 0;
  virtual void indices(double x, double y, int &i, int &j);

  virtual double x_min();
  virtual double x_max();

  virtual double y_min();
  virtual double y_max();

  virtual double dx();
  virtual double dy();

  double *x, *y, *z;
  int Mx, My;
  gsl_interp_accel *x_accel, *y_accel;
  double one_over_dx, one_over_dy;

  vector<Cell> path;
};

// A piecewise-bilinear DEM using surface elevation data on a regular grid.
class DEM_Bilinear : public DEM {
public:
  DEM_Bilinear(double *x, int Mx, double *y, int My, double *z);
  virtual ~DEM_Bilinear() {}
  virtual void eval(double x, double y, double f[], double jac[]);
  virtual double elevation(double x, double y);
};

int function(double t, const double y[], double f[], void* params);
int jacobian(double t, const double y[], double *dfdy, double dfdt[], void *params);

int read_dem(MPI_Comm com, int rank,
             const char *filename,
             vector<double> &X,
             vector<double> &Y,
             vector<double> &Z);

int write_mask(MPI_Comm com, int rank,
               const char *filename,
               vector<double> &X,
               vector<double> &Y,
               double *Z);

int flowline_gnuplot(gsl_odeiv2_system system, gsl_odeiv2_step *step,
                     double x0, double y0, const char *color);

int flowline_mask(gsl_odeiv2_system system, gsl_odeiv2_step *step,
                  double x0, double y0, int marker, bool &new_terminus, double *mask);
