
// GSL stuff
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_interp.h>

// A piecewise-bilinear DEM using surface elevation data on a regular grid.
class DEM_Bilinear {
public:
  DEM_Bilinear(double *x, int Mx, double *y, int My, double *z);
  ~DEM_Bilinear();
  void eval(double x, double y,
            double f[], double jac[]);
private:
  double one_over_dx, one_over_dy;
  double *x, *y, *z;
  int Mx, My;
  gsl_interp_accel *x_accel, *y_accel;
};

int function(double t, const double y[], double f[], void* params);
int jacobian(double t, const double y[], double *dfdy, double dfdt[], void *params);
