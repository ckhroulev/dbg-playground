
// GSL stuff
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_interp.h>

// Purely virtual DEM interface class.
class DEM {
public:
  DEM(double *x, int Mx, double *y, int My, double *z);
  virtual ~DEM();

  virtual void eval(double x, double y, double f[], double jac[]) = 0;
protected:
  double *x, *y, *z;
  int Mx, My;
  gsl_interp_accel *x_accel, *y_accel;
  double one_over_dx, one_over_dy;
};

// A piecewise-bilinear DEM using surface elevation data on a regular grid.
class DEM_Bilinear : public DEM {
public:
  DEM_Bilinear(double *x, int Mx, double *y, int My, double *z);
  virtual ~DEM_Bilinear() {}
  virtual void eval(double x, double y, double f[], double jac[]);
};

int function(double t, const double y[], double f[], void* params);
int jacobian(double t, const double y[], double *dfdy, double dfdt[], void *params);

int flowline(gsl_odeiv2_system system, gsl_odeiv2_step *step,
             double x0, double y0, const char *color);
