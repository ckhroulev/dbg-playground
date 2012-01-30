// GSL stuff
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

class Context {
  Context(double *x, int Mx, double *y, int My);
  ~Context();
public:
  double one_over_dx, one_over_dy;
  double *x, double *y;
  int Mx, My;
  gsl_interp_accel *x_accel, *y_accel;
};

Context::Context(double *my_x, int my_Mx, double *my_y, int my_My) {
  x_accel = gsl_interp_accel_alloc();
  y_accel = gsl_interp_accel_alloc();

  x = my_x;
  Mx = my_Mx;
  one_over_dx = 1.0 / (x[1] - x[0]);

  y = my_y;
  My = my_My;
  one_over_dy = 1.0 / (y[1] - y[0]);
}

Context::~Context() {
  gsl_interp_accel_free(x_accel);
  gsl_interp_accel_free(y_accel);
}

int function(double t, const double y[], // inputs
             double f[],                 // output
             void* params) {             // optional input/output

  f[0] = 1.0;
  f[1] = 1.0;

  return GSL_SUCCESS;
}

int jacobian(double t, const double y[],  // inputs
             double *dfdy, double dfdt[], // outputs
             void *params) {              // optional input/output
  // J(i,j) = dfdy[i * dimension + j]
  // J(i,j) = \frac{\partial f_i}{\partial y_j}

  dfdy[0 * 2 + 0] = 0;
  dfdy[0 * 2 + 1] = 0;
  dfdy[1 * 2 + 0] = 0;
  dfdy[1 * 2 + 1] = 0;

  dfdt[0] = 0;
  dfdt[1] = 0;

  return GSL_SUCCESS;
}
