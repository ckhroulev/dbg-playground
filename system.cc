#include "drainagebasin.hh"

DEM_Bilinear::DEM_Bilinear(double *my_x, int my_Mx, double *my_y, int my_My, double *my_z) {
  x_accel = gsl_interp_accel_alloc();
  y_accel = gsl_interp_accel_alloc();

  x = my_x;
  Mx = my_Mx;

  y = my_y;
  My = my_My;

  z = my_z;

  // We assume that the grid is equally-spaced.
  one_over_dx = 1.0 / (x[1] - x[0]);
  one_over_dy = 1.0 / (y[1] - y[0]);
}

DEM_Bilinear::~DEM_Bilinear() {
  gsl_interp_accel_free(x_accel);
  gsl_interp_accel_free(y_accel);
}

void DEM_Bilinear::eval(double x0, double y0,
                        double f[], double jac[]) {
  int i, j;

  if (x[0] > x0 || x0 >= x[Mx-1] ||
      y[0] > y0 || y0 >= y[My-1]) {

    // printf("x = %f, y = %f\n", x0, y0);

    if (f != NULL) {
      f[0] = f[1] = 0;
    }

    if (jac != NULL) {
      jac[0] = jac[1] = jac[2] = jac[3] = 0;
    }

    return;
  }

  i = gsl_interp_accel_find(x_accel, x, Mx, x0);
  j = gsl_interp_accel_find(y_accel, y, My, y0);

  gsl_matrix_view z_view = gsl_matrix_view_array(z, My, Mx);
  gsl_matrix * m = &z_view.matrix;
  double
    A = gsl_matrix_get(m, i,     j),
    B = gsl_matrix_get(m, i,     j + 1),
    C = gsl_matrix_get(m, i + 1, j + 1),
    D = gsl_matrix_get(m, i + 1, j);

  if (f != NULL) {
    double
      alpha = one_over_dx * (x[i] - x0),
      beta  = one_over_dy * (y[j] - y0);

    f[0] = - one_over_dx * ((1 - beta)  * (D - A) + beta  * (C - B));
    f[1] = - one_over_dy * ((1 - alpha) * (B - A) + alpha * (C - D));
  }

  if (jac != NULL) {
    // J(i,j) = dfdy[i * dimension + j]
    // J(i,j) = \frac{\partial f_i}{\partial y_j}

    double jac_factor = one_over_dx * one_over_dy * (A + C - B - D);
    jac[0 * 2 + 0] = 0;
    jac[0 * 2 + 1] = jac_factor;
    jac[1 * 2 + 0] = jac_factor;
    jac[1 * 2 + 1] = 0;
  }
}

int function(double t, const double y[], // inputs
             double f[],                 // output
             void* params) {             // optional input/output

  DEM_Bilinear *dem = (DEM_Bilinear*)params;
  dem->eval(y[0], y[1], f, NULL);

  return GSL_SUCCESS;
}

int jacobian(double t, const double y[],  // inputs
             double *dfdy, double dfdt[], // outputs
             void *params) {              // optional input/output

  DEM_Bilinear *dem = (DEM_Bilinear*)params;
  dem->eval(y[0], y[1], NULL, dfdy);

  dfdt[0] = 0;
  dfdt[1] = 0;

  return GSL_SUCCESS;
}
