#include "drainagebasin.hh"
#include <cmath>

DEM::DEM(double *my_x, int my_Mx, double *my_y, int my_My, double *my_z) {
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

DEM::~DEM() {
  gsl_interp_accel_free(x_accel);
  gsl_interp_accel_free(y_accel);
}

double DEM::x_min() {
  return x[0];
}

double DEM::x_max() {
  return x[Mx-1];
}


double DEM::y_min() {
  return y[0];
}

double DEM::y_max() {
  return y[My-1];
}

double DEM::dx() {
  return x[1] - x[0];
}

double DEM::dy() {
  return y[1] - y[0];
}

void DEM::indices(double x0, double y0, int &i, int &j) {
  i = gsl_interp_accel_find(x_accel, x, Mx, x0);
  j = gsl_interp_accel_find(y_accel, y, My, y0);
}

DEM_Bilinear::DEM_Bilinear(double *my_x, int my_Mx, double *my_y, int my_My, double *my_z)
  :DEM(my_x, my_Mx, my_y, my_My, my_z) {
}

double DEM_Bilinear::elevation(double x0, double y0) {
  int i, j;

  i = gsl_interp_accel_find(x_accel, x, Mx, x0);
  j = gsl_interp_accel_find(y_accel, y, My, y0);

  // Pretend that outside the grid the surface is perfectly flat:
  if (i < 0 ||
      i + 1 > Mx - 1 ||
      j < 0 ||
      j + 1 > My - 1) {

    return 0;
  }

  // Get the surface elevation at grid corners (arranged as follows):
  //
  // B-----C
  // |     |
  // |     |
  // A-----D
  //
  gsl_matrix_view z_view = gsl_matrix_view_array(z, Mx, My);
  gsl_matrix * m = &z_view.matrix;
  double
    A = gsl_matrix_get(m, i,     j),
    B = gsl_matrix_get(m, i,     j + 1),
    C = gsl_matrix_get(m, i + 1, j + 1),
    D = gsl_matrix_get(m, i + 1, j);

  double alpha = x0 - x[i],
    beta = y0 - y[j];

  return (1 - alpha)*(1 - beta)*A + (1 - alpha)*beta*B +
    alpha*beta*C + alpha*(1 - beta)*D;
}

void DEM_Bilinear::eval(double x0, double y0,
                        double f[], double jac[]) {
  int i, j;

  i = gsl_interp_accel_find(x_accel, x, Mx, x0);
  j = gsl_interp_accel_find(y_accel, y, My, y0);

  // Pretend that outside the grid the surface is perfectly flat:
  if (i < 0 ||
      i + 1 > Mx - 1 ||
      j < 0 ||
      j + 1 > My - 1) {

    if (f != NULL) {
      f[0] = f[1] = 0;
    }

    if (jac != NULL) {
      jac[0] = jac[1] = jac[2] = jac[3] = 0;
    }

    return;
  }

  // Get the surface elevation at grid corners (arranged as follows):
  //
  // B-----C
  // |     |
  // |     |
  // A-----D
  //
  gsl_matrix_view z_view = gsl_matrix_view_array(z, Mx, My);
  gsl_matrix * m = &z_view.matrix;
  double
    A = gsl_matrix_get(m, i,     j),
    B = gsl_matrix_get(m, i,     j + 1),
    C = gsl_matrix_get(m, i + 1, j + 1),
    D = gsl_matrix_get(m, i + 1, j);

  double gamma = one_over_dx * one_over_dy * (A + C - B - D);

  if (f != NULL) {
    f[0] = -((D - A) * one_over_dx + (y0 - y[j]) * gamma);
    f[1] = -((B - A) * one_over_dy + (x0 - x[i]) * gamma);
  }

  if (jac != NULL) {
    // J(i,j) = dfdy[i * dimension + j]
    // J(i,j) = \frac{\partial f_i}{\partial y_j}

    jac[0 * 2 + 0] = 0;
    jac[0 * 2 + 1] = -gamma;
    jac[1 * 2 + 0] = -gamma;
    jac[1 * 2 + 1] = 0;
  }
}
