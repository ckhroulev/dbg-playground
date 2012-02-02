#include "drainagebasin.hh"

DEM::DEM(double *my_x, int my_Mx, double *my_y, int my_My,
         double *my_z, double *my_thk) {

  x  = my_x;
  Mx = my_Mx;

  y  = my_y;
  My = my_My;

  z   = my_z;
  thk = my_thk;

  // We assume that the grid is uniform.
  x_spacing = x[1] - x[0];
  y_spacing = y[1] - y[0];

  one_over_dx = 1.0 / x_spacing;
  one_over_dy = 1.0 / y_spacing;
}

DEM::~DEM() {
}

double DEM::get_x(int i) {
  return x[i];
}

double DEM::get_y(int j) {
  return y[j];
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
  return x_spacing;
}

double DEM::dy() {
  return y_spacing;
}

int DEM::get_Mx() {
  return Mx;
}

int DEM::get_My() {
  return My;
}

void DEM::evaluate(const double *position, double *elevation, double *thickness,
                   double *f, double *jac) {
  int ierr, i, j;
  double A, B, C, D, delta_x, delta_y;

  ierr = this->find_cell(position, i, j);

  // Pretend that outside the grid the surface is perfectly flat, the elevation
  // of the sea level (0) and ice-free.
  if (ierr != 0) {

    if (elevation != NULL)
      *elevation = 0;

    if (thickness != NULL)
      *thickness = 0;

    if (f != NULL)
      f[0] = f[1] = 0;

    if (jac != NULL)
      jac[0] = jac[1] = jac[2] = jac[3] = 0;

    return;
  }

  delta_x = position[0] - x[i];
  delta_y = position[1] - y[j];

  double
    alpha = one_over_dx * delta_x,
    beta  = one_over_dy * delta_y;

  this->get_corner_values(i, j, z, A, B, C, D);

  double gamma = one_over_dx * one_over_dy * (A + C - B - D);

  // surface elevation
  if (elevation != NULL) {
    *elevation = ( (1 - alpha) * (1 - beta) * A +
                   (1 - alpha) *      beta  * B +
                   alpha       *      beta  * C +
                   alpha       * (1 - beta) * D );
  }

  // minus the gradient
  if (f != NULL) {
    f[0] = -((D - A) * one_over_dx + delta_x * gamma);
    f[1] = -((B - A) * one_over_dy + delta_y * gamma);
  }

  // the jacobian of the above
  if (jac != NULL) {
    jac[0 * 2 + 0] = 0;
    jac[0 * 2 + 1] = -gamma;
    jac[1 * 2 + 0] = -gamma;
    jac[1 * 2 + 1] = 0;
  }

  // ice thickness
  if (thickness != NULL) {
    this->get_corner_values(i, j, thk, A, B, C, D);

    *thickness = ( (1 - alpha) * (1 - beta) * A +
                   (1 - alpha) *      beta  * B +
                   alpha       *      beta  * C +
                   alpha       * (1 - beta) * D );
  }

} // end of DEM::evaluate()
