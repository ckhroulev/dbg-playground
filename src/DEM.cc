#include "drainagebasin.hh"

DEM::DEM(double *my_x, int my_Mx, double *my_y, int my_My,
         double *my_z) {

  x  = my_x;
  Mx = my_Mx;

  y  = my_y;
  My = my_My;

  z   = my_z;

  // We assume that the grid is uniform.
  x_spacing = x[1] - x[0];
  y_spacing = y[1] - y[0];

  one_over_dx = 1.0 / x_spacing;
  one_over_dy = 1.0 / y_spacing;
}

DEM::~DEM() {
}

void DEM::evaluate(const double *position, double *elevation, double *f) {
  int ierr, i, j;
  double A, B, C, D, delta_x, delta_y;

  ierr = this->find_cell(position, i, j);

  // Pretend that outside the grid the surface is perfectly flat, the elevation
  // of the sea level (0) and ice-free.
  if (ierr != 0) {

    if (elevation != NULL)
      *elevation = 0;

    if (f != NULL)
      f[0] = f[1] = 0;

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

} // end of DEM::evaluate()
