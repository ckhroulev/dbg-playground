#include "drainagebasin.hh"
#include <vector>
#include <cmath>

using namespace std;

int function(double t, const double y[], // inputs
             double f[],                 // output
             void* params) {

  ((DEM*)params)->evaluate(y, NULL, NULL, f, NULL);

  return GSL_SUCCESS;
}

int jacobian(double t, const double y[],  // inputs
             double *dfdy, double dfdt[], // outputs
             void *params) {

  ((DEM*)params)->evaluate(y, NULL, NULL, NULL, dfdy);

  dfdt[0] = 0;
  dfdt[1] = 0;

  return GSL_SUCCESS;
}

int streamline(gsl_odeiv2_system system,
               gsl_odeiv2_step *step,
               int i_start, int j_start,
               double *mask) {
  DEM *dem = (DEM*)system.params;

  gsl_matrix_view mask_view = gsl_matrix_view_array(mask, dem->get_Mx(), dem->get_My());
  gsl_matrix * m = &mask_view.matrix;

  // stop if the current cell already has a value assigned
  if (gsl_matrix_get(m, i_start, j_start) > 0)
    return 0;

  int steps_per_cell = 10.0,
    n_max = (dem->get_Mx() + dem->get_My()) * steps_per_cell,
    i_end, j_end, status;

  double grid_spacing = dem->dx() < dem->dy() ? dem->dx() : dem->dy(),
    position[2], err[2];

  position[0] = dem->get_x(i_start);
  position[1] = dem->get_y(j_start);

  double gradient[2], thickness, gradient_magnitude, step_size;
  dem->evaluate(position, NULL, &thickness, NULL, NULL);

  if (thickness < 1)
    return 0;

  for (int counter = 0; counter < n_max; ++counter) {

    if (position[0] <= dem->x_min() || position[0] >= dem->x_max() ||
        position[1] <= dem->y_min() || position[1] >= dem->y_max())
      break;

    dem->evaluate(position, NULL, &thickness, gradient, NULL);

    gradient_magnitude = sqrt(gradient[0]*gradient[0] + gradient[1]*gradient[1]);

    // Stop if the magnitude of the gradient is very low
    if (gradient_magnitude < 1e-16 || thickness < 1)
      break;

    step_size = (grid_spacing / steps_per_cell) / gradient_magnitude;

    // take a step
    status = gsl_odeiv2_step_apply(step,
                                   0,         // starting time (irrelevant)
                                   step_size, // step size
                                   position, err, NULL, NULL, &system);

    if (status != GSL_SUCCESS) {
      printf ("error, return value=%d\n", status);
      break;
    }
  }

  status = dem->find_cell(position, i_end, j_end);

  if (status != 0) {
    // the streamline exited the domain
    return 0;
  }

  double downstream_mask_value = gsl_matrix_get(m, i_end, j_end);

  if (downstream_mask_value > 0) {
    gsl_matrix_set(m, i_start, j_start,
                   downstream_mask_value);
  }

  return 0;
}
