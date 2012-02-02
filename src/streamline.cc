#include "drainagebasin.hh"
#include <vector>
#include <cmath>
#include <deque>

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
               double *old_mask, double *new_mask) {
  DEM *dem = (DEM*)system.params;

  deque<int> is, js;

  gsl_matrix_view old_mask_view = gsl_matrix_view_array(old_mask, dem->get_Mx(), dem->get_My());
  gsl_matrix * m = &old_mask_view.matrix;

  gsl_matrix_view new_mask_view = gsl_matrix_view_array(new_mask, dem->get_Mx(), dem->get_My());
  gsl_matrix * o = &new_mask_view.matrix;

  // stop if the current cell already has a value assigned
  if (gsl_matrix_get(m, i_start, j_start) > 0)
    return 0;

  int steps_per_cell = 10.0,
    n_max = (dem->get_Mx() + dem->get_My()) * steps_per_cell,
    i_end = 0, j_end = 0, status;

  double grid_spacing = dem->dx() < dem->dy() ? dem->dx() : dem->dy(),
    position[2], err[2];

  position[0] = dem->get_x(i_start) + dem->dx() * 0.5;
  position[1] = dem->get_y(j_start) + dem->dy() * 0.5;

  double gradient[2], thickness, gradient_magnitude, step_size;
  dem->evaluate(position, NULL, &thickness, NULL, NULL);

  if (thickness < 10)
    return 0;

  for (int counter = 0; counter < n_max; ++counter) {

    if (position[0] <= dem->x_min() || position[0] >= dem->x_max() ||
        position[1] <= dem->y_min() || position[1] >= dem->y_max())
      break;

    dem->evaluate(position, NULL, &thickness, gradient, NULL);

    gradient_magnitude = sqrt(gradient[0]*gradient[0] + gradient[1]*gradient[1]);

    // Stop if the ice thickness is small
    if (thickness < 1)
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

    status = dem->find_cell(position, i_end, j_end);

    if (status == 0) {
      is.push_front(i_end);
      js.push_front(j_end);

      if (is.size() > 5) {
        is.pop_back();
        is.pop_back();
      }
    }

  }

  int result = -5;

  // trace the streamline back:
  for (int k = 0; k < is.size(); ++k) {
    i_end = is[k];
    j_end = js[k];

    double end_mask_value = gsl_matrix_get(m, i_end, j_end);

    if (end_mask_value > 0) {
      result = end_mask_value;
      break;
    }
  }

  gsl_matrix_set(o, i_start, j_start, result);

  return 0;
}
