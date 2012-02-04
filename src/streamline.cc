#include "drainagebasin.hh"
#include <cmath>
#include <map>

using namespace std;

int function(double t, const double y[], // inputs
             double f[],                 // output
             void* params) {

  ((DEM*)params)->evaluate(y, NULL, f);

  return GSL_SUCCESS;
}

#define OLD_MASK(i,j) old_mask[(j) * dem->get_Mx() + (i)]
#define NEW_MASK(i,j) new_mask[(j) * dem->get_Mx() + (i)]

int streamline(gsl_odeiv_system system,
               gsl_odeiv_step *step,
               int i_start, int j_start,
               int steps_per_cell,
               int path_length,
               double min_elevation,
               double max_elevation,
               double *old_mask, double *new_mask) {
  DEM *dem = (DEM*)system.params;

  int counter, mask_counter = 0,
    n_max = (dem->get_Mx() + dem->get_My()) * steps_per_cell,
    i = 0, j = 0, i_old, j_old, status;

  double grid_spacing = dem->get_spacing(),
    position[2], err[2];

  double gradient[2], elevation, gradient_magnitude, step_size;

  map<int,int> values;

  int mask_value = OLD_MASK(i_start, j_start);

  // stop if the current cell already has a value assigned
  if (mask_value > 0 || mask_value == ICE_FREE)
    return 0;

  position[0] = dem->get_x(i_start) + dem->dx() * 0.5;
  position[1] = dem->get_y(j_start) + dem->dy() * 0.5;

  dem->evaluate(position, &elevation, NULL);

  // if there is no ice or we're below the minimum elevation, we're done
  if (elevation < min_elevation)
    return 0;

  // if we're above the maximum elevation, wait.
  if (elevation > max_elevation)
    return 1;

  for (counter = 0; counter < n_max; ++counter) {

    i_old = i; j_old = j;
    status = dem->find_cell(position, i, j);

    if (status != 0)
      break;

    mask_value = OLD_MASK(i, j);

    if (mask_value == ICE_FREE)   // ice-free
      break;

    if ((i != i_old || j != j_old) && (mask_value > 0)) {
      values[mask_value]++;
      mask_counter++;

      if (mask_counter == path_length)
        break;
    }

    dem->evaluate(position, NULL, gradient);

    gradient_magnitude = sqrt(gradient[0]*gradient[0] + gradient[1]*gradient[1]);

    step_size = (grid_spacing / steps_per_cell) / gradient_magnitude;

    // take a step
    status = gsl_odeiv_step_apply(step,
                                  0,         // starting time (irrelevant)
                                  step_size, // step size
                                  position, err, NULL, NULL, &system);

    if (status != GSL_SUCCESS) {
      printf ("error, return value=%d\n", status);
      break;
    }

  } // time-stepping loop (counter)

  // Find the mask value that appears more often than others.
  int result = NO_VALUE;
  int n = 0;

  map<int,int>::iterator k;
  for (k = values.begin(); k != values.end(); ++k) {
    if (k->second > n) {
      n = k->second;
      result = k->first;
    }
  }

  NEW_MASK(i_start, j_start) = result;

  if (result == 0)
    return 1;

  return 0;
}

#undef OLD_MASK
#undef NEW_MASK
