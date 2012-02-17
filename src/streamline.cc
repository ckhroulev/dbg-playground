#include "basins.hh"
#include "DEM.hh"
#include <cmath>                // sqrt
#include <map>                  // map<int,int>

using namespace std;

int function(double t, const double y[], // inputs
             double f[],                 // output
             void* params) {

  ((DEM*)params)->evaluate(y, NULL, f);

  f[0] *= -1;
  f[1] *= -1;

  return GSL_SUCCESS;
}

int streamline(gsl_odeiv_system system,
               gsl_odeiv_step *step,
               int i_start, int j_start,
               int steps_per_cell,
               int path_length,
               double min_elevation,
               double max_elevation,
               Array2D<int> &old_mask,
               Array2D<int> &new_mask) {
  DEM *dem = (DEM*)system.params;

  int mask_counter = 0,
    n_max = (dem->Mx + dem->My) * steps_per_cell,
    i = i_start, j = j_start,
    i_old, j_old,
    status;

  double step_length = dem->spacing / steps_per_cell,
    position[2],
    err[2],
    gradient[2],
    elevation,
    gradient_magnitude;

  map<int,int> values;

  int mask_value = old_mask(i_start, j_start);

  // stop if the current cell already has a value assigned
  if (mask_value > 0 || mask_value == ICE_FREE)
    return 0;

  position[0] = dem->x[i_start] + dem->dx * 0.5;
  position[1] = dem->y[j_start] + dem->dy * 0.5;

  dem->evaluate(position, &elevation, NULL);

  // if there is no ice or we're below the minimum elevation, we're done
  if (elevation < min_elevation)
    return 0;

  // if we're above the maximum elevation, wait.
  if (elevation > max_elevation)
    return 1;

  for (int step_counter = 0; step_counter < n_max; ++step_counter) {

    i_old = i; j_old = j;
    status = dem->find_cell(position, i, j);

    if (status != 0)
      break;

    mask_value = old_mask(i, j);

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

    // take a step
    status = gsl_odeiv_step_apply(step,
                                  0,         // starting time (irrelevant)
                                  step_length / gradient_magnitude, // step size (units of time)
                                  position, err, NULL, NULL, &system);

    if (status != GSL_SUCCESS) {
      printf ("error, return value=%d\n", status);
      break;
    }

  } // time-stepping loop (step_counter)

  // Find the mask value that appears more often than others.
  int most_frequent_mask_value = NO_VALUE;
  int number_of_occurences = 0;

  map<int,int>::iterator k;
  for (k = values.begin(); k != values.end(); ++k) {
    if (k->second > number_of_occurences) {
      number_of_occurences = k->second;
      most_frequent_mask_value = k->first;
    }
  }

  new_mask(i_start, j_start) = most_frequent_mask_value;

  if (most_frequent_mask_value == NO_VALUE)
    return 1;

  return 0;
}
