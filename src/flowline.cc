#include "drainagebasin.hh"
#include <vector>
#include <cmath>

using namespace std;

int flowline_mask(gsl_odeiv2_system system,
                  gsl_odeiv2_step *step,
                  double x0, double y0,
                  int marker, bool &new_terminus,
                  double *mask) {
  double y[2],          // new position
    y_old[2],                   // old position
    yerr[2],                    // error estimate
    dy[2];

  DEM *dem = (DEM*)system.params;

  double grid_spacing = dem->dx() < dem->dy() ? dem->dx() : dem->dy(),
    steps_per_cell = 10.0, n_max = (dem->Mx + dem->My) * steps_per_cell;
  int counter;

  dem->path.resize(n_max);

  y[0] = x0;
  y[1] = y0;
  y_old[0] = x0;
  y_old[1] = y0;

  double f[2], f_mag, step_size;
  dem->eval(y[0], y[1], f, NULL);

  f_mag = sqrt(f[0]*f[0] + f[1]*f[1]);

  if (f_mag < 1e-16)
    return 0;

  for (counter = 0; counter < n_max; ++counter) {

    dem->indices(y[0], y[1], dem->path[counter].i, dem->path[counter].j);

    if (y[0] <= dem->x_min() || y[0] >= dem->x_max() ||
        y[1] <= dem->y_min() || y[1] >= dem->y_max())
      break;

    dem->eval(y[0], y[1], f, NULL);

    f_mag = sqrt(f[0]*f[0] + f[1]*f[1]);

    // Stop if the magnitude of the gradient is very low or very high.
    if (f_mag < 1e-16 || f_mag > 10)
      break;

    step_size = (grid_spacing / steps_per_cell) / f_mag;

    // take a step
    int status = gsl_odeiv2_step_apply(step,
                                       0,         // starting time (irrelevant)
                                       step_size, // step size
                                       &y[0], &yerr[0], NULL, NULL, &system);

    if (status != GSL_SUCCESS) {
      printf ("error, return value=%d\n", status);
      break;
    }

    y_old[0] = y[0];
    y_old[1] = y[1];
  }

  dem->path.resize(counter);

  gsl_matrix_view mask_view = gsl_matrix_view_array(mask, dem->Mx, dem->My);
  gsl_matrix * m = &mask_view.matrix;
  new_terminus = true;

  for (int k = counter - 1; k >= 0; --k) {
    double current_value = gsl_matrix_get(m,
                                          dem->path[k].i,
                                          dem->path[k].j);

    if (current_value <= marker && k == counter - 1) {
      marker = current_value;
      new_terminus = false;
    }

    if (gsl_isnan(current_value) || current_value > marker) {
      gsl_matrix_set(m,
                     dem->path[k].i,
                     dem->path[k].j,
                     marker);
    }
  }

  return 0;
}


int flowline_gnuplot(gsl_odeiv2_system system,
                     gsl_odeiv2_step *step,
                     double x0, double y0, const char *color) {
  vector<double> y(2),          // new position
    y_old(2),                   // old position
    yerr(2),                    // error estimate
    dy(2);

  DEM *dem = (DEM*)system.params;

  double grid_spacing = dem->dx() < dem->dy() ? dem->dx() : dem->dy(),
    steps_per_cell = 10.0, n_max = (dem->Mx + dem->My) * steps_per_cell;

  y[0] = x0;
  y[1] = y0;
  y_old = y;

  double f[2], f_mag, step_size;
  dem->eval(y[0], y[1], f, NULL);

  f_mag = sqrt(f[0]*f[0] + f[1]*f[1]);

  if (f_mag < 1e-16)
    return 0;

  // start the plot:
  printf("plot \"-\" with lines linecolor rgb \"%s\" notitle\n",
         color);

  // plot the initial position
  printf ("%.5e %.5e\n", y[0], y[1]);

  for (int counter = 0; counter < n_max;++counter) {

    if (y[0] <= dem->x_min() || y[0] >= dem->x_max() ||
        y[1] <= dem->y_min() || y[1] >= dem->y_max())
      break;

    dem->eval(y[0], y[1], f, NULL);

    f_mag = sqrt(f[0]*f[0] + f[1]*f[1]);

    // Stop if the magnitude of the gradient is very low or very high.
    if (f_mag < 1e-16 || f_mag > 10)
      break;

    step_size = (grid_spacing / steps_per_cell) / f_mag;

    // make a step
    int status = gsl_odeiv2_step_apply(step,
                                       0,         // starting time (irrelevant)
                                       step_size, // step size
                                       &y[0], &yerr[0], NULL, NULL, &system);

    if (status != GSL_SUCCESS) {
      printf ("error, return value=%d\n", status);
      break;
    }

    printf ("%.5e %.5e\n", y[0], y[1]);
    dy[0] = y[0] - y_old[0];
    dy[1] = y[1] - y_old[1];

    if (sqrt(dy[0]*dy[0] + dy[1]*dy[1]) < grid_spacing / (2*steps_per_cell))
      break;

    y_old = y;
  }

  printf ("e\n");

  return 0;
}
