#include "drainagebasin.hh"
#include <vector>
#include <cmath>

using namespace std;

int flowline(gsl_odeiv2_system system,
             gsl_odeiv2_step *step,
             double x0, double y0, const char *color) {
  vector<double> y(2),          // new position
    y_old(2),                   // old position
    yerr(2),                    // error estimate
    dy(2);

  DEM *dem = (DEM*)system.params;

  double step_size = (dem->dx() < dem->dy() ? dem->dx() : dem->dy()) / 2.0;

  int counter = 0;

  y[0] = x0;
  y[1] = y0;
  y_old = y;

  printf("plot \"-\" with lines linecolor rgb \"%s\"\n",
         color);

  printf ("%.5e %.5e\n", y[0], y[1]);

  for (;;) {
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
    if (sqrt(dy[0]*dy[0] + dy[1]*dy[1]) < 8e3) {
      // y[0] += step_size;        // FIXME
      counter++;

      if (counter > 1)
        break;
    }

    if (y[0] <= dem->x_min() || y[0] >= dem->x_max() ||
        y[1] <= dem->y_min() || y[1] >= dem->y_max())
      break;

    y_old = y;
  }

  printf ("e\n");

  return 0;
}
