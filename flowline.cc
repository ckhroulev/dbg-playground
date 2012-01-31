#include "drainagebasin.hh"
#include <vector>
#include <cmath>

using namespace std;

int flowline(gsl_odeiv2_system system,
             gsl_odeiv2_step *step,
             double x0, double y0, const char *color) {
  vector<double> y(2),          // new position
    y_old(2),                   // old position
    yerr(2);                    // error estimate

  double step_size = 1e-2;

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

    if (sqrt((y[0] - y_old[0])*(y[0] - y_old[0]) +
             (y[1] - y_old[1])*(y[1] - y_old[1])) < step_size*step_size) {
      // y[0] += step_size;        // FIXME
      counter++;

      if (counter > 5)
        break;
    }

    y_old = y;

  }

  printf ("e\n");

  return 0;
}
