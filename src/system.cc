#include "drainagebasin.hh"
#include <cmath>

int function(double t, const double y[], // inputs
             double f[],                 // output
             void* params) {             // optional input/output

  ((DEM*)params)->eval(y[0], y[1], f, NULL);

  return GSL_SUCCESS;
}

int jacobian(double t, const double y[],  // inputs
             double *dfdy, double dfdt[], // outputs
             void *params) {              // optional input/output

  ((DEM*)params)->eval(y[0], y[1], NULL, dfdy);

  dfdt[0] = 0;
  dfdt[1] = 0;

  return GSL_SUCCESS;
}
