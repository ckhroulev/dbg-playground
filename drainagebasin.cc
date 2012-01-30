
// printf, etc
#include <cstdio>

// MPI stuff, NetCDF I/O
// #include <mpi.h>
// #include "PISMNC3File.hh"

// GSL stuff
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

int function(double t, const double y[], // inputs
             double f[],                 // output
             void* params) {             // optional input/output

  f[0] = 1.0;
  f[1] = 1.0;

  return GSL_SUCCESS;
}

int jacobian(double t, const double y[],  // inputs
             double *dfdy, double dfdt[], // outputs
             void *params) {              // optional input/output
  // J(i,j) = dfdy[i * dimension + j]
  // J(i,j) = \frac{\partial f_i}{\partial y_j}

  dfdy[0 * 2 + 0] = 0;
  dfdy[0 * 2 + 1] = 0;
  dfdy[1 * 2 + 0] = 0;
  dfdy[1 * 2 + 1] = 0;

  dfdt[0] = 0;
  dfdt[1] = 0;

  return GSL_SUCCESS;
}

int main(int argc, char **argv) {
  gsl_odeiv2_system system = {function, jacobian, 2, NULL};

  gsl_odeiv2_driver * d =
    gsl_odeiv2_driver_alloc_y_new (&system,
                                   gsl_odeiv2_step_bsimp,
                                   1e-6, 1e-6, 0.0);
  int i;
  double t = 0.0, t1 = 100.0;
  double y[2] = { 0.0, 0.0 };

  printf("set grid; plot \"-\" with points title \"ODE solution\"\n");
  for (i = 1; i <= 100; i++)
    {
      double ti = i * t1 / 100.0;
      int status = gsl_odeiv2_driver_apply (d, &t, ti, y);

      if (status != GSL_SUCCESS)
	{
	  printf ("error, return value=%d\n", status);
	  break;
	}

      printf ("%.5e %.5e\n", y[0], y[1]);
    }

  gsl_odeiv2_driver_free (d);
  return 0;
}
