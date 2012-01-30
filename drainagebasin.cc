#include <cstdio>

// MPI stuff, NetCDF I/O
// #include <mpi.h>
// #include "PISMNC3File.hh"


int main(int argc, char **argv) {
  gsl_odeiv2_system system = {function, jacobian, 2, NULL};

  gsl_odeiv2_driver * d =
    gsl_odeiv2_driver_alloc_y_new (&system,
                                   gsl_odeiv2_step_bsimp,
                                   1e-6, 1e-6, 0.0);
  int i;
  double t = 0.0, t1 = 100.0;
  double y[2] = { 0.0, 0.0 };

  printf("set grid; plot \"-\" with points title \"ODE system solution\"\n");
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
