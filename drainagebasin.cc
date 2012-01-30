#include <cstdio>
#include "drainagebasin.hh"
#include <vector>
#include <cmath>

using namespace std;

// MPI stuff, NetCDF I/O
// #include <mpi.h>
// #include "PISMNC3File.hh"


int main(int argc, char **argv) {

  vector<double> X(2), Y(2), Z(4);

  X[0] = -1;
  X[1] = 1;
  Y[0] = -1;
  Y[1] = 1;

  gsl_matrix_view z_view = gsl_matrix_view_array(&Z[0], 2, 2);
  gsl_matrix * m = &z_view.matrix;

  gsl_matrix_set(m, 0, 0,  1);
  gsl_matrix_set(m, 0, 1, -1);
  gsl_matrix_set(m, 1, 1,  1);
  gsl_matrix_set(m, 1, 0, -1);

  DEM_Bilinear dem(&X[0], 2, &Y[0], 2, &Z[0]);

  gsl_odeiv2_system system = {function, jacobian, 2, &dem};

  gsl_odeiv2_step *step = gsl_odeiv2_step_alloc(gsl_odeiv2_step_bsimp, 2);

  vector<double> y(2), y_old(2), yerr(2);
  double tol = 1e-3;

  y[0] = -0.5;
  y[1] = -0.5;
  y_old = y;

  printf("set xrange [%f:%f];\n"
         "set yrange [%f:%f]\n"
         "set grid;\n"
         "plot \"-\" with linespoints title \"ODE system solution\"\n",
         X[0], X[1],
         Y[0], Y[1]);

  for (;;) {
    int status = gsl_odeiv2_step_apply(step, 0, tol, &y[0], &yerr[0], NULL, NULL, &system);

    if (status != GSL_SUCCESS)
      {
        printf ("error, return value=%d\n", status);
        break;
      }

    printf ("%.5e %.5e\n", y[0], y[1]);

    if (sqrt((y[0] - y_old[0])*(y[0] - y_old[0]) +
             (y[1] - y_old[1])*(y[1] - y_old[1])) < tol) {
      break;
    }

    y_old = y;

  }

  gsl_odeiv2_step_free (step);
  return 0;
}
