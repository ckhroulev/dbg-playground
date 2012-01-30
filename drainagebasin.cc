#include <cstdio>
#include "drainagebasin.hh"
#include <vector>
#include <cmath>

using namespace std;

// MPI stuff, NetCDF I/O
// #include <mpi.h>
// #include "PISMNC3File.hh"

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

int main(int argc, char **argv) {

  // Set up the grid (one cell):
  vector<double> X(2), Y(2), Z(4);

  X[0] = -1;
  X[1] = 1;
  Y[0] = -1;
  Y[1] = 1;

  gsl_matrix_view z_view = gsl_matrix_view_array(&Z[0], 2, 2);
  gsl_matrix * m = &z_view.matrix;

  double C = 10;

  gsl_matrix_set(m, 0, 0,  C);
  gsl_matrix_set(m, 0, 1, -C);
  gsl_matrix_set(m, 1, 1,  C);
  gsl_matrix_set(m, 1, 0, -C);

  DEM_Bilinear dem(&X[0], 2, &Y[0], 2, &Z[0]);

  gsl_odeiv2_system system = {function, jacobian, 2, &dem};

  gsl_odeiv2_step *step = gsl_odeiv2_step_alloc(
                                                // gsl_odeiv2_step_bsimp,
                                                gsl_odeiv2_step_rk2,
                                                2);

  // Set up the plot window (pipe this through gnuplot):
  printf("set xrange [%f:%f];\n"
         "set yrange [%f:%f]\n"
         "set grid;\n"
         "set multiplot\n",
         X.front(), X.back(),
         Y.front(), Y.back());

  for (int i = 0; i < 8; ++i) {
    for (int j = 0; j < 8; ++j) {
      double x = -1 + 0.25 * i,
        y = -1 + 0.25 * j;
      if (fabs(x) != fabs(y))
        flowline(system, step, x, y, "black");
    }
  }

  flowline(system, step, 0.5, 0.51, "red");

  gsl_odeiv2_step_free (step);
  return 0;
}
