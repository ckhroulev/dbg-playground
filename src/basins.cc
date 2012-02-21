#include "basins.hh"
#include "DEM.hh"

int basins(double *x, int Mx, double *y, int My, double *z, int *mask, bool output) {
  int remaining = 0;
  double elevation_step = 10,
    min_elevation = 0, max_elevation = elevation_step;

  DEM dem(x, Mx, y, My, z);

  Array2D<int> my_mask(Mx, My), new_mask(Mx, My);
  my_mask.wrap(mask);
  if (new_mask.allocate() != 0)
    return -1;

#pragma omp parallel default(shared)
  {
    gsl_odeiv_system system = {function, NULL, 2, &dem};
    gsl_odeiv_step *step = gsl_odeiv_step_alloc(gsl_odeiv_step_rkf45, 2);

#pragma omp for schedule(dynamic)
    for (int j = 0; j < My; j++)
      for (int i = 0; i < Mx; i++)
        new_mask(i, j) = my_mask(i, j);

    do {
#pragma omp for schedule(dynamic)
      for (int j = 0; j < My; j++)
        for (int i = 0; i < Mx; i++)
          my_mask(i, j) = new_mask(i, j);

#pragma omp single
      remaining = 0;

#pragma omp for schedule(dynamic) reduction(+:remaining)
      for (int j = 0; j < My; j++) { // Note: traverse in the optimal order
        for (int i = 0; i < Mx; i++) {
          remaining += streamline(system, step, i, j,
                                  2, // steps per cell
                                  5, // visit this many "assigned" cells
                                  min_elevation,
                                  max_elevation,
                                  my_mask, new_mask);
        }
      }

#pragma omp single
      {
        min_elevation = max_elevation;
        max_elevation += elevation_step;
      }

    } while (remaining > 0);

    gsl_odeiv_step_free (step);

#pragma omp for schedule(dynamic)
    for (int j = 0; j < My; j++)
      for (int i = 0; i < Mx; i++)
	my_mask(i, j) = new_mask(i, j);

  } // end of the parallel block

  return 0;
}
