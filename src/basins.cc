#include <cstring>

#include "basins.hh"
#include "DEM.hh"

int basins(double *x, int Mx, double *y, int My, double *z, double *mask) {
  int remaining, pass_counter = 1;
  double elevation_step = 10,
    min_elevation = 0, max_elevation = elevation_step;

  DEM dem(x, Mx, y, My, z);

  Array2D<double> my_mask(Mx, My), new_mask(Mx, My);
  my_mask.wrap(mask);
  if (new_mask.allocate() != 0)
    return -1;

  memcpy(new_mask.data(), my_mask.data(), Mx*My*sizeof(double));

  gsl_odeiv_system system = {function, NULL, 2, &dem};
  gsl_odeiv_step *step = gsl_odeiv_step_alloc(gsl_odeiv_step_rkf45, 2);

  do {
    remaining = 0;
    printf("Pass %04d: elevation range [%4.0f, %4.0f] m...", pass_counter,
           min_elevation, max_elevation);
    fflush(stdout);

    for (int j = 0; j < My; j++) { // traverse in the optimal order
      for (int i = 0; i < Mx; i++) {

        remaining += streamline(system, step, i, j,
                                2, // steps per cell
                                5, // visit this many "assigned" cells
                                min_elevation,
                                max_elevation,
                                my_mask, new_mask);

      }
    }

    printf(" done; %d cells left.\n", remaining);
    fflush(stdout);

    memcpy(my_mask.data(), new_mask.data(), Mx*My*sizeof(double));

    min_elevation = max_elevation;
    max_elevation += elevation_step;

    pass_counter++;
  } while (remaining > 0);

  gsl_odeiv_step_free (step);

  return 0;
}
