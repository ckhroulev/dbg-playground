#include <mpi.h>
#include <cstdio>
#include <vector>

#include "drainagebasin.hh"
#include "DEM.hh"

using namespace std;

int main(int argc, char **argv) {

  /* MPI stuff. */
  int mpi_size, mpi_rank, mpi_namelen, ierr;
  char mpi_name[MPI_MAX_PROCESSOR_NAME];
  MPI_Comm mpi_comm = MPI_COMM_WORLD;

  /* Initialize MPI. */
  MPI_Init(&argc, &argv);
  MPI_Comm_size(mpi_comm, &mpi_size);
  MPI_Comm_rank(mpi_comm, &mpi_rank);
  MPI_Get_processor_name(mpi_name, &mpi_namelen);

  printf("# mpi_name: %s size: %d rank: %d\n", mpi_name, mpi_size, mpi_rank);

  vector<double> X, Y;          // coordinates
  Array2D<double> thk(1, 1), Z(1, 1); // the right size will be set later

  ierr = read_dem(mpi_comm, mpi_rank, "dem.nc", X, Y, Z, thk);
  if (ierr != 0) {
    printf("Initialization failed.\n");
    MPI_Finalize();
    return 1;
  }

  int Mx = X.size(), My = Y.size();
  DEM *dem = new DEM(&X[0], Mx, &Y[0], My, Z.data());

  gsl_odeiv_system system = {function, NULL, 2, dem};

  gsl_odeiv_step *step = gsl_odeiv_step_alloc(gsl_odeiv_step_rkf45, 2);

  Array2D<double> mask(Mx, My), new_mask(Mx, My);

  if (mask.allocate() != 0 || new_mask.allocate() != 0) {
    printf("Memory allocation failed.\n");
    MPI_Finalize();
    return 1;
  }

  // initialize the mask
  init_mask(thk, mask, new_mask);

  int remaining, pass_counter = 1;
  double elevation_step = 10,
    min_elevation = 0, max_elevation = elevation_step;
  do {
    remaining = 0;
    fprintf(stderr, "Pass %04d: elevation range [%4.0f, %4.0f] m...", pass_counter,
            min_elevation, max_elevation);

    for (int j = 0; j < Y.size(); j++) { // traverse in the optimal order
      for (int i = 0; i < X.size(); i++) {

        remaining += streamline(system, step, i, j,
                                2, // steps per cell
                                5, // visit this many "assigned" cells
                                min_elevation,
                                max_elevation,
                                mask, new_mask);

      }
    }

    fprintf(stderr, " done; %d cells left.\n", remaining);

    memcpy(mask.data(), new_mask.data(), Mx*My*sizeof(double));

    min_elevation = max_elevation;
    max_elevation += elevation_step;

    pass_counter++;
  } while (remaining > 0);

  ierr = write_mask(mpi_comm, mpi_rank, "mask.nc", X, Y, new_mask); CHKERRQ(ierr);

  gsl_odeiv_step_free (step);
  delete dem;

  /* Shut down MPI. */
  MPI_Finalize();

  return 0;
}
