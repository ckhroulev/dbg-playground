#include <mpi.h>

#include <cstdio>
#include <vector>
#include <cmath>

#include "drainagebasin.hh"

using namespace std;

#define CHKERRQ(e) do { \
    if ((e) != 0) {                                              \
      printf("Bailing out in file %s, line %d.\n", __FILE__, __LINE__); \
      return e;                                                         \
    }                                                                   \
  } while (0)

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

  vector<double> X, Y,          // coordinates
    Z, thk;                     // elevation, interpreted as a 2D array

  ierr = read_dem(mpi_comm, mpi_rank, "dem.nc", X, Y, Z, thk);
  if (ierr != 0) {
    printf("Initialization failed.\n");
    MPI_Finalize();
    return 1;
  }

  DEM *dem = new DEM(&X[0], X.size(), &Y[0], Y.size(), &Z[0], &thk[0]);

  gsl_odeiv2_system system = {function, jacobian, 2, dem};

  gsl_odeiv2_step *step = gsl_odeiv2_step_alloc(
                                                gsl_odeiv2_step_rk2,
                                                2);

  double *mask = new double[X.size() * Y.size()];
  double *new_mask = new double[X.size() * Y.size()];

  if (mask == NULL) {
    printf("Memory allocation failed.\n");
    MPI_Finalize();
    return 1;
  }

  // initialize the mask
  init_mask(X.size(), Y.size(), &Z[0], &thk[0], mask);

  for (int k = 0; k < 1; ++k) {
    fprintf(stderr, "Pass %d...", k+1);
    for (int i = 0; i < X.size(); i++) {
      for (int j = 0; j < Y.size(); j++) {
        streamline(system, step, i, j, mask, new_mask);
      }
    }
    fprintf(stderr, " done.\n");
  }

  ierr = write_mask(mpi_comm, mpi_rank, "mask.nc", X, Y, new_mask); CHKERRQ(ierr);

  gsl_odeiv2_step_free (step);
  delete[] mask;
  delete[] new_mask;
  delete dem;

  /* Shut down MPI. */
  MPI_Finalize();

  return 0;
}
