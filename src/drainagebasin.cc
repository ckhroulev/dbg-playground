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
    Z;                          // elevation, interpreted as a 2D array

  ierr = read_dem(mpi_comm, mpi_rank, "dem.nc", X, Y, Z);
  if (ierr != 0) {
    printf("Initialization failed.\n");
    MPI_Finalize();
    return 1;
  }

  DEM_Bilinear *dem = new DEM_Bilinear(&X[0], X.size(), &Y[0], Y.size(), &Z[0]);

  gsl_odeiv2_system system = {function, jacobian, 2, dem};

  gsl_odeiv2_step *step = gsl_odeiv2_step_alloc(
                                                gsl_odeiv2_step_rk2,
                                                2);

  // Set up the plot window (pipe this through gnuplot):
  // printf("set xrange [%f:%f];\n"
  //        "set yrange [%f:%f]\n"
  //        "set size ratio %f;\n"
  //        "set format xy \"\"\n"
  //        "set multiplot\n",
  //        X.front(), X.back(),
  //        Y.front(), Y.back(),
  //        (Y.back() - Y.front()) / (X.back() - X.front()));

  double *mask = new double[X.size() * Y.size()];
  if (mask == NULL) {
    printf("Memory allocation failed.\n");
    MPI_Finalize();
    return 1;
  }

  for (int n = 0; n < X.size() * Y.size(); ++n)
    mask[n] = GSL_NAN;

  int k = 1;
  bool new_terminus = false;
  for (int i = 0; i < X.size(); i++) {
    for (int j = 0; j < Y.size(); j++) {
      flowline_mask(system, step, X[i], Y[j], k, new_terminus, mask);
      if (new_terminus)
        k++;
    }
  }

  ierr = write_mask(mpi_comm, mpi_rank, "mask.nc", X, Y, mask); CHKERRQ(ierr);

  gsl_odeiv2_step_free (step);
  delete[] mask;
  delete dem;

  /* Shut down MPI. */
  MPI_Finalize();

  return 0;
}
