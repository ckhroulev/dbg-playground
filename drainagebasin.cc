#include <mpi.h>

#include <cstdio>
#include <vector>
#include <cmath>

#include "drainagebasin.hh"
#include "PISMNC3File.hh"

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

  printf("mpi_name: %s size: %d rank: %d\n", mpi_name, mpi_size, mpi_rank);

  // Set up the grid (read from a file)
  vector<double> X, Y, Z;
  unsigned int Mx, My;
  PISMNC3File nc(mpi_comm, mpi_rank);

  ierr = nc.open("dem.nc", NC_NOWRITE); CHKERRQ(ierr);

  ierr = nc.inq_dimlen("x", Mx); CHKERRQ(ierr);
  ierr = nc.inq_dimlen("y", My); CHKERRQ(ierr);

  vector<unsigned int> start(1), count(1);
  start[0] = 0;
  count[0] = Mx;
  X.resize(Mx);
  ierr = nc.get_vara_double("x", start, count, &X[0]); CHKERRQ(ierr);

  count[0] = My;
  Y.resize(My);
  ierr = nc.get_vara_double("y", start, count, &Y[0]); CHKERRQ(ierr);

  start.resize(2);
  count.resize(2);
  start[0] = 0;
  start[1] = 0;
  count[0] = Mx;
  count[1] = My;
  Z.resize(Mx*My);
  ierr = nc.get_vara_double("usurf", start, count, &Z[0]); CHKERRQ(ierr);

  ierr = nc.close(); CHKERRQ(ierr);

  return 0;                     // FIXME

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

  // for (int i = 0; i < 8; ++i) {
  //   for (int j = 0; j < 8; ++j) {
  //     double x = -1 + 0.25 * i,
  //       y = -1 + 0.25 * j;
  //     if (fabs(x) != fabs(y))
  //       flowline(system, step, x, y, "black");
  //   }
  // }

  flowline(system, step, 0.5, 0.51, "red");

  gsl_odeiv2_step_free (step);

  /* Shut down MPI. */
  MPI_Finalize();

  return 0;
}
