#include <mpi.h>
#include <cstdio>
#include <vector>

#include "drainagebasin.hh"
#include "DEM.hh"
#include "basins.hh"

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
  Array2D<double> mask(Mx, My);

  if (mask.allocate() != 0) {
    printf("Initialization failed.\n");
    MPI_Finalize();
    return 1;
  }

  init_mask(thk, mask);

  ierr = basins(&X[0], Mx, &Y[0], My, Z.data(), mask.data()); CHKERRQ(ierr);

  ierr = write_mask(mpi_comm, mpi_rank, "mask.nc", X, Y, mask); CHKERRQ(ierr);

  /* Shut down MPI. */
  MPI_Finalize();

  return 0;
}
