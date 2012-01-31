#include "drainagebasin.hh"
#include "PISMNC3File.hh"

int read_dem(MPI_Comm com, int rank,
             const char *filename,
             vector<double> &X,
             vector<double> &Y,
             vector<double> &Z) {

  int ierr;
  unsigned int Mx, My;
  PISMNC3File nc(com, rank);
  vector<unsigned int> start(1), count(1);

  ierr = nc.open("dem.nc", NC_NOWRITE);
  if (ierr != NC_NOERR) {
    printf("Can't open '%s'.\n", filename);
    return 1;
  }

  ierr = nc.inq_dimlen("x", Mx);
  if (ierr != NC_NOERR) {
    printf("Can't get the length of the 'x' dimension.\n");
    ierr = nc.close();
    return 1;
  }

  ierr = nc.inq_dimlen("y", My);
  if (ierr != NC_NOERR) {
    printf("Can't get the length of the 'y' dimension.\n");
    ierr = nc.close();
    return 1;
  }

  start[0] = 0;
  count[0] = Mx;
  X.resize(Mx);

  ierr = nc.get_vara_double("x", start, count, &X[0]);
  if (ierr != NC_NOERR) {
    printf("Can't read the 'x' variable.\n");
    ierr = nc.close();
    return 1;
  }

  start[0] = 0;
  count[0] = My;
  Y.resize(My);

  ierr = nc.get_vara_double("y", start, count, &Y[0]);
  if (ierr != NC_NOERR) {
    printf("Can't read the 'y' variable.\n");
    ierr = nc.close();
    return 1;
  }

  start.resize(2);
  count.resize(2);
  start[0] = 0;
  start[1] = 0;
  count[0] = Mx;
  count[1] = My;
  Z.resize(Mx*My);

  ierr = nc.get_vara_double("usurf", start, count, &Z[0]);
  if (ierr != NC_NOERR) {
    printf("Can't read the 'usurf' variable.\n");
    ierr = nc.close();
    return 1;
  }

  ierr = nc.close();

  return 0;
}
