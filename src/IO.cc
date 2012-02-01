#include "drainagebasin.hh"
#include "PISMNC3File.hh"

int read_dem(MPI_Comm com, int rank,
             const char *filename,
             vector<double> &X,
             vector<double> &Y,
             vector<double> &Z,
             vector<double> &thk) {

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

  thk.resize(Mx*My);
  ierr = nc.get_vara_double("thk", start, count, &thk[0]);
  if (ierr != NC_NOERR) {
    printf("Can't read the 'usurf' variable.\n");
    ierr = nc.close();
    return 1;
  }

  ierr = nc.close();

  return 0;
}

int write_mask(MPI_Comm com, int rank,
               const char *filename,
               vector<double> &X,
               vector<double> &Y,
               double *Z) {

  int ierr;
  unsigned int Mx = X.size(), My = Y.size();
  PISMNC3File nc(com, rank);
  vector<unsigned int> start(1), count(1);

  ierr = nc.create("mask.nc");
  if (ierr != NC_NOERR) {
    printf("Can't create '%s'.\n", filename);
    return 1;
  }

  ierr = nc.def_dim("x", Mx);
  if (ierr != NC_NOERR) {
    printf("Can't create the 'x' dimension.\n");
    ierr = nc.close();
    return 1;
  }

  ierr = nc.def_dim("y", My);
  if (ierr != NC_NOERR) {
    printf("Can't create the 'y' dimension.\n");
    ierr = nc.close();
    return 1;
  }

  vector<string> dims(1);
  dims[0] = "x";
  ierr = nc.def_var("x", NC_DOUBLE, dims);
  if (ierr != NC_NOERR) {
    printf("Can't create the 'x' variable.\n");
    ierr = nc.close();
    return 1;
  }

  dims[0] = "y";
  ierr = nc.def_var("y", NC_DOUBLE, dims);
  if (ierr != NC_NOERR) {
    printf("Can't create the 'y' variable.\n");
    ierr = nc.close();
    return 1;
  }

  dims.resize(2);
  dims[0] = "x";
  dims[1] = "y";
  ierr = nc.def_var("mask", NC_DOUBLE, dims);
  if (ierr != NC_NOERR) {
    printf("Can't create the 'mask' variable.\n");
    ierr = nc.close();
    return 1;
  }

  ierr = nc.enddef();

  start[0] = 0;
  count[0] = Mx;
  ierr = nc.put_vara_double("x", start, count, &X[0]);
  if (ierr != NC_NOERR) {
    printf("Can't write the 'x' variable.\n");
    ierr = nc.close();
    return 1;
  }

  start[0] = 0;
  count[0] = My;
  ierr = nc.put_vara_double("y", start, count, &Y[0]);
  if (ierr != NC_NOERR) {
    printf("Can't write the 'y' variable.\n");
    ierr = nc.close();
    return 1;
  }

  start.resize(2);
  count.resize(2);
  start[0] = 0;
  start[1] = 0;
  count[0] = Mx;
  count[1] = My;
  ierr = nc.put_vara_double("mask", start, count, Z);
  if (ierr != NC_NOERR) {
    printf("Can't write the 'mask' variable.\n");
    ierr = nc.close();
    return 1;
  }

  ierr = nc.close();

  return 0;
}
