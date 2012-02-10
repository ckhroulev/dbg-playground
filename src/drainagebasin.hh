// MPI and C++ stuff
#include <mpi.h>
#include <vector>
using namespace std;

#include "Array2D.hh"

#define CHKERRQ(e) do { \
    if ((e) != 0) {                                              \
      printf("Bailing out in file %s, line %d.\n", __FILE__, __LINE__); \
      return e;                                                         \
    }                                                                   \
  } while (0)

int read_dem(MPI_Comm com, int rank,
             const char *filename,
             vector<double> &X,
             vector<double> &Y,
             Array2D<double> &Z,
             Array2D<double> &thk);

int write_mask(MPI_Comm com, int rank,
               const char *filename,
               vector<double> &X,
               vector<double> &Y,
               Array2D<double> &Z);

void init_mask(Array2D<double> &thickness,
               Array2D<double> &mask);
