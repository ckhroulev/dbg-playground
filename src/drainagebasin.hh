// GSL stuff
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>

#include <mpi.h>
#include <vector>
#include <cmath>

using namespace std;

struct Node {
  double elevation, mask, new_mask;
};

enum MASK_VALUES {NO_VALUE = -2, ICE_FREE = -1};

class DEM {
public:
  DEM(double *x, int Mx, double *y, int My, Node *dem);
  ~DEM();

  void evaluate(const double *position, double *elevation, double *f);

  inline int find_cell(const double *position, int &i, int &j) {
    i = floor((position[0] - x[0]) * one_over_dx);
    j = floor((position[1] - y[0]) * one_over_dy);

    // bail if we ended up outside the grid
    if (i < 0 || i + 1 > Mx - 1 || j < 0 || j + 1 > My - 1) {
      i = j = -1;
      return 1;
    }

    return 0;
  }

  inline double get_x(int i)
  { return x[i]; }

  inline double get_y(int j)
  { return y[j]; }

  inline double x_min()
  { return x[0]; }

  inline double x_max()
  { return x[Mx-1]; }

  inline double y_min()
  { return y[0]; }

  inline double y_max()
  { return y[My-1]; }

  inline double dx()
  { return x_spacing; }

  inline double dy()
  { return y_spacing; }

  inline int get_Mx()
  { return Mx; }

  inline int get_My()
  { return My; }

  inline double get_spacing()
  { return spacing; }

protected:
  double *x, *y;
  Node *dem;
  int Mx, My;
  double spacing, x_spacing, y_spacing, one_over_dx, one_over_dy;

  inline void get_corner_values(int i, int j, Node *data,
                                double &A, double &B, double &C, double &D) {
    // Get the surface elevation at grid corners (arranged like so):
    //
    //   ^ y
    //   |
    //   |
    //   B-----C
    //   |     |
    //   | *   |   x
    // --A-----D---->
    //   |
#define DATA(i,j) ((data)[(j)*Mx + (i)].elevation)
    A = DATA(i,     j);
    B = DATA(i,     j + 1);
    C = DATA(i + 1, j + 1);
    D = DATA(i + 1, j);
#undef DATA(i,j)
  }
};

int function(double t, const double y[], double f[], void* params);

int read_dem(MPI_Comm com, int rank,
             const char *filename,
             vector<double> &X,
             vector<double> &Y,
             vector<double> &Z,
             vector<double> &thk);

int write_mask(MPI_Comm com, int rank,
               const char *filename,
               vector<double> &X,
               vector<double> &Y,
               double *Z);

void init_mask(int Mx, int My,
               double *thickness,
               double *mask,
               double *tmp);

int streamline(gsl_odeiv_system system,
               gsl_odeiv_step *step,
               int i, int j,
               int steps_per_cell,
               int path_length,
               double min_elevation,
               double max_elevation,
               double *mask,
               double *output);
