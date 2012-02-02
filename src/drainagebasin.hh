
// GSL stuff
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_math.h>       // GSL_NAN

#include <mpi.h>
#include <vector>

using namespace std;

class DEM {
public:
  DEM(double *x, int Mx, double *y, int My, double *z, double *thickness);
  ~DEM();

  void evaluate(const double *position, double *elevation, double *thickness,
                double *f, double *jac);

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

  double x_min();
  double x_max();

  double y_min();
  double y_max();

  double dx();
  double dy();

  double get_x(int i);
  double get_y(int j);

  int get_Mx();
  int get_My();

protected:
  double *x, *y, *z, *thk;
  int Mx, My;
  double x_spacing, y_spacing, one_over_dx, one_over_dy;

  inline void get_corner_values(int i, int j, double *data,
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
    gsl_matrix_view data_view = gsl_matrix_view_array(data, Mx, My);
    gsl_matrix * m = &data_view.matrix;

    A = gsl_matrix_get(m, i,     j);
    B = gsl_matrix_get(m, i,     j + 1);
    C = gsl_matrix_get(m, i + 1, j + 1);
    D = gsl_matrix_get(m, i + 1, j);
  }
};

int function(double t, const double y[], double f[], void* params);
int jacobian(double t, const double y[], double *dfdy, double dfdt[], void *params);

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
               double min_elevation,
               double max_elevation,
               double *mask,
               double *output);
