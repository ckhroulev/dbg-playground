#include <cmath>

class DEM {
public:
  DEM(double *x, int Mx, double *y, int My, double *z);
  ~DEM() {}

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

  double *x, *y, *z;
  double spacing, dx, dy;
  int Mx, My;
protected:
  double one_over_dx, one_over_dy;
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
#define DATA(i,j) ((data)[(j)*Mx + (i)])
    A = DATA(i,     j);
    B = DATA(i,     j + 1);
    C = DATA(i + 1, j + 1);
    D = DATA(i + 1, j);
#undef DATA(i,j)
  }
};
