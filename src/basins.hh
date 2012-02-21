#ifndef _BASINS_H_
#define _BASINS_H_

// GSL stuff
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>

#include "Array2D.hh"

enum MASK_VALUES {NO_VALUE = -2, ICE_FREE = -1};

int init_mask(int Mx, int My, double *thickness, int* mask);

int basins(double *x, int Mx, double *y, int My, double *z, int *mask, bool output);

int function(double t, const double y[], double f[], void* params);

struct coloring_context {
  gsl_odeiv_system system;
  gsl_odeiv_step *step;
  int steps_per_cell, path_length;
  double min_elevation, max_elevation;
};

int streamline(coloring_context ctx,
               int i, int j,
               Array2D<int> &mask,
               Array2D<int> &output);

#endif /* _BASINS_H_ */
