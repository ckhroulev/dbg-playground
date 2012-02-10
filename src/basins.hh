#ifndef _BASINS_H_
#define _BASINS_H_

// GSL stuff
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>

#include "Array2D.hh"

enum MASK_VALUES {NO_VALUE = -2, ICE_FREE = -1};

int basins(double *x, int Mx, double *y, int My, double *z, double *mask);

int function(double t, const double y[], double f[], void* params);

int streamline(gsl_odeiv_system system,
               gsl_odeiv_step *step,
               int i, int j,
               int steps_per_cell,
               int path_length,
               double min_elevation,
               double max_elevation,
               Array2D<double> &mask,
               Array2D<double> &output);

#endif /* _BASINS_H_ */
