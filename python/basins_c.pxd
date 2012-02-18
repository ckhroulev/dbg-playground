# -*- mode: python -*-

cdef extern from "../src/basins.hh":
    bint basins(double *x, int Mx, double *y, int My, double *z, int *mask, bint output)
    bint init_mask(int Mx, int My, double *thickness, int *mask)

