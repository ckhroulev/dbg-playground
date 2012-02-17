# -*- mode: python -*-

cdef extern from "../src/basins.hh":
    bint basins(double *x, int Mx, double *y, int My, double *z, int *mask, bint output)
