# -*- mode: python -*-

cdef extern from "foobar.h":
    bint foobar(double* data, int length, double* output)

