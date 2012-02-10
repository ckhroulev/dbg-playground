# -*- mode: python -*-

cdef extern from "cfoo.h":
    bint foo(double *data, int len)
