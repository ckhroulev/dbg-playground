# -*- mode: python -*-
import numpy as np
cimport numpy as np

cimport basins

ctypedef np.double_t double_t

def basins(np.ndarray[dtype=double_t, ndim=1] x,
           np.ndarray[dtype=double_t, ndim=1] y,
           np.ndarray[dtype=double_t, ndim=2, mode="c"] z,
           np.ndarray[dtype=double_t, ndim=2, mode="c"] mask, bool inplace):

    cdef np.ndarray[dtype=double_t, ndim=2, mode="c"] output

    if x.shape != mask.shape:
        raise ValueError("arguments z and mask have to have the same shape")

    if y.size != z.shape[0]:
        raise ValueError("the size of y has to match the number of rows in z")

    if x.size != z.shape[1]:
        raise ValueError("the size of x has to match the number of columns in z")

    if inplace:
        output = mask
    else:
        output = mask.copy()

    basins.basins(<double*>x.data, x.size,
                  <double*>y.data, y.size
                  <double*>z.data,
                  <double*>output.data)

    return output
