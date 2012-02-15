# -*- mode: python -*-
import numpy as np
cimport numpy as np

cimport basins

ctypedef np.double_t double_t

def basins(np.ndarray[dtype=double_t, ndim=1] x,
           np.ndarray[dtype=double_t, ndim=1] y,
           np.ndarray[dtype=double_t, ndim=2, mode="c"] z,
           np.ndarray[dtype=double_t, ndim=2, mode="c"] mask, int inplace):
    """
    arguments:
    - x, y: 1D arrays with coordinates
    - z: surface elevation, a 2D array
    - mask: mask, integers (FIXME)
    - inplace: boolean; True if the mask is to be modified in place
    """
    cdef np.ndarray[dtype=double_t, ndim=2, mode="c"] output

    if not (z.shape[0] == mask.shape[0] and z.shape[1] == mask.shape[1]):
        raise ValueError("arguments z and mask have to have the same shape: got (%d,%d) and (%d,%d)" %
                         (z.shape[0], z.shape[1], mask.shape[0], mask.shape[1]))

    if y.size != z.shape[0]:
        raise ValueError("the size of y has to match the number of rows in z")

    if x.size != z.shape[1]:
        raise ValueError("the size of x has to match the number of columns in z")

    if inplace:
        output = mask
    else:
        output = mask.copy()

    basins.basins(<double*>x.data, x.size,
                  <double*>y.data, y.size,
                  <double*>z.data,
                  <double*>output.data)

    return output
