# -*- mode: python -*-
cimport numpy as np

cimport basins_c

ctypedef np.double_t double_t
ctypedef np.int_t int_t

def basins(np.ndarray[dtype=double_t, ndim=1] x,
           np.ndarray[dtype=double_t, ndim=1] y,
           np.ndarray[dtype=double_t, ndim=2, mode="c"] z,
           np.ndarray[dtype=int_t, ndim=2, mode="c"] mask,
           copy = False, print_output = False):
    """
    arguments:
    - x, y: 1D arrays with coordinates
    - z: surface elevation, a 2D array
    - mask: mask, integers, a 2D array
    - inplace: boolean; True if the mask is to be modified in place
    """
    cdef np.ndarray[dtype=int_t, ndim=2, mode="c"] output

    # z and mask are typed, so z.shape is not a Python object.
    # This means that we have to compare z.shape[0,1] to mask.shape[0,1] 'by hand'.
    if not (z.shape[0] == mask.shape[0] and z.shape[1] == mask.shape[1]):
        raise ValueError("arguments z and mask have to have the same shape: got (%d,%d) and (%d,%d)" %
                         (z.shape[0], z.shape[1], mask.shape[0], mask.shape[1]))

    if y.size != z.shape[0]:
        raise ValueError("the size of y has to match the number of rows in z")

    if x.size != z.shape[1]:
        raise ValueError("the size of x has to match the number of columns in z")

    if copy:
        output = mask.copy()
    else:
        output = mask

    basins_c.basins(<double*>x.data, x.size,
                    <double*>y.data, y.size,
                    <double*>z.data,
                    <int*>output.data,
                    print_output)

    return output
