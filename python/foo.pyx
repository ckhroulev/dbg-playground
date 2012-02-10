# -*- mode: python -*-
import numpy as np
cimport numpy as np

cimport cfoo

ctypedef np.double_t double_t

def foo(np.ndarray[dtype=double_t, ndim=2, mode="c"] A):
    cdef np.ndarray[dtype=double_t, ndim=2, mode="c"] output

    output = A.copy()
    cfoo.foo(<double*>output.data, output.size)

    return output
