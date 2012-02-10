# -*- mode: python -*-

import numpy as np
cimport numpy as np
cimport cfoo

ctypedef np.double_t double_t

def foo_py(A):
    cdef np.ndarray[double_t, ndim=2, mode="c"] A_c
    A_c = np.ascontiguousarray(A, dtype="d")
    cfoo.foo(<double*>A_c.data, A_c.size)

    return A_c
