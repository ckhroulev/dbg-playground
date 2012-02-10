# -*- mode: python -*-

import numpy as np
cimport numpy as np
cimport cfoo

ctypedef np.double_t double_t

def foo_py(np.ndarray A):
    if A.dtype != np.float64:
        raise ValueError("only numpy.float64 arrays are supported")

    if A.flags.c_contiguous == False:
        raise ValueError("only C-contiguous arrays are supported")

    if A.ndim != 2:
        raise ValueError("only 2D arrays are supported")

    cfoo.foo(<double*>A.data, A.size)

    return A
