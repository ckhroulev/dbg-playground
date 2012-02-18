# -*- mode: python -*-
cimport numpy as np
import numpy as np

np.import_array()

cimport foobar_c

ctypedef np.double_t double_t

def foobar(np.ndarray[dtype=double_t, mode="c"] data):
    cdef np.ndarray[dtype=double_t, mode="c"] output
    cdef int length

    length = data.size

    output = np.zeros_like(data)

    foobar_c.foobar(<double*>data.data, length,
                    <double*>output.data)

    return output
