#!/usr/bin/env python

setup = """\
try:
	from netCDF3 import Dataset as NC
except:
	from netCDF4 import Dataset as NC
import sys
import numpy as np
import dbg

# read the DEM data
nc = NC(sys.argv[1])
x = np.array(nc.variables['x'][:], dtype=np.double)
y = np.array(nc.variables['y'][:], dtype=np.double)
thk = np.array(np.squeeze(nc.variables['thk'][:]), dtype=np.double)
z = np.array(np.squeeze(nc.variables['usurf'][:]), dtype=np.double)

# initialize the mask
mask = dbg.initialize_mask(thk)
print "Mask initialization: done"
"""

import timeit
N = 10
t = timeit.Timer(setup=setup, stmt="dbg.upslope_area(x,y,z,mask,copy=True)")
times = t.repeat(repeat=5, number=N)

print map(lambda(x): x/N, times)
