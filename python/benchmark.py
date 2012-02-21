#!/usr/bin/env python

setup = """\
try:
	from netCDF3 import Dataset as NC
except:
	from netCDF4 import Dataset as NC
import sys
import numpy as np
import pylab as plt
import time
import basins

# read the DEM data
nc = NC(sys.argv[1])
x = np.array(nc.variables['x'][:], dtype=np.double)
y = np.array(nc.variables['y'][:], dtype=np.double)
thk = np.array(np.squeeze(nc.variables['thk'][:]), dtype=np.double)
z = np.array(np.squeeze(nc.variables['usurf'][:]), dtype=np.double)

# initialize the mask
tic = time.clock()
mask = basins.init_mask(thk)
toc = time.clock()
print "Mask initialization took %f seconds." % (toc - tic)
"""

import timeit
N = 10
t = timeit.Timer(setup=setup, stmt="basins.basins(x,y,z,mask,copy=True)")
times = t.repeat(repeat=5, number=N)

print map(lambda(x): x/N, times)

# plt.figure(1)
# plt.pcolormesh(x, y, mask)
# plt.contour(x, y, z, colors='black')
# plt.axis('tight')
# plt.axes().set_aspect('equal')

# plt.show()
