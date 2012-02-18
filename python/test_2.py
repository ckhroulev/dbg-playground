#!/usr/bin/env python

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

tic = time.clock()
basins.basins(x, y, z, mask)
toc = time.clock()
print "Drainage basin computation took %f seconds." % (toc - tic)

plt.figure(2)
plt.pcolormesh(x, y, mask)
plt.contour(x, y, z, colors='black')
plt.axis('tight')
plt.axes().set_aspect('equal')

plt.show()
