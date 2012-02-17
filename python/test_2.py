#!/usr/bin/env python

from netCDF4 import Dataset as NC
import sys
import numpy as np
import pylab as plt
import time

# read the DEM data
nc = NC(sys.argv[1])
x = np.array(nc.variables['x'][:], dtype=np.double)
y = np.array(nc.variables['y'][:], dtype=np.double)
thk = np.squeeze(nc.variables['thk'][:])
z = np.array(np.squeeze(nc.variables['usurf'][:]), dtype=np.double)

# initialize the mask
mask = np.zeros_like(thk, dtype=np.double) - 1           # ice free

ii, jj = np.meshgrid([-1, 0, 1], [-1, 0, 1])

counter = 1
for i in range(1, x.size - 1):
    for j in range(1, y.size - 1):
        if thk[j, i] > 1:
            if np.any(thk[jj + j, ii + i] < 1):
                mask[j, i] = counter
                counter += 1
            else:
                mask[j, i] = -2

import basins

tic = time.clock()
basins.basins(x, y, z, mask, True, False)
toc = time.clock()

print "This computation took %f seconds." % (toc - tic)

plt.figure(2)
plt.pcolormesh(x, y, mask)
plt.contour(x, y, z, colors='black')
plt.axis('tight')
plt.axes().set_aspect('equal')

plt.show()
