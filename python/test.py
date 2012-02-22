#!/usr/bin/env python

try:
	from netCDF3 import Dataset as NC
except:
	from netCDF4 import Dataset as NC

import sys
import numpy as np
import pylab as plt
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

dbg.upslope_area(x, y, z, mask)
print "Drainage basin computation: done"

mask2 = np.zeros(thk.shape, dtype=np.float64) - 1
mask2[thk > 5] = 0

dbg.accumulated_flow(x, y, z, mask2, n_samples = 4)

mask2[mask2 == -1] = 0
mask2 += 0.01

plt.figure(1)
plt.pcolormesh(x, y, mask)
plt.contour(x, y, z, colors='black')
plt.axis('tight')
plt.axes().set_aspect('equal')
plt.title("upslope areas")

plt.figure(2)
plt.pcolormesh(x, y, np.log10(mask2))
plt.colorbar()
plt.contour(x, y, z, colors='black')
plt.axis('tight')
plt.axes().set_aspect('equal')
plt.title("accumulated flow")

plt.show()
