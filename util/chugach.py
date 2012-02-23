#!/usr/bin/env python
from netCDF4 import Dataset as NC

import numpy as np
import matplotlib.pyplot as plt

import sys

try:
    import netCDF4 as NC
except:
    import netCDF3 as NC

nc = NC.Dataset("/Users/constantine/chugach/chugachmountains.nc")

x = nc.variables['x'][:]
y = nc.variables['y'][:]
dem = nc.variables['dem'][:]

import shapefile
f = shapefile.Reader("/Users/constantine/chugach/outlines/outlines_zurich_sp.shp")

plt.figure(1)
# plt.pcolormesh(x,y,dem)

import matplotlib.nxutils as nx

def inc_mask(x, y, mask, pts):
    for j in xrange(y.size):
        for i in xrange(x.size):
            mask[j,i] += nx.pnpoly(x[i], y[j], pts)

# mask = np.zeros(dem.shape, dtype=np.int32)

for s in f.shapes():
    pts = np.array(s.points)

    if len(s.parts) == 1:
        # inc_mask(x, y, mask, pts)

        xx, yy = (pts[:,0], pts[:,1])
        plt.fill(xx, yy, color="black", alpha=0.25)
        continue

    for i in xrange(len(s.parts) - 1):
        start = s.parts[i]
        end = s.parts[i + 1]
        # inc_mask(x, y, mask, pts[start:end,:])

        xx, yy = (pts[start:end,0], pts[start:end,1])
        plt.fill(xx, yy, color="black", alpha=0.25)

plt.title(r"Chugach mountains")
plt.axis(xmin=x.min(), xmax=x.max(), ymin=y.min(), ymax=y.max())
plt.show()