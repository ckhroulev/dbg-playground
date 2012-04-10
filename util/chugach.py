#!/usr/bin/env python
from netCDF4 import Dataset as NC

import numpy as np
import matplotlib.pyplot as plt

import sys

try:
    import netCDF4 as NC
except:
    import netCDF3 as NC

from rasterize import *

nc = NC.Dataset("chugachmountains.nc")

x = nc.variables['x'][:]
y = nc.variables['y'][:]
dem = nc.variables['dem'][:]
dem[dem <= 0] = 0

x_min, x_max = x[0], x[-1]
y_min, y_max = y[0], y[-1]

img = rasterize("outlines/outlines_zurich_sp.shp", dem.shape[1], dem.shape[0], "white", "black",
                x_range = [x_min, x_max], y_range = [y_min, y_max])

img.save("chugach.png")

mask = np.asarray(img)

from PISMNC import PISMDataset as OUT

nc = OUT("mask.nc", 'w')

nc.create_dimensions(x, y)

nc.write_2d_field("mask", mask)
nc.write_2d_field("usurf", dem)

nc.close()

# plt.figure(1)
# plt.pcolormesh(x,y,dem)


def plot_shapes():
    import shapefile
    f = shapefile.Reader("outlines/outlines_zurich_sp.shp")

    plt.figure(2)
    for s in f.shapes():
        pts = np.array(s.points)

        if len(s.parts) == 1:
            xx, yy = (pts[:,0], pts[:,1])

            # dbg.increment_mask(x, y, xx, yy, mask)
            # sys.stdout.write(".")

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

#plot_shapes()
