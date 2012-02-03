#!/usr/bin/env python
from mpl_toolkits.basemap import Basemap
from netCDF4 import Dataset as NC

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import sys

import shapefile
f = shapefile.Reader("Greenland_basins_geog/Greenland_basins.shp")

print "Trying to open '%s'" % sys.argv[1]
nc = NC(sys.argv[1], 'r')
nc_mask = NC(sys.argv[2], 'r')

# we need to know longitudes and latitudes corresponding to our grid
lon = nc.variables['lon'][:]
lat = nc.variables['lat'][:]

# x and y *in the dataset* are only used to determine plotting domain
# dimensions
x = nc.variables['x'][:]
y = nc.variables['y'][:]
width = x.max() - x.min()
height = y.max() - y.min()

# load data and mask out ice-free areas
usurf = np.squeeze(nc.variables['usurf'][:])

mask  = np.squeeze(nc_mask.variables['mask'][:])
mask  = np.ma.array(mask, mask=(mask < 0))

m = Basemap(width=width,      # width in projection coordinates, in meters
            height=height,      # height
            resolution='l',     # coastline resolution, can be 'l' (low), 'h'
                                # (high) and 'f' (full)
            projection='stere', # stereographic projection
            lat_ts=71,          # latitude of true scale
            lon_0=-41,          # longitude of the plotting domain center
            lat_0=72)           # latitude of the plotting domain center

# draw the Blue Marble background (requires PIL, the Python Imaging Library)
#m.bluemarble()

# convert longitudes and latitudes to x and y:
xx,yy = m(lon, lat)

m.pcolormesh(xx, yy, mask)
m.contour(xx, yy, usurf, colors="black")

# draw parallels and meridians. The labels argument specifies where to draw
# ticks: [left, right, top, bottom]
# m.drawparallels(np.arange(-55.,90.,5.), labels = [1, 0, 0, 0])
# m.drawmeridians(np.arange(-120.,30.,10.), labels = [0, 0, 0, 1])

# plot basins from a shapefile
for s in f.shapes():
    pts = np.array(s.points)
    xx, yy = m(pts[:,0], pts[:,1])
    m.plot(xx, yy, color="white")

plt.title(r"Greenland surface elevation and drainage basins")
plt.show()
#plt.savefig('greenland_drainage_basins.png')
