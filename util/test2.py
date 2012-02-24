#!/usr/bin/env python

try:
    from netCDF3 import Dataset as NC
except:
    from netCDF4 import Dataset as NC

import sys
import numpy as np
import pylab as plt
import dbg

def get_terminus():
    from matplotlib.widgets import Cursor
    def tellme(s):
        print s
        plt.title(s,fontsize=16)
        plt.draw()

    plt.setp(plt.gca(),autoscale_on=False)

    cursor = Cursor(plt.axes(), useblit=True, color='white', linewidth=1 )

    happy = False
    while not happy:
        pts = []
        while len(pts) < 4:
            tellme('Select 4 corners of the terminus region')
            pts = np.asarray( plt.ginput(4, timeout=-1) )
            if len(pts) < 4:
                tellme('Too few points, starting over')
                time.sleep(1) # Wait a second

        ph = plt.fill(pts[:,0], pts[:,1], 'white', lw = 2, alpha=0.5)

        tellme('Done? Press any key if yes, mouse click to reset')

        happy = plt.waitforbuttonpress()

        # Get rid of fill
        if not happy:
            for p in ph: p.remove()

        return pts

# read the DEM data
nc = NC(sys.argv[1])
x = np.array(nc.variables['x'][:], dtype=np.double)
y = np.array(nc.variables['y'][:], dtype=np.double)
thk = np.array(np.squeeze(nc.variables['thk'][:]), dtype=np.double)
z = np.array(np.squeeze(nc.variables['usurf'][:]), dtype=np.double)

# initialize the mask
mask = dbg.initialize_mask(thk)
print "Mask initialization: done"

new_mask = dbg.upslope_area(x, y, z, mask, copy=True)
print "Drainage basin computation: done"

plt.figure(1)
plt.pcolormesh(x, y, new_mask)
plt.contour(x, y, z, colors='black')
plt.axis('tight')
plt.axes().set_aspect('equal')

pts = get_terminus()

import matplotlib.nxutils as nx

def get_colors(mask, pts):
    result = []
    for j in range(y.size):
        for i in range(x.size):
            if mask[j,i] > 0 and nx.pnpoly(x[i], y[j], pts):
                result.append(mask[j,i])
    return result

colors = get_colors(new_mask, pts)

for j in xrange(y.size):
    for i in xrange(x.size):
        if new_mask[j,i] in colors:
            new_mask[j,i] = 1
        else:
            if new_mask[j,i] > 0:
                new_mask[j,i] = 0
            else:
                new_mask[j,i] = -1

plt.figure(2)
plt.pcolormesh(x, y, new_mask)
plt.contour(x, y, z, colors='black')
plt.axis('tight')
plt.axes().set_aspect('equal')
plt.fill(pts[:,0], pts[:,1], 'white', lw = 2, alpha=0.5)

plt.show()
