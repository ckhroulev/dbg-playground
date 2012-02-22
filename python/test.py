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

plt.figure(1)
plt.pcolormesh(x, y, mask)
plt.contour(x, y, z, colors='black')
plt.axis('tight')
plt.axes().set_aspect('equal')

pts = get_terminus()

import matplotlib.nxutils as nx

def correct_mask(pts):
        for j in range(y.size):
            for i in range(x.size):
                if mask[j,i] > 0:
                    if nx.pnpoly(x[i], y[j], pts):
                        mask[j,i] = 2
                    else:
                        mask[j,i] = 1

correct_mask(pts)

dbg.upslope_area(x, y, z, mask)
print "Drainage basin computation: done"

plt.figure(1)
plt.pcolormesh(x, y, mask)
plt.contour(x, y, z, colors='black')
plt.axis('tight')
plt.axes().set_aspect('equal')
plt.title("upslope areas")

# mask2 = np.zeros(thk.shape, dtype=np.float64) - 1
# mask2[thk > 5] = 0
# mask2[mask != 2] = -1

# dbg.accumulated_flow(x, y, z, mask2, n_samples = 4)

# mask2[mask2 == -1] = 0
# mask2 += 0.01

# plt.figure(2)
# plt.pcolormesh(x, y, np.log10(mask2))
# plt.colorbar()
# plt.contour(x, y, z, colors='black')
# plt.axis('tight')
# plt.axes().set_aspect('equal')
# plt.title("accumulated flow")

plt.show()
