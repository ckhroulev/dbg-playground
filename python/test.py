#!/usr/bin/env python
import numpy as np
import pylab as plt

Mx = 401
My = 201
Lx = 20e4
Ly = 20e4

x   = np.linspace(-Lx, Lx, Mx)
y   = np.linspace(-Ly, Ly, My)
thk   = np.zeros((My, Mx))

for i in range(Mx):
    xx = x[i]
    for j in range(My):
        yy = y[j]
        r = np.sqrt(xx*xx + yy*yy)
        if r < 15e4:
            thk[j, i] = np.sqrt(15e4*15e4 - xx*xx - yy*yy)/150.0

xx,yy = np.meshgrid(x, y)
z = 0.025 * xx + thk

# initialize the mask
mask = np.zeros_like(thk) - 2           # mark everything as "no value"
mask[thk <= 1] = -1                     # mark ice free areas as such

ii = np.zeros_like(mask, dtype=np.bool)
np.logical_and(thk < 200, thk > 0, ii)  # indices of the margin
np.logical_and(ii, z < -2500, ii)

mask[ii] = 10

plt.figure(1)
plt.contour(x, y, z, 40)
plt.colorbar()

plt.figure(2)
plt.pcolormesh(x, y, mask)
plt.colorbar()

import basins

basins.basins(x, y, z, mask, True, True)

plt.figure(3)
plt.pcolormesh(x, y, mask)
plt.colorbar()

plt.show()
