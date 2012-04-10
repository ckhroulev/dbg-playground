#!/usr/bin/env python

import dbg

from PISMNC import PISMDataset as NC

nc = NC("mask.nc")

mask = nc.variables['mask'][:]
usurf = nc.variables['usurf'][:]
x = nc.variables['x'][:]
y = nc.variables['y'][:]

thk = mask.copy()
thk[mask == 0] = 100
thk[mask == 255] = 0

mask = dbg.initialize_mask(thk)

nc = NC("flow.nc", 'w')

nc.create_dimensions(x, y)

nc.write_2d_field("mask", mask)
nc.write_2d_field("thk", thk)

nc.close()
