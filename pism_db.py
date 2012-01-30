#!/usr/bin/env python

import numpy as np
from netCDF4 import Dataset as NC
import scipy.integrate as si

nc = NC(sys.argv[1], 'r')

usurf = nc.variables['usurf'][:]
x = nc.variables['x'][:]
y = nc.variables['y'][:]

