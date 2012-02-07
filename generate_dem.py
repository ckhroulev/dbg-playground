from PISMNC import PISMDataset as NC
import numpy as np

Mx = 201
My = 201
Lx = 20e4
Ly = 20e4

x   = np.linspace(-Lx, Lx, Mx)
y   = np.linspace(-Ly, Ly, My)
z   = np.zeros((My, Mx))

nc = NC("dem.nc", 'w')

nc.create_dimensions(x, y)

for i in range(Mx):
    xx = x[i]
    for j in range(My):
        yy = y[j]
        r = np.sqrt(xx*xx + yy*yy)
        if r < 15e4:
            z[j, i] = np.sqrt(15e4*15e4 - xx*xx - yy*yy)/150.0

nc.write_2d_field("usurf", z)
nc.write_2d_field("thk",   z)

nc.close()
