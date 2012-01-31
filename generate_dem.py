from PISMNC import PISMDataset as NC
import numpy as np

Mx = 101
My = 101
Lx = 20
Ly = 20

x = np.linspace(-Lx, Lx, Mx)
y = np.linspace(-Ly, Ly, My)
z = np.zeros((My, Mx))

nc = NC("dem.nc", 'w')

nc.create_dimensions(x, y)

for i in range(Mx):
    xx = x[i]
    for j in range(My):
        yy = y[j]
        z[j, i] = (xx + yy)/10.0 + np.abs(np.sin(np.pi*xx/Lx)*np.sin(np.pi*yy/Ly))

nc.write_2d_field("usurf", z)

nc.close()
