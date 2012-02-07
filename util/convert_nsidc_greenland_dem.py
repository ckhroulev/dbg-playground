from PISMNC import PISMDataset as NC
from numpy import fromfile, int32, r_

shape = (2782, 2611)
x = r_[0:shape[1]:1]
y = r_[0:shape[0]:1]

dem = fromfile("NSIDC_Grn1km_wgs84_elev_cm.dat", dtype=">i", count=-1)
dem = dem.reshape(shape) / 100.0        # convert to meters

nc = NC("NSIDC_Greenland.nc", 'w')
nc.create_dimensions(x, y)
nc.write_2d_field("usurf", dem)
nc.write_2d_field("thk",   dem)
nc.close()

    # char polar_stereographic ;
    #     polar_stereographic:Northernmost_Northing = -628500. ;
    #     polar_stereographic:Southernmost_Northing = -3410500. ;
    #     polar_stereographic:Easternmost_Easting = 1720500. ;
    #     polar_stereographic:Westernmost_Easting = -890500. ;
