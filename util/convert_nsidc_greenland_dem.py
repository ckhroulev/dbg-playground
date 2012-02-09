from PISMNC import PISMDataset as NC
from numpy import fromfile, int32, r_

shape = (2782, 2611)
x = r_[0:shape[1]:1]
y = r_[0:shape[0]:1]

#download data from http://nsidc.org/data/docs/daac/nsidc0304_0305_glas_dems.gd.html

dem = fromfile("../data/NSIDC_Grn1km_wgs84_elev_cm.dat", dtype=">i", count=-1)
dem = dem.reshape(shape) / 100.0        # convert to meters

nc = NC("NSIDC_Greenland.nc", 'w')
nc.create_dimensions(x, y)
nc.write_2d_field("usurf", dem)
nc.write_2d_field("thk",   dem)
nc.close()

err_distance = fromfile("../data/NSIDC_Grn1km_dist_mm.dat", dtype=">i", count=-1)
err_distance = err_distance.reshape(shape) / 100.0        # convert to meters

nc = NC("NSIDC_Greenland_err_distance.nc", 'w')
nc.create_dimensions(x, y)
nc.write_2d_field("err_distance", err_distance)
nc.close()

    # char polar_stereographic ;
    #     polar_stereographic:Northernmost_Northing = -628500. ;
    #     polar_stereographic:Southernmost_Northing = -3410500. ;
    #     polar_stereographic:Easternmost_Easting = 1720500. ;
    #     polar_stereographic:Westernmost_Easting = -890500. ;
