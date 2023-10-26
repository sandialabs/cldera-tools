#!/usr/bin/env python3

import argparse, sys, pathlib
import numpy as np
from netCDF4 import Dataset

###############################################################################
def parse_command_line(args, description):
###############################################################################
    parser = argparse.ArgumentParser(
        usage="""
  {0} -i file_in -m file_out
    OR
  {0} --help\n\n
""".format(pathlib.Path(args[0]).name),
        description=description,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument("-i","--input_filename", required=True,
                        help="NC file with lat/lon coords")
    parser.add_argument("-m","--mask_filename", required=True,
                        help="Name of output NC file with masks")
    parser.add_argument("-o","--overwrite",action='store_true',
                        help="Overwrite existing file")
    parser.add_argument("-d","--degrees", action='store_true',
                        help="whether lat/lon in input file are in degrees")

    return parser.parse_args(args[1:])

###############################################################################
def create_zonal_masks(input_filename,mask_filename,overwrite,degrees):
###############################################################################

    f_in = pathlib.Path(input_filename).resolve().absolute()
    f_out = pathlib.Path(mask_filename).resolve().absolute()

    if f_out.exists():
        if not overwrite:
            print (f"Cannot create output mask file. File exists. Run with -o to overwrite")
            return False
        if f_out.is_dir():
            print (f"Cannot create output mask file. File is an existing directory.")
            return False
    
    print (f"input lat/lon file: {str(f_in)}")
    print (f"output mask file: {str(f_out)}")

    ds_in = Dataset(f_in,'r')
    lat = ds_in.variables['lat']
    lon = ds_in.variables['lon']
    n = len(lat)
    mask = -1*np.ones((n),dtype=np.int32)
    rad2deg = 180.0 / np.pi
    for i in range(n):
        def get_mask_id (latv,lonv):
            if latv < -66.5:
                return 0 # Polar S
            elif latv < -35.0: # 
                return 1 # Temperate S
            elif latv < -23.5:
                return 2 # Subtropical S
            elif latv < 23.5:
                return 3 # Tropical
            elif latv < 35.0:
                return 4 # Subtropical N
            elif latv < 66.5:
                return 5 # Temperate N
            elif latv < 90.0:
                return 6 # Polar N
            else:
                return -1
        latv = np.array([lat[i]])
        lonv = np.array([lon[i]])
        if not degrees:
            latv = latv*rad2deg
            lonv = lonv*rad2deg
        mask[i] = get_mask_id(latv,lonv)
        if mask[i] < 0:
            raise SystemExit("Error: zone not found in lat = " + str(latv))


    ds_out = Dataset(f_out,'w',persist=True, format='NETCDF4')
    ncol_name = lat.dimensions[0]
    ncol = ds_in.dimensions[ncol_name].size
    ds_out.createDimension(ncol_name,ncol)
    lat_out = ds_out.createVariable("lat","f8",lat.dimensions)
    lon_out = ds_out.createVariable("lon","f8",lon.dimensions)
    mask_out = ds_out.createVariable("mask","i4",lon.dimensions)

    lat_out[:] = lat[:]
    lon_out[:] = lon[:]
    mask_out[:] = mask[:]

    ds_out.sync()
    ds_out.close()
    ds_in.close()
    return True

###############################################################################
def _main_func(description):
###############################################################################
    success = create_zonal_masks(**vars(parse_command_line(sys.argv, description)))

    print("OVERALL STATUS: {}".format("PASS" if success else "FAIL"))

    sys.exit(0 if success else 1)

###############################################################################
if (__name__ == "__main__"):
    _main_func(__doc__)
