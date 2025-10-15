#!/usr/bin/env python3

import argparse, sys, pathlib
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

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
def create_sai_masks(input_filename,mask_filename,overwrite,degrees):
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

    # define a mask id function
    def get_mask_id (latv,lonv):
        tol_lat = 1.5 # any column within X degrees lat of target site (rectangle)
        tol_lon = 1.5 # any column within X degrees lon of target site (rectangle)
        target_lats = [-50.0, -30.0, -15.0, 15.0, 30.0, 50.0] # we want 6 potential sites
        for ilat,target_lat in enumerate(target_lats):
            if abs(latv - target_lat) < tol_lat and abs(lonv-180.0) < tol_lon:
                return ilat+1
        return 0 # non injection sites
        
    for i in range(n):
        latv = np.array([lat[i]])
        lonv = np.array([lon[i]])
        if not degrees:
            latv = latv*rad2deg
            lonv = lonv*rad2deg
        mask[i] = get_mask_id(latv,lonv)
        if mask[i] < 0:
            raise SystemExit("Error: zone not found in lat = " + str(latv))

    # double-check how many grid cells are assigned to each zone (sanity check)
    counts = np.zeros((7))
    for j in range(7):
        for i in range(n):
            if mask[i] == j:
                counts[j] = counts[j] + 1
    print("Counts for each mask location: ",counts)

    ds_out = Dataset(f_out,'w',persist=True, format='NETCDF4')
    ncol_name = lat.dimensions[0]
    ncol = ds_in.dimensions[ncol_name].size
    ds_out.createDimension(ncol_name,ncol)
    lat_out = ds_out.createVariable("lat","f8",lat.dimensions)
    lon_out = ds_out.createVariable("lon","f8",lon.dimensions)
    mask_out = ds_out.createVariable("mask","i4",lon.dimensions)

    lat_out[:] = lat[:]
    lat_out.units = 'degrees_north'
    lon_out[:] = lon[:]
    lon_out.units = 'degrees_west'
    mask_out[:] = mask[:]
    mask_out.units = "0 for no injection, positive for injection site"
    mask_out.description = "Native grid masking of injection sites"

    # plot the mask we created
    fig = plt.figure(figsize=(10, 6))
    ax = plt.axes(projection=ccrs.EckertIII())
    ax.coastlines()
    ax.set_global()
    ax.gridlines(xlocs=range(-180,180,60), ylocs=range(-90,90,30),draw_labels=True)

    contour = ax.tricontourf(lon[:], lat[:], mask, s=5, transform=ccrs.PlateCarree())
    cbar = plt.colorbar(contour, ax=ax, orientation='vertical', pad=0.05)
    cbar.set_label('Mask id')

    plt.title('SAI Injection Site Map')
    plt.show()

    ds_out.sync()
    ds_out.close()
    ds_in.close()

    return True

###############################################################################
def _main_func(description):
###############################################################################
    success = create_sai_masks(**vars(parse_command_line(sys.argv, description)))

    print("OVERALL STATUS: {}".format("PASS" if success else "FAIL"))

    sys.exit(0 if success else 1)

###############################################################################
if (__name__ == "__main__"):
    _main_func(__doc__)
