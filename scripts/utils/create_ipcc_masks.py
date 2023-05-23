#!/usr/bin/env python3

import argparse, sys, pathlib
import xarray as xr
import numpy as np
import regionmask as rm

###############################################################################
def parse_command_line(args, description):
###############################################################################
    parser = argparse.ArgumentParser(
        usage="""
  {0} -i file_in -o file_out
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

    return parser.parse_args(args[1:])

###############################################################################
def create_ipcc_masks(input_filename,mask_filename,overwrite):
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

    ds_in = xr.open_dataset(str(f_in))
    lat = ds_in.variables['lat']
    lon = ds_in.variables['lon']
    n = len(lat)
    mask = -1*np.ones((n),dtype=np.int32)
    for i in range(n):
        latv = np.array([lat[i]])
        lonv = np.array([lon[i]])
        mask_obj = rm.defined_regions.ar6.all.mask(lonv, latv)
        try:
            mask[i] = int(mask_obj.values[0])
        except ValueError:
            print (f"found a lat/lon pair with NaN region mask")
            print (f" - lat: {float(latv[0])}")
            print (f" - lon: {float(lonv[0])}")
            raise

    ds_out = xr.Dataset (
            {
                "lat" : lat,
                "lon" : lon,
                "mask" : mask
            }
    )
    ds_out.to_netcdf(path=f_out)
    return True

###############################################################################
def _main_func(description):
###############################################################################
    success = create_ipcc_masks(**vars(parse_command_line(sys.argv, description)))

    print("OVERALL STATUS: {}".format("PASS" if success else "FAIL"))

    sys.exit(0 if success else 1)

###############################################################################
if (__name__ == "__main__"):
    _main_func(__doc__)
