#!/usr/bin/env python3

from utils import expect
import argparse, sys, pathlib
import numpy as np
from netCDF4 import Dataset
import xarray as xr

###############################################################################
def parse_command_line(args, description):
###############################################################################
    parser = argparse.ArgumentParser(
        usage="""
  {0} -i file_in -m mask_file -o file_out
    OR
  {0} --help\n\n
""".format(pathlib.Path(args[0]).name),
        description=description,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument("-i","--input_filename", required=True,
                        help="Input NC file with variables to be computed")
    parser.add_argument("-m","--mask_filename", required=True,
                        help="Input NC file with masks")
    parser.add_argument("-o","--output_filename", default="out_region_means.nc",
                        help="Output NC file with regional means on masks")
    parser.add_argument("-O","--overwrite",action='store_true',
                        help="Overwrite existing file")
    parser.add_argument("-e","--exclude", default=[], nargs='+',
                        help="Exclude these vars")

    return parser.parse_args(args[1:])

###############################################################################
def region_mean(input_filename, mask_filename, output_filename, overwrite, exclude):
###############################################################################

    f_in = pathlib.Path(input_filename).resolve().absolute()
    f_mask = pathlib.Path(mask_filename).resolve().absolute()
    f_out = pathlib.Path(output_filename).resolve().absolute()
    expect (not f_out.exists() or overwrite,
            "Output file already exists. Run with -O flag to overwrite")

    print (f"input file: {str(f_in)}")
    print (f"input mask file: {str(f_mask)}")
    print (f"output file: {str(f_out)}")

    # Open mask file and extract mask info
    mask_data = xr.open_dataset(f_mask)
    nreg_size = mask_data['mask'].max().values + 1

    # Open input file and init reg_mean(...,nreg) for all vars which have 'ncol'
    input_data = xr.open_dataset(f_in)
    dims = {}
    vdims = {}
    dtype = {}
    reg_means = {}
    vnames = []
    for vname, v in input_data.items():
        if vname=='time' or vname in exclude or ('ncol' not in v.dims):
            continue
        vnames.append(vname)
        vdims[vname] = []
        dtype[vname] = v.dtype
        for d in v.dims:
            if d == 'ncol':
                continue
            vdims[vname].append(d)
            if d not in dims:
                dims[d] = input_data.dims[d]
        vdims[vname].append('nreg')
        mean_shape = list(v.shape)
        mean_shape.remove(input_data.dims['ncol'])
        mean_shape.append(nreg_size)
        reg_means[vname] = np.zeros(mean_shape)
    dims['nreg'] = nreg_size

    # Compute regional areas
    if 'area' in vnames:
        area_name = 'area'
    elif 'area_mean' in vnames:
        area_name = 'area_mean'
    else:
        expect (False, "Input file does not contain area")
    print (f"  computing regional {area_name}")
    reg_areas = np.zeros(nreg_size)
    for ireg in range(nreg_size):
        reg_areas[ireg] = input_data[area_name].where(mask_data['mask']==ireg).sum().values

    # Create output file
    out_data = Dataset(f_out,'w', persist=True, format='NETCDF4')
    for name,size in dims.items():
        out_data.createDimension(name,size)

    # Compute regional means
    for vname in vnames:
        out_data.createVariable(vname,dtype[vname],vdims[vname])

        print (f"  computing regional {vname}")
        num = input_data[vname] * input_data[area_name]
        for ireg in range(nreg_size):
            val = num.where(mask_data['mask']==ireg).sum(dim='ncol') / reg_areas[ireg]
            if reg_means[vname].ndim == 1:
                reg_means[vname][ireg] = val
            elif reg_means[vname].ndim == 2:
                reg_means[vname][:,ireg] = val
            elif reg_means[vname].ndim == 3:
                reg_means[vname][:,:,ireg] = val
            else:
                expect (False, "Variable rank is larger than 3")
        out_data.variables[vname][:] = reg_means[vname]
    out_data.sync()
    out_data.close()
    return True

###############################################################################
def _main_func(description):
###############################################################################
    success = region_mean(**vars(parse_command_line(sys.argv, description)))

    print("OVERALL STATUS: {}".format("PASS" if success else "FAIL"))

    sys.exit(0 if success else 1)

###############################################################################
if (__name__ == "__main__"):
    _main_func(__doc__)
