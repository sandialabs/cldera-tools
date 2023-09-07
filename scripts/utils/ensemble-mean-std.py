#!/usr/bin/env python3

from utils import expect
import argparse, sys, pathlib
import numpy as np
from netCDF4 import Dataset

###############################################################################
def parse_command_line(args, description):
###############################################################################
    parser = argparse.ArgumentParser(
        usage="""
  {0} -f ens-filename-root -n N -o output-filename
    OR
  {0} --help\n\n
""".format(pathlib.Path(args[0]).name),
        description=description,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument("-f","--filename", required=True,
                        help="Filename for ensemble NC files. Special char '#' will expand into ens member number")
    parser.add_argument("-n","--num-files", required=True,type=int,
                        help="Number of ensemble files")
    parser.add_argument("-s","--start", default=1,type=int,
                        help="Start index for ensemble members")
    parser.add_argument("-w","--width", default=0,type=int,
                        help="Width of number in ens member string (default: width of start+num_files")
    parser.add_argument("-o","--output-file",default="out.nc",
                        help="Name of output nc file")
    parser.add_argument("-O","--overwrite",action='store_true',
                        help="Overwrite existing file")
    parser.add_argument("-e","--exclude", default=[], nargs='+',
                        help="Exclude these vars")

    return parser.parse_args(args[1:])

###############################################################################
def compute_ens_mean_std(filename,num_files,start,width,output_file,overwrite,exclude):
###############################################################################

    expect (num_files>1, f"Input number of files must greater than 1 (got {num_files})")
    expect ("#" in filename, "Input filename should contain the char '#', which will expand to the number of the ens member to parse")

    # ens files have int str attached
    nwidth = width if width>0 else len(str(start+num_files))

    # Open the first ens file, to init list of ens mean/std
    one_str = str(start).zfill(nwidth)
    ens1_filename = filename.replace("#",one_str)
    ens1_file = pathlib.Path(ens1_filename).resolve().absolute()
    ens1_db = Dataset(ens1_file,'r')

    dims = {}
    vdims = {}
    dtype = {}
    mean = {}
    std  = {}
    vnames = []
    for vname in ens1_db.variables.keys():
        if vname=='time' or vname in exclude:
            continue
        vnames.append(vname)
        v = ens1_db.variables[vname]
        vdims[vname] = []
        dtype[vname] = v.dtype
        for d in v.dimensions:
            vdims[vname].append(d)
            if d not in dims:
                dims[d] = ens1_db.dimensions[d].size
        mean[vname] = np.zeros(v.shape)
        std[vname] = np.zeros(v.shape)
    ens1_db.close()

    # Compute ens mean
    for n in range(start,start+num_files):
        num = str(n).zfill(nwidth)
        print (f"updating means with ens member {num}")
        ensN_filename = filename.replace('#',num)
        ensN_file = pathlib.Path(ensN_filename).resolve().absolute()
        ensN_db = Dataset(ensN_file,'r')

        for vname in vnames:
            print (f"  updating {vname}")
            expect (vname in ensN_db.variables.keys(),
                    f"Ensemble file {ensN_file} is missing variable {vname}")

            vN = ensN_db.variables[vname]
            mean[vname][:] += vN[:]

        ensN_db.close()

    print (f"dividing means by {num_files}")
    for vmean in mean.values():
        vmean[:] /= num_files

    # Compute ens sample standard deviation
    for n in range(start,start+num_files):
        num = str(n).zfill(nwidth)
        print (f"updating std devs with ens member {num}")
        ensN_filename = filename.replace('#',num)
        ensN_file = pathlib.Path(ensN_filename).resolve().absolute()
        ensN_db = Dataset(ensN_file,'r')

        for vname in vnames:
            print (f"  updating {vname}")
            expect (vname in ensN_db.variables.keys(),
                    f"Ensemble file {ensN_file} is missing variable {vname}")

            vN = ensN_db.variables[vname]
            std[vname][:] += np.square(vN[:]-mean[vname])

        ensN_db.close()

    print (f"dividing stds by {num_files-1} and taking square root")
    for vstd in std.values():
        vstd[:] / (num_files-1)
        vstd[:] = np.sqrt(vstd)

    # Create output file
    f_out = pathlib.Path(output_file).resolve().absolute()
    expect (not f_out.exists() or overwrite,
            "Output file already exists. Run with -O flag to overwrite")

    ds_out = Dataset(f_out,'w',persist=True, format='NETCDF4')
    for name,size in dims.items():
        ds_out.createDimension(name,size)

    for vname in vnames:
        ds_out.createVariable(vname+"_mean",dtype[vname],vdims[vname])
        ds_out.createVariable(vname+"_std",dtype[vname],vdims[vname])

        vmean = ds_out.variables[vname+"_mean"]
        vstd  = ds_out.variables[vname+"_std"]

        vmean[:] = mean[vname]
        vstd[:]  = std[vname]

    ds_out.sync()
    ds_out.close()
    return True

###############################################################################
def _main_func(description):
###############################################################################
    success = compute_ens_mean_std(**vars(parse_command_line(sys.argv, description)))

    print("OVERALL STATUS: {}".format("PASS" if success else "FAIL"))

    sys.exit(0 if success else 1)

###############################################################################
if (__name__ == "__main__"):
    _main_func(__doc__)
