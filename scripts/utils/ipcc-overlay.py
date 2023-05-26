#!/usr/bin/env python3

import argparse, sys, pathlib
import netCDF4 as nc
import numpy as np
import regionmask as rm
import cartopy.crs as ccrs
import matplotlib.animation as animation
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt

###############################################################################
def parse_command_line(args, description):
###############################################################################
    parser = argparse.ArgumentParser(
        usage="""
  {0} -i file_in
    OR
  {0} --help\n\n
""".format(pathlib.Path(args[0]).name),
        description=description,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument("-i","--input_filename", required=True,
                        help="NC file with fields to plot")
    parser.add_argument("-r","--radians", action='store_true',
                        help="whether lat/lon in input file are in radians")
    parser.add_argument("-l","--lev", default=None, type=int,
                        help="If variable has a lev/ilev, slice it a this vertical level")
    parser.add_argument("-o","--output_filename", default="movie.mp4",
                        help="Name of mp4 file to generate")
    parser.add_argument("-m","--mask-file", default=None,
                        help="Location of the nc file with masks for the regions where var is defined")
    parser.add_argument("-n","--mask-name", default=None,
                        help="Name of the mask variable in mask file")
    parser.add_argument("varname",
                        help="Variable to plot")

    return parser.parse_args(args[1:])

###############################################################################
def plot_over_ipcc_regions(input_filename,radians,lev,output_filename,mask_file,mask_name,varname):
###############################################################################

    dpi = 100

    # Create figure, and set some static info/data
    fig = plt.figure()
    #  fig, ax = plt.subplots()
    ax = plt.axes(projection=ccrs.EckertIII())
    ax.coastlines()
    ax.set_global()

    #  Overlay IPCC regions
    ar6 = rm.defined_regions.ar6.all
    text_kws = dict(
        bbox=dict(color="none"),
        path_effects=[pe.withStroke(linewidth=1, foreground="w",alpha=0.5)],
        fontsize=6,
    )
    ar6.plot_regions(ax=ax, line_kws=dict(lw=0.5), text_kws=text_kws)

    # Load lat/lon, and adjust depending on units and ranges
    ds = nc.Dataset(input_filename,mode='r')
    lat = ds.variables['lat'][:]
    lon = ds.variables['lon'][:]
    if radians:
        print (" -> converting lat/lon to degrees")
        lat = lat*180/np.pi
        lon = lon*180/np.pi

    min_lat = np.amin(lat);
    min_lon = np.amin(lon);
    if min_lon >= 0:
        print (f" -> min_lon: {min_lon}. Subtracting 180, to get it in [-180,180]")
        lon = lon - 180
    if min_lat >= 0:
        print (f" -> min_lat: {min_lat}. Subtracting 90, to get it in [-90,90]")
        lat = lat - 90

    # Get time dimension length
    has_time = 'time' in ds.variables.keys()
    if not has_time:
        print ("Cannot handle datasets without time")
        raise
    time = ds.variables['time']
    nt = len(time[:])

    # Get variable layout
    var = ds.variables[varname]
    dims = var.dimensions
    ndims = len(dims)
    if ndims>3:
        raise ValueError(f"Unsupported variable layout: {dims}.")
    time_dim = dims.index('time')
    if time_dim!=0:
        raise ValueError(f"We expect time to be the slowest dim, but layout is {dims}")
    has_ncol = 'ncol' in dims
    if has_ncol:
        ncol_dim = dims.index('ncol')
        ncols = ds.dimensions['ncol'].size
    else:
        if mask_file is None or mask_name is None:
            raise ValueError(f"var does not have ncol dim. It must come from a masked integral, so -m/--mask-file and -n/--mask-name are needed\n")
        dsm = nc.Dataset(mask_file,mode='r')
        if mask_name not in dsm.variables.keys():
            raise ValueError(f"Could not locate mask var '{mask_name}' in mask file '{mask_file}'")
        mask = dsm.variables[mask_name]
        if len(mask.dimensions)!=1 or mask.dimensions[0]!='ncol':
            raise ValueError(f"Unexpected layout for mask variable: {mask.dimensions}")
        # Since we don't know how mask labels are distributed, we need to map mask labels to [0,N),
        # with N=num_mask_labels. We *ASSUME* that cldera-tools creates a mask stat where, say,
        # the 10-th entry corresponds to the 10-th smallest mask value
        ncols = dsm.dimensions['ncol'].size
        m2idx = {}
        for m in mask[:]:
            if m not in m2idx.keys():
                idx = len(m2idx)
                m2idx[m] = idx

    has_lev = False
    if 'lev' in dims or 'ilev' in dims:
        if lev is None:
            raise ValueError(f"Variable has a vertical level dimension: {var.dimensions}. "
                               "You must provide -l/--lev N argument.")
        lev_dim = dims.index('lev') if 'lev' in dims else dims.index('ilev')
        has_lev = True

    # Helper function, to get slice of data at given time index
    def get_data_at_time(n):
        if has_lev:
            data_t = np.take(var,n,time_dim)
            data = np.take(data_t,lev,lev_dim-1)
        else:
            data = np.take(var,n,time_dim)
        return data

    if has_ncol:
        filled_c = ax.tricontourf(lon, lat, get_data_at_time(0),
                       transform=ccrs.PlateCarree(),alpha=1.0)
        fig.colorbar(filled_c, orientation="horizontal")
        def plot_at_time (n):
            data = get_data_at_time(n)
            filled_c = ax.tricontourf(lon, lat, data,
                           transform=ccrs.PlateCarree(),alpha=0.8)

        ani = animation.FuncAnimation(fig,plot_at_time,nt,interval=100,repeat=False)
        writer = animation.writers['ffmpeg'](fps=10)
        ani.save(output_filename,writer=writer,dpi=dpi)
    else:
        # We need to create a "heat map"
        def create_data_2d(data):
            data2d = np.zeros((ncols))
            for i in range(0,ncols):
                mval = int(mask[i])
                idx = m2idx[mval]
                data2d[i] = data[idx]
            return data2d
        filled_c = ax.tricontourf(lon, lat, create_data_2d(get_data_at_time(0)),
                       transform=ccrs.PlateCarree(),alpha=1.0)
        fig.colorbar(filled_c, orientation="horizontal")
        def plot_at_time (n):
            data = get_data_at_time(n)
            data2d = create_data_2d(data)
            filled_c = ax.tricontourf(lon, lat, data2d,
                           transform=ccrs.PlateCarree(),alpha=0.8)

        ani = animation.FuncAnimation(fig,plot_at_time,nt,interval=100,repeat=False)
        writer = animation.writers['ffmpeg'](fps=10)
        ani.save(output_filename,writer=writer,dpi=dpi)

    #  plt.show()
    return True

###############################################################################
def _main_func(description):
###############################################################################
    success = plot_over_ipcc_regions(**vars(parse_command_line(sys.argv, description)))

    print("OVERALL STATUS: {}".format("PASS" if success else "FAIL"))

    sys.exit(0 if success else 1)

###############################################################################
if (__name__ == "__main__"):
    _main_func(__doc__)

