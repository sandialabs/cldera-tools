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
    parser.add_argument("-c","--cmp", default=None,
                        help="If variable has a 3d dim, slice it a this entry")
    parser.add_argument("-o","--output_filename", default="movie.mp4",
                        help="Name of mp4 file to generate")
    parser.add_argument("varname",
                        help="Variable to plot")

    return parser.parse_args(args[1:])

###############################################################################
def plot_over_ipcc_regions(input_filename,radians,cmp,output_filename,varname):
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

    has_time = 'time' in ds.variables.keys()
    if not has_time:
        print ("Cannot handle datasets without time")
        raise

    # Loop over time dimension for requested var
    var = ds.variables[varname]
    ndims = len(var.dimensions)
    if ndims>3:
        raise ValueError(f"Unsupported variable layout: {var.dimensions}.")

    dims = var.dimensions
    time_dim = dims.index('time')
    ncol_dim = dims.index('ncol')
    if ndims==3:
        if cmp is None:
            raise ValueError(f"Variable has 3+ dimensions: {var.dimensions}. "
                               "You must provide -c/--cmp N argument.")
        comp_dims = [dims.index(i) for i in dims if i!='time' and i!='ncol']
        if len(comp_dims)!=1:
            raise ValueError(f"Something went wrong while locating cmp dim in {dims}")
        comp_dim = comp_dims[0]
        cmp = int(cmp)

    time = ds.variables['time']
    nt = len(time[:])

    def get_data_at_time(n):
        if ndims==3:
            data_t = np.take(var,n,time_dim)
            data = np.take(data_t,cmp,comp_dim-1)
        else:
            data = np.take(var,n,time_dim)
        return data

    filled_c = ax.tricontourf(lon, lat, get_data_at_time(0),
                   transform=ccrs.PlateCarree(),alpha=1.0)
    fig.colorbar(filled_c, orientation="horizontal")
    def plot_at_time (n):
        data = get_data_at_time(n)
        print (f"shape data: {np.shape(data)}")
        filled_c = ax.tricontourf(lon, lat, data,
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

