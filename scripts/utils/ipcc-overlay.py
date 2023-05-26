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

    parser.add_argument("varname",
                        help="Variable to plot")
    parser.add_argument("-i","--input_filename", required=True,
                        help="NC file with fields to plot")
    parser.add_argument("-r","--radians", action='store_true',
                        help="whether lat/lon in input file are in radians")
    parser.add_argument("-lonp","--lon-positive", action='store_true',
                        help="whether lon is positive. If not, assume centered")
    parser.add_argument("-latp","--lat-positive", action='store_true',
                        help="whether lat is positive. If not, assume centered")

    return parser.parse_args(args[1:])

###############################################################################
def plot_over_ipcc_regions(input_filename,radians,lon_positive,lat_positive,varname):
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
        print ("convert to degrees")
        lat = lat*180/np.pi
        lon = lon*180/np.pi
    if lon_positive:
        print ("subtract 180 deg to lon")
        lon = lon - 180
    if lat_positive:
        print ("subtract 90 deg to lat")
        lat = lat - 90

    has_time = 'time' in ds.variables.keys()
    if not has_time:
        print ("Cannot handle datasets without time")
        raise

    # Loop over time dimension for requested var
    var = ds.variables[varname]
    time = ds.variables['time']
    nt = len(time[:])
    filled_c = ax.tricontourf(lon, lat, var[1,:],
                   transform=ccrs.PlateCarree(),alpha=1.0)
    fig.colorbar(filled_c, orientation="horizontal")

    def animate (n):
        data = var[n,:]
        filled_c = ax.tricontourf(lon, lat, data,
                       transform=ccrs.PlateCarree(),alpha=0.8)

    ani = animation.FuncAnimation(fig,animate,nt,interval=100,repeat=False)
    writer = animation.writers['ffmpeg'](fps=10)
    ani.save('temp.mp4',writer=writer,dpi=dpi)

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

