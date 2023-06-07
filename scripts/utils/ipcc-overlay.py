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
def expect(condition, error_msg, exc_type=SystemExit, error_prefix="ERROR:"):
###############################################################################
    """
    Similar to assert except doesn't generate an ugly stacktrace. Useful for
    checking user error, not programming error.

    >>> expect(True, "error1")
    >>> expect(False, "error2")
    Traceback (most recent call last):
        ...
    SystemExit: ERROR: error2
    """
    if not condition:
        msg = error_prefix + " " + error_msg
        raise exc_type(msg)

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
    parser.add_argument("-d","--degrees", action='store_true',
                        help="whether lat/lon in input file are in radians")
    parser.add_argument("-l","--lev", default=None, type=int,
                        help="If variable has a lev/ilev, slice it a this vertical level")
    parser.add_argument("-o","--output_filename", default=None,
                        help="Name of mp4 file to generate")
    parser.add_argument("-m","--mask-file", default=None,
                        help="Location of the nc file with masks for the regions where var is defined")
    parser.add_argument("-n","--mask-name", default=None,
                        help="Name of the mask variable in mask file")
    parser.add_argument("-t","--time",default=None,type=int,
                        help="Single time frame to use when producing an image")
    parser.add_argument("-r","--resolution", default=10, type=int,
                        help="Number of levels to use in the contour plot")
    parser.add_argument("varname",
                        help="Variable to plot")

    return parser.parse_args(args[1:])

###############################################################################
def plot_over_ipcc_regions(input_filename,degrees,lev,output_filename,mask_file,mask_name,time,resolution,varname):
###############################################################################

    dpi = 100
    expect (resolution>0,f"Invalid number of levels ({resolution}). Value must be positive")

    # Create figure, and set some static info/data
    fig = plt.figure()
    #  fig, ax = plt.subplots()
    ax = plt.axes(projection=ccrs.EckertIII())
    ax.coastlines()
    ax.set_global()
    ax.gridlines(xlocs=range(-180,180,60), ylocs=range(-90,90,30),draw_labels=True)
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
    if not degrees:
        print (" -> converting lat/lon from radians to degrees")
        lat = lat*180/np.pi
        lon = lon*180/np.pi

    min_lat = np.amin(lat);
    min_lon = np.amin(lon);
    if min_lon >= 0:
        print (f" -> min_lon: {min_lon}. Wrapping longitude around, to get it in [-180,180]")
        lon = ((lon+180) % 360) - 180
    if min_lat >= 0:
        print (f" -> min_lat: {min_lat}. Wrapping latitude around, to get it in [-90,90]")
        lon = ((lat+90) % 180) - 90

    # Get variable layout
    var = ds.variables[varname]
    dims = var.dimensions
    ndims = len(dims)
    expect(ndims<=3, f"Unsupported variable layout: {dims}.")
    has_time = 'time' in dims
    expect (time is None or has_time, f"Flag -t/--time was used for non-time-dependent variable '{varname}{dims}'.",ValueError)
    if has_time:
        time_dim = ds.variables['time']
        nt = len(time_dim[:])
        time_dim = dims.index('time')
        expect(time_dim==0,f"We expect time to be the slowest dim, but layout is {dims}",ValueError)
        expect (time is None or (time>=0 and time<nt),f"Input time slice ({time}) out of bounds: [0,{nt})")
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
        # with N=num_mask_labels. We *ASSUME* that cldera-tools creates a mask stat where the entries
        # are ordered based on the mask ids.
        ncols = dsm.dimensions['ncol'].size
        m2idx = {}
        mvals = []
        for m in mask[:]:
            if m not in mvals:
                mvals.append(m)
        mvals.sort()
        for m in mvals:
            m2idx[m] = mvals.index(m)
    has_lev = 'lev' in dims or 'ilev' in dims
    if has_lev:
        expect (lev is not None,
                "Variable has a vertical level dimension: {var.dimensions}. "
                "You must provide -l/--lev N argument.",ValueError)
        lev_dim = dims.index('lev') if 'lev' in dims else dims.index('ilev')
        has_lev = True

    # Helper function, to get slice of data at given time index
    def get_data_at_time(n):
        data = var[:]
        if has_lev:
            data = np.take(data,lev,lev_dim)
        if has_time:
            data = np.take(data,n,time_dim)
        return data

    def do_contour(data,set_colorbar=False):
        filled_c = ax.tricontourf(lon, lat, data,
                       transform=ccrs.PlateCarree(),alpha=1.0,
                       levels=np.linspace(np.amin(data)-0.0001,np.amax(data)+0.0001,resolution))
        if set_colorbar:
            fig.colorbar(filled_c, orientation="horizontal")

    if has_ncol:
        do_contour(get_data_at_time(0),True)
        def plot_at_time (n):
            data = get_data_at_time(n)
            do_contour(data)

        if time is not None or not has_time:
            plt.savefig("fig.png" if output_filename is None else output_filename)
        else:
            ani = animation.FuncAnimation(fig,plot_at_time,nt,interval=100,repeat=False)
            writer = animation.writers['ffmpeg'](fps=10)
            ani.save("movie.mp4" if output_filename is None else  output_filename,writer=writer,dpi=dpi)
    else:
        # We need to create a "heat map"
        def create_data_2d(data):
            data2d = np.zeros((ncols))
            for i in range(0,ncols):
                mval = int(mask[i])
                idx = m2idx[mval]
                data2d[i] = data[idx]
            return data2d
        do_contour(create_data_2d(get_data_at_time(0)),True)
        def plot_at_time (n):
            data = get_data_at_time(n)
            data2d = create_data_2d(data)
            do_contour(data2d)

        if time is not None or not has_time:
            plt.savefig("fig.png" if output_filename is None else output_filename)
        else:
            ani = animation.FuncAnimation(fig,plot_at_time,nt,interval=100,repeat=False)
            writer = animation.writers['ffmpeg'](fps=10)
            ani.save("movie.mp4" if output_filename is None else output_filename,writer=writer,dpi=dpi)

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

