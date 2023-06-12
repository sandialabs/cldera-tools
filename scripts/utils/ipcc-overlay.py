#!/usr/bin/env python3

import argparse, sys, pathlib
import netCDF4 as nc
import numpy as np
import regionmask as rm
import cartopy.crs as ccrs
import matplotlib.animation as animation
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import ssl
from utils import expect

ssl._create_default_https_context = ssl._create_unverified_context

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
                        help="whether lat/lon in input file are in degrees")
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
    parser.add_argument("--time-start",default=None,type=int,
                        help="When making a video, start from this time")
    parser.add_argument("--time-end",default=None,type=int,
                        help="When making a video, stop at this time")
    parser.add_argument("--time-freq",default=None,type=int,
                        help="When making a video, sample every this many time steps")
    parser.add_argument("-r","--resolution", default=10, type=int,
                        help="Number of levels to use in the contour plot")
    parser.add_argument("--log-scale", action='store_true',
                        help="Use a log scale in the plot")
    parser.add_argument("--cminmax", default=None, nargs='+',
                        help="Min/max value for colorbar in contour")
    parser.add_argument("varname",
                        help="Variable to plot")

    return parser.parse_args(args[1:])

###############################################################################
def plot_over_ipcc_regions(input_filename,degrees,lev,output_filename,mask_file,mask_name,time,time_start,time_end,time_freq,resolution,log_scale,cminmax,varname):
###############################################################################

    dpi = 100
    expect (resolution>0,f"Invalid number of levels ({resolution}). Value must be positive")

    # Create figure, and set some static info/data
    fig = plt.figure()
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
    expect ((time_start is None and time_end is None) or has_time,
            f"Flag --time-start and/or --time-end were used for non-time-dependent variable '{varname}{dims}'.",
            ValueError)
    expect (time is None or (time_start is None and time_end is None and time_freq is None),
            f"Cannot specify a single time slice as well as a time interval")
    if has_time:
        time_var = ds.variables['time']
        nt = len(time_var[:])
        time_dim = dims.index('time')
        expect(time_dim==0,f"We expect time to be the slowest dim, but layout is {dims}",ValueError)
        expect (time is None or (time>=0 and time<nt),f"Input time slice ({time}) out of bounds: [0,{nt})")
        expect (time_start is None or (time_start>=0 and time_start<nt),
                f"Input time start ({time_start}) out of bounds: [0,{nt})")
        expect (time_end is None or (time_end>=0 and time_end<nt),
                f"Input time end ({time_end}) out of bounds: [0,{nt})")
        expect (time_freq is None or time_freq>0,
                f"Time freq should be positive (got {time_freq} instead)")
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

    # Deduce the time window specs 
    time_start = (0 if time is None else time) if time_start is None else time_start
    time_end = (nt if time is None else time) if time_end is None else time_end
    time_freq = 1 if time_freq is None else time_freq
    num_frames = int((time_end - time_start) / time_freq) + 1
    print (f"time_start={time_start},time_end={time_end},time_freq={time_freq},num_frames={num_frames}")

    # Helper function, to get slice of data at given time index
    def get_data_at_time(n):
        data = var[:]
        if has_lev:
            data = np.take(data,lev,lev_dim)
        if has_time:
            data = np.take(data,n,time_dim)
        return data

    # Get the global max/min of data over the whole time interval,
    # so we can set the contour
    if cminmax is None:
        gmax = -np.Inf
        gmin =  np.Inf

        for f in range(0,num_frames):
            t = time_start + f*time_freq
            gmax = max(gmax,np.amax(get_data_at_time(t)))
            gmin = min(gmin,np.amin(get_data_at_time(t)))

        if log_scale:
            if gmin==0:
                gmin = 1e-40 if gmax==0 else 1e-40*gmax
            if gmax==0:
                gmax = 1e-40
            gmax = np.log10(gmax)
            gmin = np.log10(gmin)
    else:
        expect (len(cminmax)==2, "When using --cminmax, you must provide exactly two values")
        gmin = float(cminmax[0])
        gmax = float(cminmax[1])
        expect (gmax>=gmin, "The arguments to --cminmax must be in increasing order")

    # This fcn "massages" the data (clips, takes log) and plot contour
    def do_contour(data,set_colorbar=False):
        dmax = np.amax(data[:])
        dmin = np.amin(data[:])

        # If log scale, actually compute the log of the data
        if log_scale:
            if dmin<=0:
                for i in range(0,len(data)):
                    if data[i]<=0:
                        data[i] = 1e-40
                clip_min = 1e-40 if dmax==0 else 1e-40*dmax
                clip_max = 1e-40 if dmax==0 else dmax
                print (f"WARNING! Data range is [{dmin}, {dmax}]. Clipping to [{clip_min},{clip_max}]")
                data = np.clip(data,clip_min,clip_max)

                # Recompute bounds
                dmin = np.amin(data)
                dmax = np.amax(data)

            # Convert to log
            data = np.log10(data)

            # Recompute bounds
            dmin = np.amin(data)
            dmax = np.amax(data)

        nlevels = resolution if resolution<len(data) else len(data)

        try:
            #  print (f"gmax={gmax},gmin={
            levels = np.linspace(gmin,gmax,nlevels)
            filled_c = ax.tricontourf(lon, lat, data[:],
                           transform=ccrs.PlateCarree(),alpha=1.0,
                           levels=levels)
            if set_colorbar:
                cbar = plt.colorbar(filled_c, orientation="horizontal") 

                # For the colorbar, 5 ticks are enough
                levels = np.linspace(gmin,gmax,5)
                cbar.set_ticks(levels)
                cbar.set_ticklabels(["{:15.3f}".format(s) for s in levels])
                ticks = cbar.get_ticks()
                if log_scale:
                    cbar_label = f"log({varname})"
                else:
                    cbar_label = varname
                cbar.set_label(cbar_label)

        except ValueError as e:
            if str(e)!="Contour levels must be increasing":
                raise

    # I can't find a way to make FuncAnimation skip the call to init_func.
    # So provide a do-nothing func instead
    def ani_init ():
        return

    time_text = ax.text(0.5, 1.2,'',horizontalalignment='center',verticalalignment='top', transform=ax.transAxes)

    if has_ncol:
        # Normal plot of ncol values
        def plot_frame (n):
            tt = time_var[time_start+n*time_freq]
            print (f"plotting frame at time = {tt}")
            data = get_data_at_time(time_start + n*time_freq)
            do_contour(data,n==0)
            time_text.set_text("time = {:7.3f} days".format(tt))

        if time is not None or not has_time:
            if time is not None:
                time_text.set_text("time = {:4.3f} days".format(time_var[time]))
            do_contour(get_data_at_time(time_start),True)
            plt.savefig("fig.png" if output_filename is None else output_filename)
        else:
            ani = animation.FuncAnimation(fig,plot_frame,num_frames,interval=100,repeat=False,init_func=ani_init)
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
        def plot_frame (n):
            tt = time_var[time_start+n*time_freq]
            print (f"plotting frame at time = {tt}")
            data = get_data_at_time(time_start + n*time_freq)
            data2d = create_data_2d(data)
            do_contour(data2d,n==0)
            time_text.set_text("time = {:7.3f} days".format(tt))

        if time is not None or not has_time:
            if time is not None:
                time_text.set_text("time = {:7.3f} days".format(time_var[time]))
            do_contour(create_data_2d(get_data_at_time(time_start)),True)
            plt.savefig("fig.png" if output_filename is None else output_filename)
        else:
            ani = animation.FuncAnimation(fig,plot_frame,num_frames,interval=100,repeat=False,init_func=ani_init)
            writer = animation.writers['ffmpeg'](fps=10)
            ani.save("movie.mp4" if output_filename is None else output_filename,writer=writer,dpi=dpi)

    # Uncomment to visualize without having to open a fig/mov file
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

