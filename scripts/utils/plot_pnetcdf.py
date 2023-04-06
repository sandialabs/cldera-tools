# Import libraries
#import cartopy.crs as ccrs
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.patches import Rectangle
from scipy.interpolate import griddata
import xarray

###################################################################################################
def animate(i):
    global ax, cb
    so2 = ds['Tmean'].isel(time=i+1)
    so2_max = so2.max(dim='lev', keep_attrs=True)
    cmap = plt.colormaps['twilight']
    h = plt.tricontourf(ds['lon'], ds['lat'], so2_max, cmap=cmap, transform=ccrs.PlateCarree())
    if i > 0 and (i % 10 == 0):
        cb.remove()
        cb = plt.colorbar(h, orientation='horizontal', label=f'{so2.attrs["long_name"]}')
    return h

###################################################################################################
if __name__ == "__main__":
    '''
    Simple script to visualize output
    '''
    # Point to nc file that contains the value
    nc_file="/home/gbharpe/Programming/cpp/CLDERA/cldera_HSW_Tstats.h2.nc"
    ds = xarray.open_mfdataset(nc_file)

    # Test SO2 animation with cartopy
    global ax, cb
    fig = plt.figure(figsize=(8,4))
    ax = fig.add_subplot(111, projection=ccrs.PlateCarree(central_longitude=180))
    ax.coastlines()
    ax.set_xlim(-180,180)
    h = animate(0)
    so2 = ds['Tmean']
    cb = plt.colorbar(h, orientation='horizontal', label=f'{so2.attrs["long_name"]}')
    anim = animation.FuncAnimation(fig, animate, frames=120, interval=20)
    name = 'Tmean'
    #anim.save(name + '.mp4', fps=10, extra_args=['-vcodec', 'libx264'])
    plt.show()
