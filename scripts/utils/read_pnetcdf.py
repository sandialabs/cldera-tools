#!/usr/bin/env python3

from netCDF4 import Dataset
import numpy as np

import pathlib
import math
import matplotlib.pyplot as plt

# filename
#filename = "/home/gbharpe/Programming/cpp/CLDERA/cldera-tools/data/sample_database.nc"
#filename = "/home/gbharpe/Programming/cpp/CLDERA/cldera-tools/data/cldera_Tstats.nc"
ref_filename = "/home/gbharpe/Programming/cpp/CLDERA/cldera_HSW_Tstats.h2.nc"
#filename = "/home/gbharpe/Programming/cpp/CLDERA/cldera_stats.nc"
#filename = "/ascldap/users/gbharpe/Programming/cpp/CLDERA/HSW_ne16pg2_pnetcdf_reference_data/cldera_stats.nc"
#filename = "/home/gbharpe/Programming/cpp/CLDERA/E3SM_ne16pg2_ne16pg2_L72_FIDEAL-CLDERA_SAI.eam.h1.0001-01-01-00000.nc"
filename = "/home/gbharpe/Programming/cpp/CLDERA/E3SM_ne16pg2_ne16pg2_L72_FIDEAL-CLDERA_SAI.eam.h1.0001-06-30-00000.nc"

def print_available_variables(nc_filename):
  ds = Dataset(nc_filename,'r')

  # print dimension info
  print("Dimensions")
  for key in ds.dimensions:
    print(key)
    print(ds.dimensions[key])

  ntime = ds.dimensions["time"].size
  # time scale factor is 86400
  nlev = ds.dimensions["lev"].size
  ncol = ds.dimensions["ncol"].size

  # print variable info
  print("Variables")
  for key in ds.variables:
    print(key)
    print(ds.variables[key])

  # print more detailed info for temperature
  lat = ds.variables["lat"]
  lon = ds.variables["lon"]
  time = ds.variables["time"]

  #for icol in range(0,ncol):
  #  print("{}, {}".format(lat[icol],lon[icol]))
  print(np.unique(lat))

############################################
# main
############################################

#print_available_variables(filename)

ilev = 11
time_ratio = 4 # comparing 6-hourly differences to 2-day means and deviations

ref_ds = Dataset(ref_filename,'r')
run_ds = Dataset(filename,'r')

n_reftime = ref_ds.dimensions["time"].size
n_runtime = run_ds.dimensions["time"].size
ncol = ref_ds.dimensions["ncol"].size

i_runtime = 0

latbounds = [-20,20]

lonarray = np.array(run_ds.variables["lon"])
latarray = np.array(run_ds.variables["lat"])

# there should be 151 2-day averages
# but the other file only has 481 6-hourly instantaneous
for i_reftime in range(0,120): 
  Tmean_slice = ref_ds.variables["Tmean"][i_reftime+120,ilev,:]
  Tstd_slice = ref_ds.variables["Tstd"][i_reftime+120,ilev,:]
  # take 4 steps before updating the other
  for j in range(0,time_ratio): 
    T_slice = run_ds.variables["T"][i_runtime,ilev,:]
    # compute difference
    deviations = (T_slice - Tmean_slice)/Tstd_slice
    deviations = np.multiply(deviations,run_ds.variables["area"])
    deviations = np.where(np.isfinite(deviations), deviations, 0)
    h = plt.tricontourf(lonarray, latarray, deviations)
    cb = plt.colorbar(h, orientation='horizontal', label='')
    plt.xlabel("lon")
    plt.ylabel("lat")
    plt.title("Deviation from reference on day {}, 1991".format(run_ds.variables["time"][i_runtime]))
    plt.savefig("eruption_"+str(i_runtime).zfill(4))
    cb.remove()
    # isolate for lat bounds and accumulate into zonal deviation
    #zonal_deviation = 0
    #total_area = 1
#    for icol in range(0,ncol):
#      latval = run_ds.variables["lat"][icol]
#      if latval > latbounds[0] and latval < latbounds[1] and not math.isnan(deviations[icol]):
#        area = run_ds.variables["area"][icol]
#        total_area = total_area + area
#        zonal_deviation = zonal_deviation + (deviations[icol]*area)*(deviations[icol]*area)
#        print(deviations[icol])

    # increment runtime
    #print("Time: {}, Anomaly score: {}".format(run_ds.variables["time"][i_runtime],math.sqrt(zonal_deviation)/total_area))
    i_runtime = i_runtime+1

# for itime in range(0,ntime):
#   print(time[itime]*86400,end=' ')
#for ilev in range(0,nlev):
#  for icol in range(0,ncol):
#    print(T[9,ilev,icol],end=', ')
    # print("{}, {}, {}".format(T[0,0,icol],lat[icol],lon[icol]))

