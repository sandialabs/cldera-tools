#!/usr/bin/env python3
import glob
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
import sys

# Set fontsize and plotting text to latex
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=14)

###############################################################################
def load_files():
    # Open input files
    print("Loading files...")
    f_cf = sorted(glob.glob('/gpfs/cldera/data/E3SM-cldera-tools/CLDERA-historical/limvar/v2.LR.WCYCL20TR.pmcpu.limvar.ens6.cf/cldera_stats*'))
    f_er = sorted(glob.glob('/gpfs/cldera/data/E3SM-cldera-tools/CLDERA-historical/limvar/v2.LR.WCYCL20TR.pmcpu.limvar.ens6/cldera_stats*'))
    cf_data = xr.open_mfdataset(f_cf, concat_dim='time', combine='nested')
    er_data = xr.open_mfdataset(f_er, concat_dim='time', combine='nested')
    # f_cf = '/gpfs/cldera/data/E3SM-cldera-tools/CLDERA-historical/limvar/v2.LR.WCYCL20TR.pmcpu.limvar.ens6.cf/cldera_stats-1991-06-01-00000.nc'
    # f_er = '/gpfs/cldera/data/E3SM-cldera-tools/CLDERA-historical/limvar/v2.LR.WCYCL20TR.pmcpu.limvar.ens6/cldera_stats-1991-06-01-00000.nc'
    # cf_data = xr.open_dataset(f_cf)
    # er_data = xr.open_dataset(f_er)
    return cf_data, er_data


###############################################################################
def plot_diff(cf_data, er_data, var, reg):
    print("Plotting {}, region {:d}...".format(var,reg))
    num_time = cf_data.dims['time'] # cf dataset is smaller
    # avg_type = '1M'
    # cf_data = cf_data.set_index('Timestamp').resample(time=avg_type).mean(dim='time')
    # er_data = er_data.set_index('Timestamp').resample(time=avg_type).mean(dim='time')
    time = cf_data['time']
    cf = cf_data[var][:num_time, reg]
    er = er_data[var][:num_time, reg]
    var_diff = er - cf

    # Compute time average
    print('Computing time average...')
    avg_days = 1
    days_tot = 365*7
    # days_tot = 300
    steps_per_day = 48
    days = np.array(range(0, days_tot, avg_days))
    QOIavg_er = np.zeros(len(days))
    QOIavg_cf = np.zeros(len(days))
    for iday in range(len(days)):
        for istep in range(avg_days*steps_per_day):
            if (er[iday * avg_days*steps_per_day + istep] > 1e100 or cf[iday * avg_days*steps_per_day + istep] > 1e100):
                continue
            QOIavg_er[iday] += er[iday * avg_days*steps_per_day + istep]
            QOIavg_cf[iday] += cf[iday * avg_days*steps_per_day + istep]
        QOIavg_er[iday] /= float(avg_days*steps_per_day)
        QOIavg_cf[iday] /= float(avg_days*steps_per_day)
    QOI_diff = QOIavg_er - QOIavg_cf

    # Plot
    print("Creating plot...")
    fig, axs = plt.subplots(2, figsize=(12,12))
    # axs[0].plot(time, cf, linewidth=1.5, label='CF')
    # axs[0].plot(time, er, linewidth=1.5, label='ER')
    axs[0].plot(days, QOIavg_cf, linewidth=1.5, label='CF')
    axs[0].plot(days, QOIavg_er, linewidth=1.5, label='ER')
    axs[0].legend()
    axs[0].set_ylabel(var)
    axs[1].plot(days, QOI_diff, linewidth=1.5, label='Daily')
    axs[1].plot(time, var_diff, linewidth=1.5, alpha=0.25, label='Inst')
    axs[1].set_ylabel('Difference')
    axs[1].set_xlabel('Days')
    axs[1].legend()
    axs[1].set_ylabel(var)
    axs[0].set_title('Region ' + str(reg))
    for ax in axs:
        ax.set_xlim(0, days[-1])
        ax.set_axisbelow(True)
        ax.grid(linestyle='--', alpha=0.5)
    #plt.show()
    fig_name = 'CTOOLS_{}_r{:d}.png'.format(var,reg)
    plt.savefig(fig_name, dpi=300)

###############################################################################
if __name__ == "__main__":
    # Inputs
    #vars = ['AEROD_v_surf','FLNT_surf','T_50','FSDSC_surf','TREFHT_surf']
    vars = ['SO2_50',]
    #vars = ['AEROD_v_surf',]
    # vars = ['TREFHT_surf',]
    #regs = [38, 3, 44, 48, 55]
    regs = [38,]
    # regs = [48,]
    #regs = [55,]
    cf_data, er_data = load_files()
    for var in vars:
        for reg in regs:
            plot_diff(cf_data, er_data, var, reg)
