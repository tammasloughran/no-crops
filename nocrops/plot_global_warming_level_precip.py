#!/usr/bin/env python3
# Plot the difference of no-crops and with crops for the global warmining level simulations.
import glob
import os

import scipy.stats as stats
import cartopy.crs as ccrs
import cdo_decorators as cdod
import ipdb
import matplotlib.pyplot as plt
import matplotlib as mpl
import netCDF4 as nc
import numpy as np
from cdo import Cdo
import warnings
warnings.filterwarnings('ignore')

cdo = Cdo()
cdo.debug = False

VARIABLES = {
        #'tas':'fld_s03i236',
        'pr':'fld_s05i216',
        }
EXPERIMENTS = {
        'PI-GWL-t6':'GWL-NoCrops-B2030',
        'PI-GWL-B2035':'GWL-NoCrops-B2035',
        'PI-GWL-B2040':'GWL-NoCrops-B2040',
        'PI-GWL-B2045':'GWL-NoCrops-B2045',
        'PI-GWL-B2050':'GWL-NoCrops-B2050',
        'PI-GWL-B2055':'GWL-NoCrops-B2055',
        'PI-GWL-B2060':'GWL-NoCrops-B2060',
        #'PI-GWL-B2060_duplicate':'GWL-EqFor-B2060',
        }
ARCHIVE_DIR = '/g/data/p66/tfl561/archive_data'
RAW_CMIP_DIR = '/g/data/p73/archive/non-CMIP/ACCESS-ESM1-5'
RAW_NOCROP_DIR = '/g/data/p66/tfl561/ACCESS-ESM'
LAND_FRAC_CODE = 'fld_s03i395'
TILE_FRAC_CODE = 'fld_s03i317'
EXAMPLE_FLIE = f'{ARCHIVE_DIR}/GWL-NoCrops-B2030/cLeaf_GWL-NoCrops-B2030_0500-0700.nc'
G_IN_PG = 10**15
DPI = 200
LAST30 = str(-30*12)
COLORS = {
        'GWL-NoCrops-B2030':'#62EA00',
        'GWL-NoCrops-B2035':'#24CC00',
        'GWL-NoCrops-B2040':'#079F2A',
        'GWL-NoCrops-B2045':'#00786B',
        'GWL-NoCrops-B2050':'#055992',
        'GWL-NoCrops-B2055':'#1140AB',
        'GWL-NoCrops-B2060':'#1E31B6',
        'GWL-EqFor-B2060':'deepskyblue',
        }


def moving_average(x:np.ndarray, window:int)->np.ndarray:
    """Calculate a moving average using a window of size `window`.
    """
    return np.convolve(x, np.ones(window), mode='same')/window


def yearly_mean_from_monthly(data:np.ndarray)->np.ndarray:
    """Calculate a yearly mean on a numpy array of monthly data.
    The 0th dimension must be time and divisible by 12.
    """
    if np.mod(data.shape[0], 12)!=0:
        raise ValueError("Not enough months in 0th dimension.")
    toshape = list(data.shape)
    toshape.pop(0)
    toshape.insert(0, 12)
    toshape.insert(0, int(data.shape[0]/12))
    fraction_of_year = np.array([31,28,31,30,31,30,31,31,30,31,30,31])/365.0
    return np.average(data.reshape(toshape), axis=1, weights=fraction_of_year)


def plot_last30(data_to_plot, title, save):
    fig2 = plt.figure()
    ax = fig2.add_subplot(1, 1, 1, projection=ccrs.Robinson())
    discrete_bins = mpl.colors.BoundaryNorm(boundaries=np.arange(-3.1, 3.2, 0.2), ncolors=256)
    colors = ax.pcolormesh(lons, lats, data_to_plot,
            cmap='seismic_r',
            norm=discrete_bins,
            transform=ccrs.PlateCarree(),
            )
    ax.coastlines()
    plt.colorbar(colors,
            label='$\Delta$ Precip. (mm/day)',
            ticks=[-3,-2.5,-2,-1.5,-1,-.5,0,.5,1,1.5,2,2.5,3],
            orientation='horizontal',
            pad=0.05,
            )
    plt.title(title)
    plt.savefig(f'plots/{save}', dpi=DPI)


def plot_australia_last30(data_to_plot, title, save):
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    discrete_bins = mpl.colors.BoundaryNorm(boundaries=np.arange(-3.1, 3.2, 0.2), ncolors=256)
    australia_lonlat = [110,155,-45,-10]
    ax.set_extent(australia_lonlat, crs=ccrs.PlateCarree())
    colors = ax.pcolormesh(lons, lats, data_to_plot,
            cmap='seismic_r',
            norm=discrete_bins,
            transform=ccrs.PlateCarree(),
            )
    ax.coastlines()
    plt.colorbar(colors,
            label='$\Delta$ Precip. (mm/day)',
            ticks=[-3,-2.5,-2,-1.5,-1,-.5,0,.5,1,1.5,2,2.5,3],
            orientation='horizontal',
            pad=0.05,
            )
    plt.title(title)
    plt.savefig(f'plots/{save}', dpi=DPI)


def apply_time_series_to_plot(years, line):
    smooth = moving_average(line, window=10)
    plt.figure(1)
    plt.plot(years, line,
            color=COLORS[nocrop_exp],
            alpha=0.3,
            )
    plt.plot(years, smooth,
            color=COLORS[nocrop_exp],
            label=exp,
            )
    plt.hlines(y=0, xmin=years[0], xmax=years[-1], color='black')
    plt.xlabel('Time (years)')
    plt.ylabel('$\Delta$ Precip. (mm/day)')
    plt.xlim(left=400)
    plt.legend(frameon=False)


load_from_npy = True

# Load lats and lons.
ncfile = nc.Dataset(EXAMPLE_FLIE, 'r')
lats = ncfile.variables['lat'][:]
lons = ncfile.variables['lon'][:]

data = {}
data_tmean = {}
plt.figure()
for gwl_exp,nocrop_exp in EXPERIMENTS.items():
    if 'duplicate' in gwl_exp: gwl_exp = 'PI-GWL-B2060'
    for exp in [gwl_exp, nocrop_exp]:
        print(f"{exp}")
        data[exp] = {}
        data_tmean[exp] = {}


        @cdod.cdo_mulc('86400')
        @cdod.cdo_fldmean(weights='TRUE') # Global mean.
        def load_global_sum(var, input:str)->np.ma.MaskedArray:
            """Load global sum of a carbon pool variable.
            """
            ncfile = cdo.copy(input=input, returnCdf=True, options='-L')
            return ncfile.variables[var][:]


        @cdod.cdo_mulc('86400')
        @cdod.cdo_seltimestep(f'{LAST30}/-1') # last 30 years
        @cdod.cdo_timmean # Temporal average
        def load_last30(var, input:str)->np.ma.MaskedArray:
            """Load last 30 year mean of a carbon pool variable.
            """
            ncfile = cdo.copy(input=input, returnCdf=True, options='-L')
            return ncfile.variables[var][:]


        # Load the data.
        period = '0500-0700'
        if load_from_npy:
            for var in VARIABLES.keys():
                data[exp][var] = np.load(f'data/{var}_{exp}_global_mean.npy')
                data_tmean[exp][var] = np.load(f'data/{var}_{exp}_last30.npy')
        else:
            for var in VARIABLES.keys():
                print(f"Loading {var}")
                data[exp][var] = load_global_sum(var,
                        input=f'{ARCHIVE_DIR}/{exp}/{var}_{exp}_{period}.nc',
                        )
                np.save(f'data/{var}_{exp}_global_mean.npy', data[exp][var].data)
                data_tmean[exp][var] = load_last30(var,
                        input=f'{ARCHIVE_DIR}/{exp}/{var}_{exp}_{period}.nc',
                        )
                np.save(f'data/{var}_{exp}_last30.npy', data_tmean[exp][var].data)

    # Plot the difference time series of temperature.
    if len(data[nocrop_exp]['pr'])==1223:
        data[nocrop_exp]['pr'] = np.append(data[nocrop_exp]['pr'], [np.nan])
    line = data[nocrop_exp]['pr'].squeeze() - data[gwl_exp]['pr'].squeeze()
    nmonths = data[exp]['pr'].squeeze().shape[0]
    years = np.linspace(400, 400+nmonths*(1/12), nmonths)
    line = yearly_mean_from_monthly(line)
    years = yearly_mean_from_monthly(years)
    apply_time_series_to_plot(years, line)

    # Plot the map of the difference for the last 30 years
    data_to_plot = data_tmean[nocrop_exp]['pr'].squeeze() - data_tmean[gwl_exp]['pr'].squeeze()
    plot_last30(data_to_plot, f'{nocrop_exp} - {gwl_exp}', f'pr_{nocrop_exp}_last30.png')

    # Plot the map of the difference for the last 30 years in Australia.
    plot_australia_last30(
            data_to_plot,
            f'{nocrop_exp} - {gwl_exp}',
            f'pr_{nocrop_exp}_australia_last30.png',
            )

# Significance testing of all experiments. I have tried standard ttest, kolmogorov-smirnov test
# and Wilcoxon test. I think Wilcoxon test is the most appropriate one for paired data.
shape = data_tmean[nocrop_exp]['pr'].squeeze().shape
data_to_plot = np.ones((len(EXPERIMENTS.keys()),)+shape)*np.nan
dataa = np.ones((len(EXPERIMENTS.keys()),)+shape)*np.nan
datab = np.ones((len(EXPERIMENTS.keys()),)+shape)*np.nan
i = 0
for gwl_exp,nocrop_exp in EXPERIMENTS.items():
    if 'duplicate' in gwl_exp: gwl_exp = 'PI-GWL-B2060'
    data_to_plot[i] = data_tmean[nocrop_exp]['pr'].squeeze() - data_tmean[gwl_exp]['pr'].squeeze()
    dataa[i] = data_tmean[nocrop_exp]['pr'].squeeze()
    datab[i] = data_tmean[gwl_exp]['pr'].squeeze()
    i += 1
statistic, pvalue = stats.wilcoxon(dataa, datab, method='exact', axis=0)
data_to_plot = data_to_plot.mean(axis=0)
# Mask out regions the exceed the 5% significance threshold.
data_to_plot[pvalue>0.05] = np.nan

# Plotting
plot_last30(data_to_plot, 'Mean of all experiments', 'pr_all_experiments_last30_sig.png')

plot_australia_last30(
        data_to_plot,
        'Mean of all experiments',
        'pr_all_experiments_australia_last30_sig.png',
        )

# Save figures.
plt.figure(1)
plt.savefig(f'plots/pr_GWL_gloabl_mean.png', dpi=DPI)
plt.show()

