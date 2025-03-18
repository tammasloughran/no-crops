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
        'tas':'fld_s03i236',
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
LAST100 = str(-100*12)
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


def my_ks_test(samp1:np.ndarray, samp2:np.ndarray)->tuple:
    """Custom 2-dimensional Kolmogorov-Smirnov test for the goodness of fit of 2 populations.
    """
    nlats = samp1.shape[-2]
    nlons = samp1.shape[-1]
    statistic = np.ones((nlats,nlons))*np.nan
    pvalue = np.ones((nlats,nlons))*np.nan
    for j in range(nlats):
        for i in range(nlons):
            results = stats.ks_2samp(samp1[:,j,i], samp2[:,j,i])
            statistic[j,i] = results.statistic
            pvalue[j,i] = results.pvalue
    return statistic, pvalue


def plot_globe(data:np.ndarray, title:str, save:str)->None:
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.Robinson())
    discrete_bins = mpl.colors.BoundaryNorm(boundaries=np.arange(-2.1, 2.2, 0.2), ncolors=256)
    colors = ax.pcolormesh(lons, lats, data,
            cmap='seismic',
            norm=discrete_bins,
            transform=ccrs.PlateCarree(),
            )
    ax.coastlines()
    plt.colorbar(colors,
            label='$\Delta$ temperature ($^{\circ}$C)',
            ticks=[-2,-1.5,-1,-.5,0,.5,1,1.5,2],
            orientation='horizontal',
            pad=0.05,
            )
    plt.title(title)
    plt.savefig(f'plots/{save}', dpi=DPI)

def plot_australia(data:np.ndarray, title:str, save:str)->None:
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    australia_lonlat = [110,155,-45,-10]
    ax.set_extent(australia_lonlat, crs=ccrs.PlateCarree())
    discrete_bins = mpl.colors.BoundaryNorm(boundaries=np.arange(-2.1, 2.2, 0.2), ncolors=256)
    colors = ax.pcolormesh(lons, lats, data,
            cmap='seismic',
            norm=discrete_bins,
            transform=ccrs.PlateCarree(),
            )
    ax.coastlines()
    plt.colorbar(colors,
            label='$\Delta$ Surface air temperature ($^{\circ}$C)',
            orientation='horizontal',
            ticks=[-2,-1.5,-1,-.5,0,.5,1,1.5,2],
            pad=0.05,
            )
    plt.title(title)
    plt.savefig(f'plots/{save}', dpi=DPI)


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


        @cdod.cdo_fldmean(weights='TRUE') # Global sum.
        def load_global_sum(var, input:str)->np.ma.MaskedArray:
            """Load global sum of a carbon pool variable.
            """
            ncfile = cdo.copy(input=input, returnCdf=True, options='-L')
            return ncfile.variables[var][:]


        @cdod.cdo_seltimestep(f'{LAST100}/-1') # last 100 years
        @cdod.cdo_timmean # Temporal average
        def load_last100(var, input:str)->np.ma.MaskedArray:
            """Load last 100 year mean of a carbon pool variable.
            """
            ncfile = cdo.copy(input=input, returnCdf=True, options='-L')
            return ncfile.variables[var][:]


        # Load the data.
        period = '0500-0700'
        if load_from_npy:
            for var in VARIABLES.keys():
                data[exp][var] = np.load(f'data/{var}_{exp}_global_mean.npy')
                data_tmean[exp][var] = np.load(f'data/{var}_{exp}_last100.npy')
        else:
            for var in VARIABLES.keys():
                print(f"Loading {var}")
                data[exp][var] = load_global_sum(var,
                        input=f'{ARCHIVE_DIR}/{exp}/{var}_{exp}_{period}.nc',
                        )
                np.save(f'data/{var}_{exp}_global_mean.npy', data[exp][var].data)
                data_tmean[exp][var] = load_last100(var,
                        input=f'{ARCHIVE_DIR}/{exp}/{var}_{exp}_{period}.nc',
                        )
                np.save(f'data/{var}_{exp}_last100.npy', data_tmean[exp][var].data)

    # Plot the difference time series of temperature.
    plt.figure(1)
    if len(data[nocrop_exp]['tas'])==1223:
        data[nocrop_exp]['tas'] = np.append(data[nocrop_exp]['tas'], [np.nan])
    line = data[nocrop_exp]['tas'].squeeze() - data[gwl_exp]['tas'].squeeze()
    nmonths = data[exp]['tas'].squeeze().shape[0]
    years = np.linspace(400, 400+nmonths*(1/12), nmonths)
    line = yearly_mean_from_monthly(line)
    smooth = moving_average(line, window=50)
    years = yearly_mean_from_monthly(years)
    smooth_years = years[5:-5]
    plt.plot(years, line,
            color=COLORS[nocrop_exp],
            alpha=0.3,
            )
    #plt.plot(years[25:-25], smooth[50:-49],
    plt.plot(years[:-25], smooth[:-25],
            color=COLORS[nocrop_exp],
            label=exp[-4:],
            )
    plt.hlines(y=0, xmin=years[0], xmax=years[-1], color='black')
    plt.xlabel('Time (years)')
    plt.ylabel('$\Delta$ Temperature ($^{\circ}$C)')
    plt.xlim(left=400)
    plt.legend(frameon=False)

    # Plot the map of the difference for the last 100 years
    difference = data_tmean[nocrop_exp]['tas'].squeeze() - data_tmean[gwl_exp]['tas'].squeeze()
    plot_globe(difference, f'{nocrop_exp} - {gwl_exp}', f'tas_{nocrop_exp}_last100.png')
    # Plot the map of the difference for the last 100 years in Australia.
    plot_australia(
            difference,
            f'{nocrop_exp} - {gwl_exp}',
            f'tas_{nocrop_exp}_australia_last100.png',
            )

# Save the line plot
plt.figure(1)
plt.savefig(f'plots/tas_GWL_gloabl_mean.png', dpi=DPI)

# Significance testing of all experiments. I have tried standard ttest, kolmogorov-smirnov test
# and Wilcoxon test. I think Wilcoxon test is the most appropriate one for paired data.
shape = data_tmean[nocrop_exp]['tas'].squeeze().shape
data_to_plot = np.ones((len(EXPERIMENTS.keys()),)+shape)*np.nan
dataa = np.ones((len(EXPERIMENTS.keys()),)+shape)*np.nan
datab = np.ones((len(EXPERIMENTS.keys()),)+shape)*np.nan
i = 0
for gwl_exp,nocrop_exp in EXPERIMENTS.items():
    if 'duplicate' in gwl_exp: gwl_exp = 'PI-GWL-B2060'
    data_to_plot[i] = data_tmean[nocrop_exp]['tas'].squeeze() - data_tmean[gwl_exp]['tas'].squeeze()
    dataa[i] = data_tmean[nocrop_exp]['tas'].squeeze()
    datab[i] = data_tmean[gwl_exp]['tas'].squeeze()
    i += 1
statistic, pvalue = stats.wilcoxon(dataa, datab, method='exact', axis=0)
data_to_plot = data_to_plot.mean(axis=0)
# Mask out regions the exceed the 5% significance threshold.
data_to_plot[pvalue>0.05] = np.nan

# Plotting
plot_globe(data_to_plot, 'Mean of all experiments', 'tas_all_experiments_last100_sig.png')

plot_australia(
        data_to_plot,
        'Mean of all experiments',
        'tas_all_experiments_australia_last100_sig.png',
        )

# Show figures.
plt.show()

