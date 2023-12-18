#!/usr/bin/env python3
# Plot the difference of no-crops and with crops for the global warmining level simulations.
import glob
import os

import scipy.stats as stats
import cartopy.crs as ccrs
import cdo_decorators as cdod
import ipdb
import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np
from cdo import Cdo
import warnings
warnings.filterwarnings('ignore')

cdo = Cdo()
cdo.debug = True

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
        'PI-GWL-B2060_duplicate':'GWL-EqFor-B2060',
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


import numpy as _np
def yearly_mean_from_monthly(data:_np.ndarray)->_np.ndarray:
    """Calculate a yearly mean on a numpy array of monthly data.
    The 0th dimension must be time and divisible by 12.
    """
    if _np.mod(data.shape[0], 12)!=0:
        raise ValueError("Not enough months in 0th dimension.")
    toshape = list(data.shape)
    toshape.pop(0)
    toshape.insert(0, 12)
    toshape.insert(0, int(data.shape[0]/12))
    fraction_of_year = _np.array([31,28,31,30,31,30,31,31,30,31,30,31])/365.0
    return _np.average(data.reshape(toshape), axis=1, weights=fraction_of_year)


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


        @cdod.cdo_seltimestep(f'{LAST30}/-1') # last 30 years
        @cdod.cdo_timmean # Temporal average
        def load_last30(var, input:str)->np.ma.MaskedArray:
            """Load last 30 year mean of a carbon pool variable.
            """
            ncfile = cdo.copy(input=input, returnCdf=True, options='-L')
            return ncfile.variables[var][:]


        # Load the data.
        if 'B2030' in exp or 't6' in exp:
            period = '0500-0700'
        else:
            period = '0500-0601'
        if load_from_npy:
            for var in VARIABLES.keys():
                data[exp][var] = np.load(f'data/{var}_{exp}_global_sum.npy') # [g(C)]
                data_tmean[exp][var] = np.load(f'data/{var}_{exp}_last20.npy')
        else:
            for var in VARIABLES.keys():
                print(f"Loading {var}")
                #data[exp][var] = load_global_sum(var,
                #        input=f'{ARCHIVE_DIR}/{exp}/{var}_{exp}_{period}.nc',
                #        )
                #np.save(f'data/{var}_{exp}_global_sum.npy', data[exp][var].data) # [g(C)]
                #np.save(f'data/{var}_{exp}_last20.npy', data_tmean[exp][var].data) # [g(C)]
                data_tmean[exp][var] = load_last30(var,
                        input=f'{ARCHIVE_DIR}/{exp}/{var}_{exp}_{period}.nc',
                        )
                np.save(f'data/{var}_{exp}_firstyear.npy', data_tmean[exp][var].data) # [g(C)]

    # Plot the difference time series of temperature.
    plt.figure(1)
    if len(data[nocrop_exp]['tas'])==1223:
        data[nocrop_exp]['tas'] = np.append(data[nocrop_exp]['tas'], [np.nan])
    line = data[nocrop_exp]['tas'].squeeze() - data[gwl_exp]['tas'].squeeze()
    nmonths = data[exp]['tas'].squeeze().shape[0]
    years = np.linspace(400, 400+nmonths*(1/12), nmonths)
    line = yearly_mean_from_monthly(line)
    years = yearly_mean_from_monthly(years)
    plt.plot(years, line,
            color=COLORS[nocrop_exp],
            label=exp,
            )
    plt.xlabel('Time (years)')
    plt.ylabel('$\Delta$ Temperature ($^{\circ}$C)')
    plt.xlim(left=400)
    #plt.ylim(bottom=0)
    plt.legend(frameon=False)

    # Plot the map of the difference for the last 30 years
    difference = data_tmean[nocrop_exp]['tas'].squeeze() - data_tmean[gwl_exp]['tas'].squeeze()
    fig2 = plt.figure()
    ax = fig2.add_subplot(1, 1, 1, projection=ccrs.Robinson())
    colors = ax.pcolormesh(lons, lats, difference,
            cmap='seismic',
            vmin=-4,
            vmax=4,
            transform=ccrs.PlateCarree(),
            )
    ax.coastlines()
    plt.colorbar(colors,
            label='$\Delta$ temperature ($^{\circ}$C)',
            orientation='horizontal',
            pad=0.05,
            )
    plt.title(f'{nocrop_exp} - {gwl_exp}')
    plt.savefig(f'plots/tas_{nocrop_exp}_last30.png', dpi=DPI)

    # Plot the map of the difference for the last 30 years in Australia.
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    australia_lonlat = [110,155,-45,-10]
    ax.set_extent(australia_lonlat, crs=ccrs.PlateCarree())
    data_to_plot = data_tmean[nocrop_exp]['tas'].squeeze() - data_tmean[gwl_exp]['tas'].squeeze()
    colors = ax.pcolormesh(lons, lats, data_to_plot,
            cmap='seismic',
            vmin=-4,
            vmax=4,
            transform=ccrs.PlateCarree(),
            )
    ax.coastlines()
    plt.colorbar(colors,
            label='$\Delta$ Surface air temperature ($^{\circ}$C)]',
            orientation='horizontal',
            pad=0.05,
            )
    plt.title(f'{nocrop_exp} - {gwl_exp}')
    plt.savefig(f'plots/tas_{nocrop_exp}_australia_last30.png', dpi=DPI)

# Significance testing of all experiments.
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
t_statistic, pvalue = stats.ttest_ind(dataa, datab, axis=0, equal_var=False)
data_to_plot = data_to_plot.mean(axis=0)
data_to_plot[pvalue>0.05] = np.nan
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
colors = ax.pcolormesh(lons, lats, data_to_plot.squeeze(),
        cmap='seismic',
        vmin=-4,
        vmax=4,
        transform=ccrs.PlateCarree(),
        )
ax.coastlines()
plt.colorbar(colors,
        label='$\Delta$ Surface air temperature ($^{\circ}$C)]',
        orientation='horizontal',
        pad=0.05,
        )
plt.savefig(f'plots/tas_{nocrop_exp}_last30_sig.png', dpi=DPI)

# Save figures.
plt.figure(1)
plt.savefig(f'plots/tas_GWL_gloabl_mean.png', dpi=DPI)
plt.show()

