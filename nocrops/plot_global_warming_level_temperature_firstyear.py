#!/usr/bin/env python3
# Plot the difference of no-crops and with crops for the global warmining level simulations.
import glob
import os

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

data_tmean = {}
plt.figure()
for gwl_exp,nocrop_exp in EXPERIMENTS.items():
    if 'duplicate' in gwl_exp: gwl_exp = 'PI-GWL-B2060'
    for exp in [gwl_exp, nocrop_exp]:
        print(f"{exp}")
        data_tmean[exp] = {}


        #@cdod.cdo_seltimestep(f'{LAST100}/-1') # last 100 years
        @cdod.cdo_seltimestep(f'1/12') # last year
        @cdod.cdo_timmean # Temporal average
        def load_last100(var, input:str)->np.ma.MaskedArray:
            """Load last 100 year mean of a carbon pool variable.
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
                data_tmean[exp][var] = np.load(f'data/{var}_{exp}_firstyear.npy')
        #for var in VARIABLES.keys():
            #data_tmean[exp][var] = load_last100(
            #        var,
            #        input=f'{ARCHIVE_DIR}/{exp}/{var}_{exp}_{period}.nc',
            #        )
            #np.save(f'data/{var}_{exp}_last100.npy', data_tmean[exp][var].data) # [g(C)]
            #np.save(f'data/{var}_{exp}_firstyear.npy', data_tmean[exp][var].data) # [g(C)]

    # Plot the map of the difference for the last 100 years
    fig2 = plt.figure()
    ax = fig2.add_subplot(1, 1, 1, projection=ccrs.Robinson())
    colors = ax.pcolormesh(lons, lats,
            data_tmean[nocrop_exp]['tas'].squeeze() - data_tmean[gwl_exp]['tas'].squeeze(),
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
    plt.savefig(f'plots/tas_{nocrop_exp}_firstyear.png', dpi=DPI)

    # Plot the map of the difference for the last 100 years in Australia.
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
    plt.savefig(f'plots/tas_{nocrop_exp}_australia_firstyear.png', dpi=DPI)

# Aggregate
ensemble_mean = 0
n = 0
for gwl_exp,nocrop_exp in EXPERIMENTS.items():
    if 'duplicate' in gwl_exp: gwl_exp = 'PI-GWL-B2060'
    ensemble_mean += data_tmean[nocrop_exp]['tas'].squeeze() - data_tmean[gwl_exp]['tas'].squeeze()
    n += 1
ensemble_mean /= n

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1, projection=ccrs.Robinson())
colors = ax.pcolormesh(lons, lats, ensemble_mean,
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
plt.title(f'First year experiment mean')
plt.savefig(f'plots/tas_{nocrop_exp}_firstyear_ensemble_mean.png', dpi=DPI)
#plt.show()
