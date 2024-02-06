#!/usr/bin/env python3
# Plot the carbon pools from the single pft no-crops forestation experiments.
import glob
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np
from cartopy import crs as ccrs

from cdo import Cdo; cdo = Cdo()
import cdo_decorators as cdod
[sys.path.append(i) for i in ['.', '..']]
from contrib.tools import global_mean

cdo.debug = False

EXPERIMENTS = [
        'PI-GWL-t6',
        'GWL-NoCrops-B2030',
        'GWL-50pct-B2030',
        'GWL-25pct-B2030',
        'GWL-10pct-B2030',
        ]
ARCHIVE_DIR = '/g/data/p66/tfl561/archive_data'
TASFIG = 1
MAPFIG = 2
TILE_FRAC_CODE = 'fld_s03i317'
RAW_CMIP_DIR = '/g/data/p73/archive/non-CMIP/ACCESS-ESM1-5'
RAW_NOCROP_DIR = '/g/data/p66/tfl561/ACCESS-ESM'
EXAMPLE = '/g/data/p66/tfl561/archive_data/GWL-NoCrops-B2040/tas_GWL-NoCrops-B2040_0500-0601.nc'
LABELS = {
        'PI-GWL-t6':'Global warming level',
        'GWL-NoCrops-B2030':'100%',
        'GWL-10pct-B2030':'10%',
        'GWL-25pct-B2030':'25%',
        'GWL-50pct-B2030':'50%',
        }
COLORS = {
        'PI-GWL-t6':'black',
        'GWL-NoCrops-B2030':'blue',
        'GWL-10pct-B2030':'lightgreen',
        'GWL-25pct-B2030':'forestgreen',
        'GWL-50pct-B2030':'darkgreen',
        }


@cdod.cdo_fldmean(weights='TRUE')
@cdod.cdo_yearmonmean
def cdo_load(var:str, input:str)->np.ma.MaskedArray:
    """Load global mean of a variable using CDO."""
    return cdo.copy(input=input, options='-L', returnCdf=True).variables[var][:].squeeze()


@cdod.cdo_selyear('570/599')
@cdod.cdo_timmean
def cdo_load_map(var:str, input:str)->np.ma.MaskedArray:
    """Load map of a variable using CDO."""
    return cdo.copy(input=input, options='-L', returnCdf=True).variables[var][:].squeeze()


def moving_average(x:np.ndarray, window:int)->np.ndarray:
    """Calculate a moving average using a window of size `window`.
    """
    return np.convolve(x, np.ones(window), mode='valid')/window


# Load grid data.
examplenc = nc.Dataset(EXAMPLE, 'r')
lats = examplenc.variables['lat'][:]
lons = examplenc.variables['lon'][:]
examplenc.close()

# Set up figures.
plt.figure(TASFIG)

data = {}
maps = {}
for exper in EXPERIMENTS:
    # Load the tas data
    print("Loading TAS", exper)
    files = sorted(glob.glob(f'{ARCHIVE_DIR}/{exper}/tas_{exper}_*'))

    # Load data
    data[exper] = cdo_load(input=files[0], var='tas')
    maps[exper] = cdo_load_map(input=files[0], var='tas')

    # Plot
    plt.figure(TASFIG)
    length = len(data[exper])
    dates = np.arange(length)
    if not exper=='PI-GWL-t6':
        plt.plot(dates, data[exper] - data['PI-GWL-t6'][:length],
                color=COLORS[exper],
                alpha=0.4,
                )
        plt.plot(dates[15:-14], moving_average(data[exper] - data['PI-GWL-t6'][:length], 30),
                color=COLORS[exper],
                label=LABELS[exper],
                )

plt.figure(TASFIG)
plt.legend(frameon=False)
plt.ylabel('$\Delta$TAS ($^{\circ}$C)')
plt.xlabel('Year')
plt.savefig('plots/tas_partial_forestation.png', dpi=200)


def plot_map(data:np.ndarray, title:str)->None:
    """Plot a global map pf surface air temperature.
    """
    ax = plt.axes(projection=ccrs.Robinson())
    discrete_bins = mpl.colors.BoundaryNorm(boundaries=np.arange(-3.25, 3.5, 0.5), ncolors=256)
    shading = plt.pcolormesh(lons, lats, data,
            cmap='bwr',
            edgecolors='face',
            linewidth=0.2,
            norm=discrete_bins,
            transform=ccrs.PlateCarree(),
            )
    plt.colorbar(shading,
            label='$\Delta$TAS $^{\circ}$C',
            orientation='horizontal',
            pad=0.05,
            ticks=np.arange(-3, 3.5, 0.5),
            )
    ax.coastlines()
    plt.title(title)


def plot_australia(data:np.ndarray, title:str)->None:
    """Plot a map of Asutralia temperature anomalies.
    """
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_extent([110, 160, -45, -10], crs=ccrs.PlateCarree())
    discrete_bins = mpl.colors.BoundaryNorm(boundaries=np.arange(-2.1, 2.2, 0.2), ncolors=256)
    shading = ax.pcolormesh(lons, lats, data,
            cmap='bwr',
            edgecolors='face',
            linewidth=0.2,
            norm=discrete_bins,
            transform=ccrs.PlateCarree(),
            )
    cbar = plt.colorbar(shading,
            label='TAS ($^{\circ}$C)',
            orientation='horizontal',
            pad=0.05,
            ticks=[-2, -1.5, -1, -.5, .5, 1, 1.5, 2],
            )
    ax.coastlines()
    plt.title(title)


plt.figure(2)
plot_map(maps['GWL-10pct-B2030'] - maps['PI-GWL-t6'], 'Forestation on 10% of crops')
plt.savefig('plots/tas_10pct_map.png', dpi=200)
plt.figure(3)
plot_map(maps['GWL-25pct-B2030'] - maps['PI-GWL-t6'], 'Forestation on 25% of crops')
plt.savefig('plots/tas_25pct_map.png', dpi=200)
plt.figure(4)
plot_map(maps['GWL-50pct-B2030'] - maps['PI-GWL-t6'], 'Forestation on 50% of crops')
plt.savefig('plots/tas_50pct_map.png', dpi=200)


plt.figure(5)
plot_australia(maps['GWL-10pct-B2030'] - maps['PI-GWL-t6'], 'Forestation on 10% of crops')
plt.savefig('plots/tas_10pct_australia.png', dpi=200)
plt.figure(6)
plot_australia(maps['GWL-25pct-B2030'] - maps['PI-GWL-t6'], 'Forestation on 25% of crops')
plt.savefig('plots/tas_25pct_australia.png', dpi=200)
plt.figure(7)
plot_australia(maps['GWL-50pct-B2030'] - maps['PI-GWL-t6'], 'Forestation on 50% of crops')
plt.savefig('plots/tas_50pct_australia.png', dpi=200)

plt.show()
