#!/usr/bin/env python3
# Plot the tas from the single pft no-crops forestation experiments.
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from cartopy import crs as ccrs
import netCDF4 as nc
import glob
import sys
[sys.path.append(i) for i in ['.', '..']]
from contrib.tools import global_mean
from cdo import Cdo
cdo = Cdo()
cdo.debug = False
import cdo_decorators as cdod


EXPERIMENTS = [
        'PI-GWL-t6',
        'GWL-NoCrops-B2030',
        'GWL-EGNL-B2030',
        'GWL-EGBL-B2030',
        'GWL-DCBL-B2030',
        ]
ARCHIVE_DIR = '/g/data/p66/tfl561/archive_data'
TASFIG = 1
MAPFIG = 2
TILE_FRAC_CODE = 'fld_s03i317'
RAW_CMIP_DIR = '/g/data/p73/archive/non-CMIP/ACCESS-ESM1-5'
RAW_NOCROP_DIR = '/g/data/p66/tfl561/ACCESS-ESM'
EXAMPLE = '/g/data/p66/tfl561/archive_data/GWL-NoCrops-B2040/tas_GWL-NoCrops-B2040_0500-0700.nc'
LABELS = {
        'PI-GWL-t6':'Global warming level',
        'GWL-NoCrops-B2030':'Mixed forest',
        'GWL-EGNL-B2030':'Evergreen needle leaf',
        'GWL-EGBL-B2030':'Evergreen broad leaf',
        'GWL-DCBL-B2030':'Deciduous broad leaf',
        }
COLORS = {
        'PI-GWL-t6':'black',
        'GWL-NoCrops-B2030':'blue',
        'GWL-EGNL-B2030':'lightgreen',
        'GWL-EGBL-B2030':'darkgreen',
        'GWL-DCBL-B2030':'orange',
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
        plt.plot(dates[15:-14], moving_average(data[exper] - data['PI-GWL-t6'][:length], 50),
                color=COLORS[exper],
                label=LABELS[exper],
                )


plt.figure(TASFIG)
plt.legend(frameon=False)
plt.ylabel('$\Delta$TAS ($^{\circ}$C)')
plt.xlabel('Year')
plt.savefig('plots/tas_single_pft_forestation.png', dpi=200)


def plot_map(data:np.ndarray, title:str)->None:
    """Plot a global map pf surface air temperature.
    """
    ax = plt.axes(projection=ccrs.Robinson())
    discrete_bins = mpl.colors.BoundaryNorm(boundaries=np.arange(-3.25, 3.5, 0.5), ncolors=256)
    shading = plt.pcolormesh(lons, lats, data,
            cmap='bwr',
            norm=discrete_bins,
            transform=ccrs.PlateCarree(),
            linewidth=0.2,
            edgecolors='face',
            )
    ax.coastlines()
    plt.colorbar(shading,
            label='$\Delta$TAS $^{\circ}$C',
            ticks=np.arange(-3, 3.5, 0.5),
            orientation='horizontal',
            pad=0.05,
            )
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
plot_map(maps['GWL-EGBL-B2030'] - maps['PI-GWL-t6'], 'Evergreen broadleaf')
plt.savefig('plots/tas_EGBL_map.png', dpi=200)
plt.figure(3)
plot_map(maps['GWL-EGNL-B2030'] - maps['PI-GWL-t6'], 'Evergreen needleleaf')
plt.savefig('plots/tas_EGNL_map.png', dpi=200)
plt.figure(4)
plot_map(maps['GWL-DCBL-B2030'] - maps['PI-GWL-t6'], 'Deciduous broadleaf')
plt.savefig('plots/tas_DCBL_map.png', dpi=200)

plt.figure(5)
plot_australia(maps['GWL-EGBL-B2030'] - maps['PI-GWL-t6'], 'Evergreen broadleaf')
plt.savefig('plots/tas_EGBL_australia.png', dpi=200)
plt.figure(6)
plot_australia(maps['GWL-EGNL-B2030'] - maps['PI-GWL-t6'], 'Evergreen needleleaf')
plt.savefig('plots/tas_EGNL_australia.png', dpi=200)
plt.figure(7)
plot_australia(maps['GWL-DCBL-B2030'] - maps['PI-GWL-t6'], 'Deciduous broadleaf')
plt.savefig('plots/tas_DCBL_australia.png', dpi=200)

plt.show()
