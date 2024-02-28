#!/usr/bin/env python3
# Plot the carbon pools from the single pft no-crops forestation experiments.
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
import ipdb


EXPERIMENTS = [
        'PI-GWL-t6',
        'GWL-NoCrops-B2030',
        'GWL-50pct-B2030',
        'GWL-25pct-B2030',
        'GWL-10pct-B2030',
        ]
ARCHIVE_DIR = '/g/data/p66/tfl561/archive_data'
WOODFIG = 1
LAST30 = str(-30*12)
TILE_FRAC_CODE = 'fld_s03i317'
RAW_CMIP_DIR = '/g/data/p73/archive/non-CMIP/ACCESS-ESM1-5'
RAW_NOCROP_DIR = '/g/data/p66/tfl561/ACCESS-ESM'
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


def plot_australia(data:np.ndarray, title:str)->None:
    """Plot a map of Asutralia temperature anomalies.
    """
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_extent([110, 160, -45, -10], crs=ccrs.PlateCarree())
    discrete_bins = mpl.colors.BoundaryNorm(
            boundaries=np.arange(-0.65, 0.65+0.1, 0.1),
            ncolors=256
            )
    shading = ax.pcolormesh(lons, lats, data,
            cmap='bwr',
            edgecolors='face',
            linewidth=0.2,
            norm=discrete_bins,
            transform=ccrs.PlateCarree(),
            )
    cbar = plt.colorbar(shading,
            label='$\Delta$cLand [Pg(C)]',
            orientation='horizontal',
            pad=0.05,
            ticks=np.arange(-0.6, 0.6+0.1, 0.1),
            )
    ax.coastlines()
    plt.title(title)


def plot_map(data:np.ndarray, title:str)->None:
    """Plot a global map pf surface air temperature.
    """
    ax = plt.axes(projection=ccrs.Robinson())
    discrete_bins = mpl.colors.BoundaryNorm(
            boundaries=np.arange(-0.65, 0.65+0.1, 0.1),
            ncolors=256
            )
    shading = plt.pcolormesh(lons, lats, data,
            cmap='bwr',
            norm=discrete_bins,
            transform=ccrs.PlateCarree(),
            linewidth=0.2,
            edgecolors='face',
            )
    ax.coastlines()
    plt.colorbar(shading,
            label='$\Delta$cLand [Pg(C)]',
            ticks=np.arange(-0.6, 0.6+0.1, 0.1),
            orientation='horizontal',
            pad=0.05,
            )
    plt.title(title)


example_file = f'/g/data/p66/tfl561/archive_data/GWL-EGNL-B2030/cLand_GWL-EGNL-B2030_0500-0602.nc'
ncfile = nc.Dataset(example_file)
lats = ncfile.variables['lat'][:]
lons = ncfile.variables['lon'][:]
ncfile.close()

plt.figure(WOODFIG)
data = dict()
maps = dict()
for exper in EXPERIMENTS:
    fractions_file = glob.glob(f'{ARCHIVE_DIR}/{exper}/frac_*.nc')[0]


    @cdod.cdo_mul(input2=fractions_file) # This experiment's cover fractions
    @cdod.cdo_mul(input2='data/cell_area.nc') # -> g(C)
    @cdod.cdo_fldsum
    @cdod.cdo_vertsum()
    @cdod.cdo_divc('1e15') # -> Pg(C)
    @cdod.cdo_yearmonmean
    def cdo_load(var:str, input:str)->np.ma.MaskedArray:
        """Load global sum of a variable using CDO.
        """
        return cdo.copy(
                input=input,
                options='-L',
                returnCdf=True,
                ).variables[var][:].squeeze()


    @cdod.cdo_mul(input2=fractions_file) # This experiment's cover fractions
    @cdod.cdo_mul(input2='data/cell_area.nc') # -> g(C)
    @cdod.cdo_vertsum()
    @cdod.cdo_seltimestep(f'{LAST30}/-1') # last 20 years
    @cdod.cdo_timmean # Temporal average
    @cdod.cdo_divc('1e15') # -> Pg(C)
    def cdo_load_last30(var:str, input:str)->np.ma.MaskedArray:
        """Load last 30 year avearge of a variable using CDO.
        """
        return cdo.copy(
                input=input,
                options='-L',
                returnCdf=True,
                ).variables[var][:].squeeze()


    # Load the cLand
    print("Loading cLand", exper)
    files = sorted(glob.glob(f'{ARCHIVE_DIR}/{exper}/cLand_*'))

    # Load data
    data[exper] = cdo_load(input=files[0], var='cLand')
    maps[exper] = cdo_load_last30(input=files[0], var='cLand')

    # Plot
    length = len(data[exper])
    dates = np.arange(length)#/12.0
    if not exper=='PI-GWL-t6':
        plt.plot(dates, data[exper] - data['PI-GWL-t6'][:length],
                color=COLORS[exper],
                label=LABELS[exper],
                )

plt.legend(frameon=False)
plt.ylabel('$\Delta$cLand [Pg(C)]')
plt.xlabel('Year')
plt.savefig('plots/cLand_partial_forestation.png', dpi=200)

plt.figure()
plot_map(maps['GWL-10pct-B2030'] - maps['PI-GWL-t6'], title='10%')
plt.savefig('cLand_10pct_last30.png', dpi=200)
plt.figure()
plot_map(maps['GWL-25pct-B2030'] - maps['PI-GWL-t6'], title='25%')
plt.savefig('cLand_25pct_last30.png', dpi=200)
plt.figure()
plot_map(maps['GWL-50pct-B2030'] - maps['PI-GWL-t6'], title='50%')
plt.savefig('cLand_50pct_last30.png', dpi=200)

plt.figure()
plot_australia(maps['GWL-10pct-B2030'] - maps['PI-GWL-t6'], title='10%')
plt.savefig('cLand_10pct_australia.png', dpi=200)
plt.figure()
plot_australia(maps['GWL-25pct-B2030'] - maps['PI-GWL-t6'], title='25%')
plt.savefig('cLand_25pct_australia.png', dpi=200)
plt.figure()
plot_australia(maps['GWL-50pct-B2030'] - maps['PI-GWL-t6'], title='50%')
plt.savefig('cLand_50pct_australia.png', dpi=200)

plt.show()
