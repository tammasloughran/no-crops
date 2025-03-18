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


EXPERIMENTS = [
        'PI-GWL-t6',
        'GWL-NoCrops-B2030',
        'GWL-EGNL-B2030',
        'GWL-EGBL-B2030',
        'GWL-DCBL-B2030',
        ]
ARCHIVE_DIR = '/g/data/p66/tfl561/archive_data'
WOODFIG = 1
LAST100 = str(-100*12)
TILE_FRAC_CODE = 'fld_s03i317'
RAW_CMIP_DIR = '/g/data/p73/archive/non-CMIP/ACCESS-ESM1-5'
RAW_NOCROP_DIR = '/g/data/p66/tfl561/ACCESS-ESM'
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
VARIABLE = 'cLand'


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


example_file = f'/g/data/p66/tfl561/archive_data/GWL-EGNL-B2030/cLand_GWL-EGNL-B2030_0500-0700.nc'
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
    @cdod.cdo_seltimestep(f'{LAST100}/-1') # last 100 years
    @cdod.cdo_timmean # Temporal average
    @cdod.cdo_divc('1e15') # -> Pg(C)
    def cdo_load_last100(var:str, input:str)->np.ma.MaskedArray:
        """Load last 100 year avearge of a variable using CDO.
        """
        return cdo.copy(
                input=input,
                options='-L',
                returnCdf=True,
                ).variables[var][:].squeeze()

    # Load the data
    print(f"Loading {VARIABLE}", exper)
    files = sorted(glob.glob(f'{ARCHIVE_DIR}/{exper}/{VARIABLE}_*'))

    # Load data
    data[exper] = cdo_load(input=files[0], var=VARIABLE)
    maps[exper] = cdo_load_last100(input=files[0], var=VARIABLE)

    # Plot
    length = len(data[exper])
    dates = np.arange(length)#/12.0
    if not exper=='PI-GWL-t6':
        plt.plot(dates, data[exper] - data['PI-GWL-t6'][:length],
                color=COLORS[exper],
                label=LABELS[exper],
                )

plt.legend(frameon=False)
plt.ylabel(f'{VARIABLE} [Pg(C)]')
plt.xlabel('Year')
plt.savefig(f'plots/{VARIABLE}_single_pft_forestation.png', dpi=200)

plt.figure()
plot_map(maps['GWL-EGNL-B2030'] - maps['PI-GWL-t6'], title='Evergreen needle leaf')
plt.savefig('plots/cLand_EGNL_last100.png', dpi=200)
plt.figure()
plot_map(maps['GWL-EGBL-B2030'] - maps['PI-GWL-t6'], title='Evergreen broad leaf')
plt.savefig('plots/cLand_EGBL_last100.png', dpi=200)
plt.figure()
plot_map(maps['GWL-DCBL-B2030'] - maps['PI-GWL-t6'], title='Deciduous broad leaf')
plt.savefig('plots/cLand_DCBL_last100.png', dpi=200)

plt.figure()
plot_australia(maps['GWL-EGNL-B2030'] - maps['PI-GWL-t6'], title='Evergreen needle leaf')
plt.savefig('plots/cLand_EGNL_australia.png', dpi=200)
plt.figure()
plot_australia(maps['GWL-EGBL-B2030'] - maps['PI-GWL-t6'], title='Evergreen broad leaf')
plt.savefig('plots/cLand_EGBL_australia.png', dpi=200)
plt.figure()
plot_australia(maps['GWL-DCBL-B2030'] - maps['PI-GWL-t6'], title='Deciduous broad leaf')
plt.savefig('plots/cLand_DCBL_australia.png', dpi=200)

plt.show()
