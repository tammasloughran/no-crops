#!/usr/bin/env python3
# Plot the carbon pools from the single pft no-crops forestation experiments.
import glob
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np
from cartopy import crs as ccrs

[sys.path.append(i) for i in ['.', '..']]
from cdo import Cdo

from contrib.tools import global_mean

cdo = Cdo()
cdo.debug = False
import cdo_decorators as cdod
import ipdb
import scipy.stats as stats

EXPERIMENTS = [
        'PI-GWL-t6',
        'GWL-NoCrops-B2030',
        'GWL-50pct-B2030',
        'GWL-25pct-B2030',
        'GWL-10pct-B2030',
        'GWL-50pc-B2030-02',
        'GWL-25pc-B2030-02',
        'GWL-10pc-B2030-02',
        ]
ARCHIVE_DIR = '/g/data/p66/tfl561/archive_data'
WOODFIG = 1
LAST100 = str(-100*12)
TILE_FRAC_CODE = 'fld_s03i317'
RAW_CMIP_DIR = '/g/data/p73/archive/non-CMIP/ACCESS-ESM1-5'
RAW_NOCROP_DIR = '/g/data/p66/tfl561/ACCESS-ESM'
LABELS = {
        'PI-GWL-t6':'Global warming level',
        'GWL-NoCrops-B2030':'100%',
        'GWL-10pct-B2030':'10%',
        'GWL-10pc-B2030-02':'10%',
        'GWL-25pct-B2030':'25%',
        'GWL-25pc-B2030-02':'25%',
        'GWL-50pct-B2030':'50%',
        'GWL-50pc-B2030-02':'50%',
        }
COLORS = {
        'PI-GWL-t6':'black',
        'GWL-NoCrops-B2030':'blue',
        'GWL-10pct-B2030':'lightgreen',
        'GWL-10pc-B2030-02':'lightgreen',
        'GWL-25pct-B2030':'forestgreen',
        'GWL-25pc-B2030-02':'forestgreen',
        'GWL-50pct-B2030':'darkgreen',
        'GWL-50pc-B2030-02':'darkgreen',
        }


def regression_point(reg, x:float):
    """Return a point on a regression line.
    """
    return reg.intercept + reg.slope*x


example_file = f'/g/data/p66/tfl561/archive_data/GWL-EGNL-B2030/cLand_GWL-EGNL-B2030_0500-0700.nc'
ncfile = nc.Dataset(example_file)
lats = ncfile.variables['lat'][:]
lons = ncfile.variables['lon'][:]
ncfile.close()

plt.figure(WOODFIG)
data = dict()
for exper in EXPERIMENTS:
    fractions_file = glob.glob(f'{ARCHIVE_DIR}/{exper}/frac_*.nc')[0]


    @cdod.cdo_mul(input2=fractions_file) # This experiment's cover fractions
    @cdod.cdo_mul(input2='data/cell_area.nc') # -> g(C)
    @cdod.cdo_sellonlatbox('110','155','-45','-10') # Australia region.
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


    # Load the cLand
    print("Loading cLand", exper)
    files = sorted(glob.glob(f'{ARCHIVE_DIR}/{exper}/cLand_*'))

    # Load data
    data[exper] = cdo_load(input=files[0], var='cLand')

    # Plot
    length = len(data[exper])
    dates = np.arange(length)#/12.0
    if not exper=='PI-GWL-t6':
        plt.plot(dates, data[exper] - data['PI-GWL-t6'][:length],
                color=COLORS[exper],
                #label=LABELS[exper],
                alpha=0.2,
                )

for exper in ['50','25','10']:
    e1 = f'GWL-{exper}pct-B2030'
    e2 = f'GWL-{exper}pc-B2030-02'
    length = len(data[e2])
    ens_mean = (data[e1][:length] + data[e2])/2
    dates = np.arange(length)
    plt.plot(dates, ens_mean - data['PI-GWL-t6'][:length], color=COLORS[e1], label=LABELS[e1])
plt.plot(
        np.arange(len(data['PI-GWL-t6'])),
        data['GWL-NoCrops-B2030'] - data['PI-GWL-t6'],
        color='blue',
        label='100%',
        )

plt.legend(frameon=False)
plt.ylabel('$\Delta$cLand [Pg(C)]')
plt.xlabel('Year')
plt.hlines(0, xmin=dates[0], xmax=dates[-1], color='black')
plt.xlim(left=dates[0], right=dates[-1])
plt.title("Australia sum total land carbon")
plt.savefig('plots/cLand_partial_forestation_australia_tseries.png', dpi=200)

# Plot the ensemble mean.
plt.figure()
for exper in ['50','25','10']:
    e1 = f'GWL-{exper}pct-B2030'
    e2 = f'GWL-{exper}pc-B2030-02'
    length = len(data[e2])
    ens_mean = (data[e1][:length] + data[e2])/2
    dates = np.arange(length)
    plt.plot(
            dates,
            ens_mean - data['PI-GWL-t6'][:length],
            color=COLORS[e1],
            label=LABELS[e1],
            )
plt.plot(
        np.arange(len(data['PI-GWL-t6'])),
        data['GWL-NoCrops-B2030'] - data['PI-GWL-t6'],
        color='blue',
        label='100%',
        )
plt.hlines(0, xmin=dates[0], xmax=dates[-1], color='black')
plt.legend(frameon=False)
plt.ylabel('$\Delta$cLand [Pg(C)]')
plt.xlabel('Year')
plt.xlim(left=dates[0], right=dates[-1])
plt.title("Australia sum total land carbon")

plt.show()

# Do a regression of Area vs C uptake
c_uptake = []
area = []
for x in EXPERIMENTS[1:]:
    area.append(int(LABELS[x][:-1]))
    c_uptake.append(np.mean(data[x][-100:] - data['PI-GWL-t6'][-100:]))
plt.figure()
plt.scatter(area, c_uptake)
regression = stats.linregress(area, c_uptake)
plt.plot([0,100], [regression_point(regression, 0), regression_point(regression, 100)])
plt.xlabel('Forestation on croplands (%)')
plt.ylabel('$\Delta$ CLand [Pg (C)]')
plt.show()

