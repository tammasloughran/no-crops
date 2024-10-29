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
        'GWL-NoCr-B2030-02',
        'GWL-50pct-B2030',
        'GWL-25pct-B2030',
        'GWL-10pct-B2030',
        'GWL-50pc-B2030-02',
        'GWL-25pc-B2030-02',
        'GWL-10pc-B2030-02',
        ]
ARCHIVE_DIR = '/g/data/p66/tfl561/archive_data'
PRFIG = 1
MAPFIG = 2
SEC_PER_DAY = str(60*60*24)
TILE_FRAC_CODE = 'fld_s03i317'
RAW_CMIP_DIR = '/g/data/p73/archive/non-CMIP/ACCESS-ESM1-5'
RAW_NOCROP_DIR = '/g/data/p66/tfl561/ACCESS-ESM'
EXAMPLE = '/g/data/p66/tfl561/archive_data/GWL-NoCrops-B2040/pr_GWL-NoCrops-B2040_0500-0700.nc'
LABELS = {
        'PI-GWL-t6':'Global warming level',
        'GWL-NoCrops-B2030':'100%',
        'GWL-NoCr-B2030-02':'100%',
        'GWL-10pct-B2030':'10%',
        'GWL-25pct-B2030':'25%',
        'GWL-50pct-B2030':'50%',
        'GWL-10pc-B2030-02':'10%',
        'GWL-25pc-B2030-02':'25%',
        'GWL-50pc-B2030-02':'50%',
        }
COLORS = {
        'PI-GWL-t6':'black',
        'GWL-NoCrops-B2030':'blue',
        'GWL-NoCr-B2030-02':'blue',
        'GWL-10pct-B2030':'lightgreen',
        'GWL-25pct-B2030':'forestgreen',
        'GWL-50pct-B2030':'darkgreen',
        'GWL-10pc-B2030-02':'lightgreen',
        'GWL-25pc-B2030-02':'forestgreen',
        'GWL-50pc-B2030-02':'darkgreen',
        }


@cdod.cdo_sellonlatbox('110','155','-45','-10') # Australia region.
@cdod.cdo_selseason('JJA')
@cdod.cdo_fldmean(weights='TRUE')
@cdod.cdo_mulc(SEC_PER_DAY)
@cdod.cdo_yearmonmean
def cdo_load(var:str, input:str)->np.ma.MaskedArray:
    """Load global mean of a variable using CDO."""
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
plt.figure(PRFIG)

data = {}
for exper in EXPERIMENTS:
    # Load the pr data
    print("Loading Precip", exper)
    files = sorted(glob.glob(f'{ARCHIVE_DIR}/{exper}/pr_{exper}_*'))

    # Load data
    data[exper] = cdo_load(input=files[0], var='pr')

ensmean = {}
ensmean['GWL-NoCrops-B2030'] = (data['GWL-NoCrops-B2030'][:-10] + data['GWL-NoCr-B2030-02'])/2
ensmean['GWL-50pct-B2030'] = (data['GWL-50pct-B2030'][:-9] + data['GWL-50pc-B2030-02'])/2
ensmean['GWL-25pct-B2030'] = (data['GWL-25pct-B2030'][:-10] + data['GWL-25pc-B2030-02'])/2
ensmean['GWL-10pct-B2030'] = (data['GWL-10pct-B2030'][:-10] + data['GWL-10pc-B2030-02'])/2

# Plot
for exper in [
        'GWL-NoCrops-B2030',
        'GWL-50pct-B2030',
        'GWL-25pct-B2030',
        'GWL-10pct-B2030',
        ]:
    length = len(ensmean[exper])
    dates = np.arange(400, 400 + length)
    plt.plot(dates, ensmean[exper] - data['PI-GWL-t6'][:length],
            color=COLORS[exper],
            alpha=0.4,
            )
    plt.plot(dates[15:-14], moving_average(ensmean[exper] - data['PI-GWL-t6'][:length], 30),
            color=COLORS[exper],
            label=LABELS[exper],
            )
plt.hlines(0, dates[0], dates[-1], colors='k')
plt.legend(frameon=False)
plt.xlim(left=400, right=400+length)
plt.ylabel('$\Delta$Precip. (mm day$^{-1}$)')
plt.xlabel('Year')
plt.savefig('plots/pr_partial_forestation_australia.png', dpi=200)

plt.show()
