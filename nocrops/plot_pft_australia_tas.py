#!/usr/bin/env python3
# Plot the carbon pools from the single pft no-crops forestation experiments for australia.
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
        'GWL-EGNL-B2030',
        'GWL-EGNL-B2030-02',
        'GWL-EGBL-B2030',
        'GWL-EGBL-B2030-02',
        'GWL-DCBL-B2030',
        'GWL-DCBL-B2030-02',
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
        'GWL-EGNL-B2030-02':'',
        'GWL-EGBL-B2030':'Evergreen broad leaf',
        'GWL-EGBL-B2030-02':'',
        'GWL-DCBL-B2030':'Deciduous broad leaf',
        'GWL-DCBL-B2030-02':'',
        }
COLORS = {
        'PI-GWL-t6':'black',
        'GWL-NoCrops-B2030':'blue',
        'GWL-EGNL-B2030':'lightgreen',
        'GWL-EGNL-B2030-02':'lightgreen',
        'GWL-EGBL-B2030':'darkgreen',
        'GWL-EGBL-B2030-02':'darkgreen',
        'GWL-DCBL-B2030':'orange',
        'GWL-DCBL-B2030-02':'orange',
        }
VARIABLE = 'tas'
SEC_PER_DAY = str(60*60*24)


@cdod.cdo_sellonlatbox('110','155','-45','-10') # Australia region.
@cdod.cdo_fldmean(weights='TRUE') # Australia mean.
@cdod.cdo_yearmonmean
def cdo_load(var:str, input:str)->np.ma.MaskedArray:
    """Load global sum of a variable using CDO.
    """
    return cdo.copy(
            input=input,
            options='-L',
            returnCdf=True,
            ).variables[var][:].squeeze()


def moving_average(x:np.ndarray, window:int)->np.ndarray:
    """Calculate a moving average using a window of size `window`.
    """
    return np.convolve(x, np.ones(window), mode='same')/window


example_file = f'/g/data/p66/tfl561/archive_data/GWL-EGNL-B2030/{VARIABLE}_GWL-EGNL-B2030_0500-0700.nc'
ncfile = nc.Dataset(example_file)
lats = ncfile.variables['lat'][:]
lons = ncfile.variables['lon'][:]
ncfile.close()

plt.figure(WOODFIG)
data = {}
maps = {}
for exper in EXPERIMENTS:
    # Load the data
    print(f"Loading {VARIABLE}", exper)
    files = sorted(glob.glob(f'{ARCHIVE_DIR}/{exper}/{VARIABLE}_*'))

    # Load data
    data[exper] = cdo_load(input=files[0], var=VARIABLE)

    # Plot
    length = len(data[exper])
    dates = np.arange(length)#/12.0
    if not exper=='PI-GWL-t6':
        plt.plot(dates, data[exper] - data['PI-GWL-t6'][:length],
                color=COLORS[exper],
                alpha=0.4,
                )

for exper in [
        'GWL-EGNL-B2030',
        'GWL-EGBL-B2030',
        'GWL-DCBL-B2030',
        ]:
    plt.plot(dates[:190], moving_average((data[exper][:190]+data[exper+'-02'][:190])/2 - data['PI-GWL-t6'][:190], 50),
            color=COLORS[exper],
            label=LABELS[exper],
            )

plt.legend(frameon=False)
plt.hlines(0, dates[0], dates[-1], colors='k')
plt.xlim(left=dates[0], right=dates[-1])
plt.ylim(top=1)
plt.ylabel(f'{VARIABLE} ($^o$C)')
plt.xlabel('Year')
plt.savefig(f'plots/{VARIABLE}_single_pft_forestation_australia_tseries.png', dpi=200)

plt.show()
