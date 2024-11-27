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
LAST30 = str(-30*12)
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
VARIABLE = 'pr'
SEC_PER_DAY = str(60*60*24)


def moving_average(x:np.ndarray, window:int)->np.ndarray:
    """Calculate a moving average using a window of size `window`.
    """
    return np.convolve(x, np.ones(window), mode='same')/window


@cdod.cdo_sellonlatbox('110','155','-45','-10') # Australia region.
@cdod.cdo_selseason('JJA')
@cdod.cdo_fldmean(weights='TRUE') # Australia mean.
@cdod.cdo_mulc(SEC_PER_DAY)
@cdod.cdo_yearmonmean
def cdo_load(var:str, input:str)->np.ma.MaskedArray:
    """Load global sum of a variable using CDO.
    """
    return cdo.copy(
            input=input,
            options='-L',
            returnCdf=True,
            ).variables[var][:].squeeze()


load_npy = True

example_file = f'/g/data/p66/tfl561/archive_data/GWL-EGNL-B2030/pr_GWL-EGNL-B2030_0500-0700.nc'
ncfile = nc.Dataset(example_file)
lats = ncfile.variables['lat'][:]
lons = ncfile.variables['lon'][:]
ncfile.close()

plt.figure(WOODFIG)
data = dict()
maps = dict()
for exper in EXPERIMENTS:
    # Load the data
    print(f"Loading {VARIABLE}", exper)
    files = sorted(glob.glob(f'{ARCHIVE_DIR}/{exper}/{VARIABLE}_*'))

    # Load data
    if load_npy:
        data[exper] = np.load(f'data/pr_aus_{exper}.npy')
    else:
        data[exper] = cdo_load(input=files[0], var=VARIABLE)
        np.save(f'data/pr_aus_{exper}.npy', data[exper].data)

for exper in [
        'GWL-EGNL-B2030',
        'GWL-EGBL-B2030',
        'GWL-DCBL-B2030',
        ]:
    # Plot
    length = len(data[exper])
    dates = np.arange(length)#/12.0
    plt.plot(dates, data[exper] - data['PI-GWL-t6'][:length],
            color=COLORS[exper],
            alpha=0.3,
            )
    length = len(data[exper+'-02'])
    dates = np.arange(length)#/12.0
    plt.plot(dates, data[exper+'-02'] - data['PI-GWL-t6'][:length],
            color=COLORS[exper],
            alpha=0.3,
            )
    length = len(data[exper+'-02'])
    ens_mean = (data[exper][:length] + data[exper+'-02'])/2.0
    plt.plot(dates[:length], moving_average(ens_mean - data['PI-GWL-t6'][:length], window=30),
            color=COLORS[exper],
            label=LABELS[exper],
            )
length = len(data['GWL-NoCrops-B2030'])
dates = np.arange(length)#/12.0
plt.plot(
        dates,
        moving_average(data['GWL-NoCrops-B2030'] - data['PI-GWL-t6'][:length], window=30),
        color=COLORS['GWL-NoCrops-B2030'],
        label=LABELS['GWL-NoCrops-B2030'],
        )
plt.legend(frameon=False)
plt.hlines(0, dates[0], dates[-1], colors='k')
plt.xlim(left=dates[0], right=dates[-1])
plt.ylabel(f'Precip. (mm/day)')
plt.xlabel('Year')
plt.savefig(f'plots/{VARIABLE}_single_pft_forestation_australia_tseries.png', dpi=200)

plt.show()
