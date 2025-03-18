#!/usr/bin/env python3
"""Plot figures for the GWL-NoCrAu-B2035 experiment.
This eperiment has only forestation applied to croplands in Australia.
It is not expected to have a large impact on global atmoshperic CO2 or temperatures.
However it is expected to have an impact on local temperatures.
These results are intended to be in anticipation of reviewer comments and are likely to only be
present in the supplmentary material.
"""
import cartopy.crs as ccrs
import cdo_decorators as cdod
import ipdb
import matplotlib as mpl
import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np

from cdo import Cdo; cdo = Cdo(); cdo.debug = False

ARCHIVE_DIR = '/g/data/p66/tfl561/archive_data'
ANOTHER_EXPER = f'{ARCHIVE_DIR}/GWL-NoCrops-B2030/tas_GWL-NoCrops-B2030_0500-0700.nc'
DPI = 200
EXAMPLE_FLIE = f'{ARCHIVE_DIR}/GWL-NoCrops-B2030/cLeaf_GWL-NoCrops-B2030_0500-0700.nc'
EXPER_SIM = f'{ARCHIVE_DIR}/GWL-NoCrAu-B2030/tas_GWL-NoCrAu-B2030_0500-0700.nc'
LAST100 = str(-100*12)
READ_ONLY = 'r'
REF_SIM = f'{ARCHIVE_DIR}/PI-GWL-t6/tas_PI-GWL-t6_0500-0700.nc'


@cdod.cdo_sellonlatbox('110','155','-45','-10') # Australia region.
@cdod.cdo_yearmonmean
@cdod.cdo_fldmean(weights='TRUE') # Australia mean.
def load_australia_mean(var:str, input:str)->np.ma.MaskedArray:
    """Load Australia data and spatially average.
    """
    ncfile = cdo.copy(input=input, returnCdf=True, options='-L')
    return ncfile.variables[var][:].squeeze()


@cdod.cdo_seltimestep('1') # Select the first time step.
@cdod.cdo_sellevel('9') # Select crops vegetation type.
@cdod.cdo_gtc('0') # Select cropped regions.
@cdod.cdo_setctomiss('0') # Set non-forestation regions to missing.
def make_forestation_mask(input:str)->None:
    """Make a crops region mask.
    """
    cdo.copy(input=input, output='forestation_mask.nc', options='-L')


@cdod.cdo_mul(input2='forestation_mask.nc') # mask out non-forestation.
@cdod.cdo_sellonlatbox('110', '155', '-45', '-10') # Australia region.
@cdod.cdo_yearmonmean
@cdod.cdo_fldmean(weights='TRUE') # Australia mean.
def load_forestation_mean(var:str, input:str)->np.ndarray:
    """Load data in Australian forestation regions and spatially average.
    """
    ncfile = cdo.copy(input=input, returnCdf=True, options='-L')
    return ncfile.variables[var][:].squeeze()


@cdod.cdo_fldmean(weights='TRUE') # Global mean.
@cdod.cdo_yearmonmean
def load_global_mean(var, input:str)->np.ma.MaskedArray:
    """Load global mean of a variable.
    """
    ncfile = cdo.copy(input=input, returnCdf=True, options='-L')
    return ncfile.variables[var][:].squeeze()


@cdod.cdo_seltimestep(f'{LAST100}/-1') # last 100 years
@cdod.cdo_timmean # Temporal average
def load_last100(var, input:str)->np.ma.MaskedArray:
    """Load last 100 year mean of a carbon pool variable.
    """
    ncfile = cdo.copy(input=input, returnCdf=True, options='-L')
    return ncfile.variables[var][:].squeeze()


def moving_average(x:np.ndarray, window:int)->np.ndarray:
    """Calculate a moving average using a window of size `window`.
    """
    return np.convolve(x, np.ones(window), mode='same')/window


def plot_australia(data:np.ndarray, title:str, save:str)->None:
    """Plot global data focussed on only Australia.
    """
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    australia_lonlat = [110,155,-45,-10]
    #ax.set_extent(australia_lonlat, crs=ccrs.PlateCarree())
    ax.set_extent([110,165,-55,-10], crs=ccrs.PlateCarree())
    discrete_bins = mpl.colors.BoundaryNorm(boundaries=np.arange(-2.1, 2.2, 0.2), ncolors=256)
    colors = ax.pcolormesh(lons, lats, data,
            cmap='seismic',
            norm=discrete_bins,
            transform=ccrs.PlateCarree(),
            )
    ax.coastlines()
    plt.colorbar(colors,
            label='$\Delta$ Surface air temperature ($^{\circ}$C)',
            orientation='horizontal',
            ticks=[-2,-1.5,-1,-.5,0,.5,1,1.5,2],
            pad=0.05,
            )
    plt.title(title)
    plt.savefig(f'plots/{save}', dpi=DPI)


def yearly_mean_from_monthly(data:np.ndarray)->np.ndarray:
    """Calculate a yearly mean on a numpy array of monthly data.
    The 0th dimension must be time and divisible by 12.
    """
    if np.mod(data.shape[0], 12)!=0:
        raise ValueError("Not enough months in 0th dimension.")
    toshape = list(data.shape)
    toshape.pop(0)
    toshape.insert(0, 12)
    toshape.insert(0, int(data.shape[0]/12))
    fraction_of_year = np.array([31,28,31,30,31,30,31,31,30,31,30,31])/365.0
    return np.average(data.reshape(toshape), axis=1, weights=fraction_of_year)


# Load lats and lons.
ncfile = nc.Dataset(EXAMPLE_FLIE, 'r')
lats = ncfile.variables['lat'][:]
lons = ncfile.variables['lon'][:]

# Plot the global mean temperature with the reference simulation.
# Demonstrate that there is little impact on the global mean temperature.
temperature_reference = load_global_mean('tas', input=REF_SIM)
temperature_another = load_global_mean('tas', input=ANOTHER_EXPER)
temperature_exper = load_global_mean('tas', input=EXPER_SIM)

plt.figure()
anomaly = temperature_exper - temperature_reference
length = len(anomaly)
dates = np.arange(length)
color = 'blue'
plt.plot(anomaly,
        color=color,
        label='Forestation in Australia only',
        alpha=0.4,
        )
plt.plot(moving_average(anomaly, window=50),
        color=color,
        )

color = 'green'
anomaly = temperature_another - temperature_reference
plt.plot(anomaly,
        color=color,
        label='Forestation at 2030',
        alpha=0.4,
        )
plt.plot(moving_average(anomaly, window=50),
        color=color,
        )
plt.hlines(0, dates[0], dates[-1], colors='k')
plt.xlim(left=0, right=200)
plt.legend()
plt.ylabel('Surface air temperature ($^{\circ}$C)')
plt.xlabel('Time (year)')
plt.title('Global mean temperature anomaly from forestation')
plt.savefig('plots/tas_GWL-NoCrAu-B2030_global_tseries.png', dpi=DPI)

# Plot the time series of Australia wide mean.
aus_temp_reference = load_australia_mean('tas', input=REF_SIM)
aus_temp_another = load_australia_mean('tas', input=ANOTHER_EXPER)
aus_temp_exper = load_australia_mean('tas', input=EXPER_SIM)

plt.figure()
anomaly = aus_temp_exper - aus_temp_reference
color = 'blue'
plt.plot(anomaly,
        color=color,
        label='Forestation in Australia only',
        alpha=0.4,
        )
plt.plot(moving_average(anomaly, window=50),
        color=color,
        )

anomaly = aus_temp_another - aus_temp_reference
color = 'green'
plt.plot(anomaly,
        color=color,
        label='Forestation at 2030',
        alpha=0.4,
        )
plt.plot(moving_average(anomaly, window=50),
        color=color,
        )
plt.hlines(0, dates[0], dates[-1], colors='k')
plt.xlim(left=0, right=dates[-1])
plt.legend()
plt.ylabel('Surface air temperature ($^{\circ}$C)')
plt.xlabel('Time (year)')
plt.title('Australia mean temperature anomaly from forestation')
plt.savefig('plots/tas_GWL-NoCrAu-B2030_australia_tseries.png', dpi=DPI)

# Plot the time series of afforested area mean.
# Create mask for forestation only regions.
make_forestation_mask(input=f'{ARCHIVE_DIR}/PI-GWL-t6/frac_PI-GWL-t6_0500-0700.nc')

aus_temp_reference = load_forestation_mean('tas', input=REF_SIM)
aus_temp_another = load_forestation_mean('tas', input=ANOTHER_EXPER)
aus_temp_exper = load_forestation_mean('tas', input=EXPER_SIM)

plt.figure()
anomaly = aus_temp_exper - aus_temp_reference
color = 'blue'
plt.plot(anomaly,
        color=color,
        label='Forestation in Australia only',
        alpha=0.4,
        )
plt.plot(moving_average(anomaly, window=50),
        color=color,
        )

anomaly = aus_temp_another - aus_temp_reference
color = 'green'
plt.plot(anomaly,
        color=color,
        label='Forestation at 2030',
        alpha=0.4,
        )
plt.plot(moving_average(anomaly, window=30),
        color=color,
        )
plt.hlines(0, dates[0], dates[-1], colors='k')
plt.xlim(left=0, right=dates[-1])
plt.legend()
plt.ylabel('Surface air temperature ($^{\circ}$C)')
plt.xlabel('Time (year)')
plt.title('Australia afforested area mean temperature anomaly')
plt.savefig('plots/tas_GWL-NoCrAu-B2030_australia_forests_tseries.png', dpi=DPI)

# plot the map of the end of the simulation.
temp_last100_ref = load_last100('tas', input=REF_SIM)
temp_last100_another = load_last100('tas', input=ANOTHER_EXPER)
temp_last100_exper = load_last100('tas', input=EXPER_SIM)

plot_australia(temp_last100_exper - temp_last100_ref,
        title='',
        save='tas_aus_only.png')
plt.savefig('plots/tas_GWL-NoCrAu-B2030_australia_map.png', dpi=DPI)
plt.show()
