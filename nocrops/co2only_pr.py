#!/usr/bin/env python3
"""Plot the temperature and precipitation difference between the 2030 forestation experiment
and the CO2 only experiment. This would represent the biogeochemical climate response to
forestation vs the biogeophysical response.

Experiments:
    PI-GWL-t6 - reference experiment of 2030 stabilized climate.
    GWL-CO2only-B2030 - biogeochemical climate effects of forestation only.
    GWL-NoCrops-B2030 - total climate effects of forestation.
"""
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import matplotlib as mpl
from cartopy import crs as ccrs
from cdo import Cdo
import cdo_decorators as cdod

cdo = Cdo()
cdo.debug = True

# Constants
ARCHIVE = '/g/data/p66/tfl561/archive_data'

TAS_REF_FILENAME = f'/g/data/p66/tfl561/archive_data/PI-GWL-t6/pr_PI-GWL-t6_0500-0700.nc'
TAS_BGC_FILENAME = f'{ARCHIVE}/GWL-CO2only-B2030/pr_GWL-CO2only-B2030_0500-0600.nc'
TAS_TOTAL_FILENAME = f'{ARCHIVE}/GWL-NoCrops-B2030/pr_GWL-NoCrops-B2030_0500-0700.nc'


@cdod.cdo_selyear('540/569') # Load only the last 30 years.
@cdod.cdo_timmean
def load_last_30(input:str, var:str)->np.ndarray:
    """Load a map of the last 30 years of the simulation into a numpy array.
    """
    return cdo.copy(input=input, returnCdf=True).variables[var][:].squeeze()


@cdod.cdo_fldmean(weights='TRUE') # Global mean
def load_time_series(input:str, var:str)->np.ndarray:
    """Use CDO to do an area weighted global mean on a NetCDF file and load it into a numpy array.
    """
    return cdo.copy(input=input, returnCdf=True).variables[var][:].squeeze()



def sub_trim(a:np.ndarray, b:np.ndarray)->np.ndarray:
    """Trim arrays to the same length and subtract b from a.
    """
    if a.shape[0]>b.shape[0]:
        a = a[:b.shape[0]]
    elif b.shape[0]>a.shape[0]:
        b = b[:a.shape[0]]
    return a - b


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


# Load grid specification
examplenc = nc.Dataset(TAS_TOTAL_FILENAME, 'r')
lats = examplenc.variables['lat'][:]
lons = examplenc.variables['lon'][:]
examplenc.close()

# Temperature maps as anomaly from the reference experiment.
pr_ref = load_last_30(input=TAS_REF_FILENAME, var='pr')
pr_bgc = pr_ref - load_last_30(input=TAS_BGC_FILENAME, var='pr')
pr_total = pr_ref - load_last_30(input=TAS_TOTAL_FILENAME, var='pr')
# The biogeophysical component is the total (GWL-NoCrops-B2030) minus the biogeochemical component.
pr_bgp = pr_total - pr_bgc

# Temperture time series as anomaly from the reference experiment
pr_tseries_ref = load_time_series(input=TAS_REF_FILENAME, var='pr')
bgc_simulation = load_time_series(input=TAS_BGC_FILENAME, var='pr')
pr_tseries_bgc = sub_trim(pr_tseries_ref, bgc_simulation)
total_simulation = load_time_series(input=TAS_TOTAL_FILENAME, var='pr')
pr_tseries_total = sub_trim(pr_tseries_ref, total_simulation)
pr_tseries_bgp = sub_trim(pr_tseries_total, pr_tseries_bgc)

# Plot the time series of temperature
plt.figure()
plt.plot(yearly_mean_from_monthly(pr_tseries_total[:1200]), color='black', label='Total')
plt.plot(yearly_mean_from_monthly(pr_tseries_bgc), color='blue', label='Biogeochemical')
plt.plot(yearly_mean_from_monthly(pr_tseries_bgp), color='red', label='Biogeophysical')
plt.xlabel('Time (years)')
plt.ylabel('Precip (kg m$^{-2}$ s$^{-1}$)')
plt.legend()
#plt.show()


def plot_map(data:np.ndarray, title:str)->None:
    """Plot a global map of temperature.
    """
    plt.figure()
    ax = plt.axes(projection=ccrs.Robinson())
    discrete_bins = mpl.colors.BoundaryNorm(boundaries=np.arange(-4e-5, 4.5e-5, 0.5e-5), ncolors=256)
    shading = ax.pcolormesh(lons, lats, data,
            cmap='bwr_r',
            norm=discrete_bins,
            #vmin=-2.0,
            #vmax=2.0,
            transform=ccrs.PlateCarree(),
            )
    plt.title(title)
    cbar = plt.colorbar(shading,
            #ticks=[-2, -1.5, -1, -.5, .5, 1, 1.5, 2],
            orientation='horizontal',
            label='Precip (kg m$^{-2}$ s$^{-1}$)',
            pad=0.05,
            )
    ax.coastlines()


# Plot the map of temperature
plot_map(pr_bgc, 'Biogeochemical component')
plot_map(pr_bgp, 'Biogeophysical component')


def plot_australia(data:np.ndarray, title:str)->None:
    """Plot the Australian region for temperature.
    """
    plt.figure()
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_extent([110, 160, -45, -10], crs=ccrs.PlateCarree())
    discrete_bins = mpl.colors.BoundaryNorm(boundaries=np.arange(-4e-5, 4.5e-5, 0.5e-5), ncolors=256)
    shading = ax.pcolormesh(lons, lats, data,
            cmap='bwr_r',
            norm=discrete_bins,
            #vmin=-2.0,
            #vmax=2.0,
            transform=ccrs.PlateCarree(),
            )
    plt.title(title)
    cbar = plt.colorbar(shading,
            #ticks=[-2, -1.5, -1, -.5, .5, 1, 1.5, 2],
            orientation='horizontal',
            label='Precip (kg m$^{-2}$ s$^{-1}$)',
            pad=0.05,
            )
    ax.coastlines()


plot_australia(pr_bgc, 'Biogeochemical component')
plot_australia(pr_bgp, 'Biogeophysical component')

plt.figure(1)
plt.savefig('plots/pr_time_series_bgc-bgp_breakdown.png', dpi=200)
plt.figure(2)
plt.savefig('plots/pr_global_map_bgc.png', dpi=200)
plt.figure(3)
plt.savefig('plots/pr_global_map_bgp.png', dpi=200)
plt.figure(4)
plt.savefig('plots/pr_australia_map_bgc.png', dpi=200)
plt.figure(5)
plt.savefig('plots/pr_australia_map_bgp.png', dpi=200)

plt.show()
