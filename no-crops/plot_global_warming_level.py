#!/usr/bin/env python3
# Plot the difference of no-crops and with crops for the global warmining level simulations.
import glob
import os

import cartopy.crs as ccrs
import cdo_decorators as cdod
import ipdb
import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np
from cdo import Cdo

cdo = Cdo()
cdo.debug = False

VARIABLES = {
        #'tas':'fld_s03i236',
        'cLeaf':'fld_s03i852',
        'cWood':'fld_s03i853',
        'cRoot':'fld_s03i854',
        'cMetabolic':'fld_s03i855',
        'cStructural':'fld_s03i856',
        'cCoarseWoodyDebris':'fld_s03i857',
        'cMicrobial':'fld_s03i858',
        'cSlow':'fld_s03i859',
        'cPassive':'fld_s03i860',
        }
EXPERIMENTS = {
        'PI-GWL-t6':'GWL-NoCrops-B2030',
        'PI-GWL-B2035':'GWL-NoCrops-B2035',
        'PI-GWL-B2040':'GWL-NoCrops-B2040',
        'PI-GWL-B2045':'GWL-NoCrops-B2045',
        'PI-GWL-B2050':'GWL-NoCrops-B2050',
        'PI-GWL-B2055':'GWL-NoCrops-B2055',
        'PI-GWL-B2060':'GWL-NoCrops-B2060',
        'PI-GWL-B2060_duplicate':'GWL-EqFor-B2060',
        }
ARCHIVE_DIR = '/g/data/p66/tfl561/archive_data'
RAW_CMIP_DIR = '/g/data/p73/archive/non-CMIP/ACCESS-ESM1-5'
RAW_NOCROP_DIR = '/g/data/p66/tfl561/ACCESS-ESM'
LAND_FRAC_CODE = 'fld_s03i395'
TILE_FRAC_CODE = 'fld_s03i317'
EXAMPLE_FLIE = f'{ARCHIVE_DIR}/GWL-NoCrops-B2030/cLeaf_GWL-NoCrops-B2030_0500-0700.nc'
G_IN_PG = 10**15
DPI = 200
LAST30 = str(-30*12)
COLORS = {
        'GWL-NoCrops-B2030':'#62EA00',
        'GWL-NoCrops-B2035':'#24CC00',
        'GWL-NoCrops-B2040':'#079F2A',
        'GWL-NoCrops-B2045':'#00786B',
        'GWL-NoCrops-B2050':'#055992',
        'GWL-NoCrops-B2055':'#1140AB',
        'GWL-NoCrops-B2060':'#1E31B6',
        'GWL-EqFor-B2060':'deepskyblue',
        }

load_from_npy = False

# Cell area file
cdo.gridarea(input=EXAMPLE_FLIE, output='data/cell_area.nc')

# Load lats and lons.
ncfile = nc.Dataset(EXAMPLE_FLIE, 'r')
lats = ncfile.variables['lat'][:]
lons = ncfile.variables['lon'][:]

# Create tile frac file for both experiments.
for exp in EXPERIMENTS.keys():
    if 'duplicate' in exp: exp = 'PI-GWL-B2060'
    cdo.selvar(
            TILE_FRAC_CODE,
            input=f'{RAW_CMIP_DIR}/{exp}/history/atm/netCDF/{exp}.pa-050101_mon.nc',
            output=f'data/frac_{exp}.nc',
            )
for exp in EXPERIMENTS.values():
    cdo.selvar(
            TILE_FRAC_CODE,
            input=f'{RAW_NOCROP_DIR}/{exp}/history/atm/netCDF/{exp}.pa-050101_mon.nc',
            output=f'data/frac_{exp}.nc',
            )

data = {}
data_tmean = {}
plt.figure()
for gwl_exp,nocrop_exp in EXPERIMENTS.items():
    if 'duplicate' in gwl_exp: gwl_exp = 'PI-GWL-B2060'
    for exp in [gwl_exp, nocrop_exp]:
        print(f"{exp}")
        data[exp] = {}
        data_tmean[exp] = {}


        @cdod.cdo_mul(input2=f'data/frac_{exp}.nc') # Multiply by tile fractions.
        @cdod.cdo_mul(input2='data/cell_area.nc') # Multiply by cell area.
        @cdod.cdo_vertsum(weights='FALSE') # Plant functional type sum.
        @cdod.cdo_fldsum # Global sum.
        def load_global_sum(var, input:str)->np.ma.MaskedArray:
            """Load global sum of a carbon pool variable.
            """
            ncfile = cdo.copy(input=input, returnCdf=True, options='-L')
            return ncfile.variables[var][:]


        @cdod.cdo_mul(input2=f'data/frac_{exp}.nc') # Multiply by tile fractions.
        @cdod.cdo_mul(input2='data/cell_area.nc') # Multiply by cell area.
        @cdod.cdo_vertsum(weights='FALSE') # Plant functional type sum.
        @cdod.cdo_seltimestep(f'{LAST30}/-1') # last 20 years
        @cdod.cdo_timmean # Temporal average
        def load_last20(var, input:str)->np.ma.MaskedArray:
            """Load last 20 year mean of a carbon pool variable.
            """
            ncfile = cdo.copy(input=input, returnCdf=True, options='-L')
            return ncfile.variables[var][:]


        # Load the data.
        if 'B2030' in exp or 't6' in exp:
            period = '0500-0700'
        else:
            period = '0500-0601'
        if load_from_npy:
            for var in VARIABLES.keys():
                data[exp][var] = np.load(f'data/{var}_{exp}_global_sum.npy') # [g(C)]
                data_tmean[exp][var] = np.load(f'data/{var}_{exp}_last20.npy')
        else:
            for var in VARIABLES.keys():
                print(f"Loading {var}")
                data[exp][var] = load_global_sum(
                        var,
                        input=f'{ARCHIVE_DIR}/{exp}/{var}_{exp}_{period}.nc',
                        )
                np.save(f'data/{var}_{exp}_global_sum.npy', data[exp][var].data) # [g(C)]
                data_tmean[exp][var] = load_last20(
                        var,
                        input=f'{ARCHIVE_DIR}/{exp}/{var}_{exp}_{period}.nc',
                        )
                np.save(f'data/{var}_{exp}_last20.npy', data_tmean[exp][var].data) # [g(C)]

        # Calculate cLand.
        data[exp]['cLand'] = 0
        data_tmean[exp]['cLand'] = 0
        for var in VARIABLES.keys():
            if var=='tas': continue
            data[exp]['cLand'] += data[exp][var]
            data_tmean[exp]['cLand'] += data_tmean[exp][var]

    # Find the difference for cLand.
    cLand_diff = data[nocrop_exp]['cLand'] - data[gwl_exp]['cLand']
    cLand_diff_tmean = data_tmean[nocrop_exp]['cLand'] - data_tmean[gwl_exp]['cLand']

    ## Plot the difference time series of cLand.
    plt.figure(1)
    nmonths = cLand_diff.squeeze().shape[0]
    years = np.linspace(500, 500+nmonths*(1/12), nmonths)
    plt.plot(years, cLand_diff.squeeze()/G_IN_PG, color=COLORS[nocrop_exp], label=exp) # [Pg(C)]
    plt.xlabel('Time (years)')
    plt.ylabel('$\Delta$ cLand (Pg)')
    plt.xlim(left=500)
    plt.ylim(bottom=0, top=200)
    plt.legend(frameon=False)

    # Plot the map of the difference for the last 20 years
    fig2 = plt.figure()
    ax = fig2.add_subplot(1, 1, 1, projection=ccrs.Robinson())
    colors = ax.pcolormesh(
            lons,
            lats,
            cLand_diff_tmean.squeeze()/G_IN_PG,
            cmap='seismic',
            vmin=-0.6,
            vmax=0.6,
            transform=ccrs.PlateCarree(),
            )
    ax.coastlines()
    plt.colorbar(colors, label='$\Delta$ cLand (Pg)', orientation='horizontal', pad=0.05)
    plt.title(f'{nocrop_exp} - {gwl_exp}')
    plt.savefig(f'plots/cLand_{nocrop_exp}_last20.png', dpi=DPI)

    # Plot the difference of all the carbon pools.
    plt.figure(3)
    for v in VARIABLES.keys():
        nmonths = cLand_diff.squeeze().shape[0]
        years = np.linspace(500, 500+nmonths*(1/12), nmonths)
        plot_this = data[nocrop_exp][v].squeeze() - data[gwl_exp][v].squeeze()
        plt.plot(
                years,
                plot_this/G_IN_PG,
                color=COLORS[nocrop_exp],
                #label=v,
                )
    plt.xlim(left=500)
    plt.xlabel('Time (years)')
    plt.ylabel('$\Delta$ C (Pg)')
    plt.annotate('cWood', (675,125))
    plt.annotate('cRoot', (675,40))
    plt.annotate('cCWD', (675,20))
    plt.annotate('cStructural', (675,10))
    plt.annotate('cMicrobial', (675,-10))
    plt.annotate('cSlow', (675,-25))
    plt.title("Carbon pools")

    # Plot the map of the difference for the last 20 years in Australia.
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    ax.set_extent([110,155,-45,-10], crs=ccrs.PlateCarree())
    colors = ax.pcolormesh(
            lons,
            lats,
            cLand_diff_tmean.squeeze()/G_IN_PG,
            cmap='seismic',
            vmin=-0.6,
            vmax=0.6,
            transform=ccrs.PlateCarree(),
            )
    ax.coastlines()
    plt.colorbar(colors, label='$\Delta$ cLand (Pg)]', orientation='horizontal', pad=0.05)
    plt.title(f'{nocrop_exp} - {gwl_exp}')
    plt.savefig(f'plots/cLand_{nocrop_exp}_australia_last20.png', dpi=DPI)

# Save figures.
plt.figure(3)
plt.savefig(f'plots/cpools_gloabl_sum.png', dpi=DPI)
plt.figure(1)
plt.savefig(f'plots/cLand_GWL_gloabl_sum.png', dpi=DPI)
plt.show()

