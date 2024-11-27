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

DPI = 200
G_IN_PG = 10**15

ARCHIVE_DIR = '/g/data/p66/tfl561/archive_data'
EXAMPLE_FLIE = f'{ARCHIVE_DIR}/GWL-NoCrops-B2030/cLeaf_GWL-NoCrops-B2030_0500-0700.nc'
LAND_FRAC_CODE = 'fld_s03i395'
LAST20 = str(-20*12)
RAW_CMIP_DIR = '/g/data/p73/archive/non-CMIP/ACCESS-ESM1-5'
RAW_NOCROP_DIR = '/g/data/p66/tfl561/ACCESS-ESM'
TILE_FRAC_CODE = 'fld_s03i317'

COLORS = {
        'GWL-NoCrops-B2030':'#62EA00',
        #'GWL-NoCrops-B2030':'pink',
        'GWL-NoCrops-B2035':'#24CC00',
        'GWL-NoCrops-B2040':'#079F2A',
        'GWL-NoCrops-B2045':'#00786B',
        'GWL-NoCrops-B2050':'#055992',
        'GWL-NoCrops-B2055':'#1140AB',
        'GWL-NoCrops-B2060':'#1E31B6',
        #'GWL-NoCr-B2030-02':'pink',
        'GWL-NoCr-B2030-02':'#62EA00',
        #'GWL-NoCrops-B2060':'red',
        'GWL-NoCr-B2030-02':'#62EA00',
        'GWL-NoCr-B2035-02':'#24CC00',
        'GWL-NoCr-B2040-02':'#079F2A',
        'GWL-NoCr-B2045-02':'#00786B',
        'GWL-NoCr-B2050-02':'#055992',
        'GWL-NoCr-B2055-02':'#1E31B6',
        'GWL-NoCr-B2060-02':'red',
        #'GWL-NoCr-B2060-02':'#1E31B6',
        'GWL-EqFor-B2060':'deepskyblue',
        }
EXPERIMENTS = {
        'PI-GWL-t6':'GWL-NoCrops-B2030',
        'PI-GWL-B2035':'GWL-NoCrops-B2035',
        'PI-GWL-B2040':'GWL-NoCrops-B2040',
        'PI-GWL-B2045':'GWL-NoCrops-B2045',
        'PI-GWL-B2050':'GWL-NoCrops-B2050',
        'PI-GWL-B2055':'GWL-NoCrops-B2055',
        'PI-GWL-B2060':'GWL-NoCrops-B2060',
        'PI-GWL-t6_duplicate':'GWL-NoCr-B2030-02',
        'PI-GWL-B2035_duplicate':'GWL-NoCr-B2035-02',
        'PI-GWL-B2040_duplicate':'GWL-NoCr-B2040-02',
        'PI-GWL-B2045_duplicate':'GWL-NoCr-B2045-02',
        'PI-GWL-B2050_duplicate':'GWL-NoCr-B2050-02',
        'PI-GWL-B2055_duplicate':'GWL-NoCr-B2055-02',
        'PI-GWL-B2060_duplicate':'GWL-NoCr-B2060-02',
        #'PI-GWL-B2060_duplicate':'GWL-EqFor-B2060',
        }
VARIABLES = {
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

def my_add(x1:np.ndarray, x2:np.ndarray)->np.ndarray:
    """Add x1 and x2, ignoring array size missmatch.
    """
    if len(x1)<len(x2):
        return x1 + x2[:len(x1)]
    elif len(x1)>len(x2):
        return x1[:len(x2)] + x2
    else:
        return x1 + x2


load_from_npy = True

# Cell area file
cdo.gridarea(input=EXAMPLE_FLIE, output='data/cell_area.nc')

# Load lats and lons.
ncfile = nc.Dataset(EXAMPLE_FLIE, 'r')
lats = ncfile.variables['lat'][:]
lons = ncfile.variables['lon'][:]

# Create tile frac file for both experiments.
#for exp in EXPERIMENTS.keys():
#    if 'duplicate' in exp: exp = exp[:-10]
#    cdo.selvar(
#            TILE_FRAC_CODE,
#            input=f'{RAW_CMIP_DIR}/{exp}/history/atm/netCDF/{exp}.pa-051001_mon.nc',
#            output=f'data/frac_{exp}.nc',
#            )
#for exp in EXPERIMENTS.values():
#    cdo.selvar(
#            TILE_FRAC_CODE,
#            input=f'{RAW_NOCROP_DIR}/{exp}/history/atm/netCDF/{exp}.pa-051001_mon.nc',
#            output=f'data/frac_{exp}.nc',
#            )

data = {}
data_tmean = {}
cLand_diff = {}
plt.figure()
for gwl_exp,nocrop_exp in EXPERIMENTS.items():
    if 'duplicate' in gwl_exp: gwl_exp = gwl_exp[:-10]
    for exp in [gwl_exp, nocrop_exp]:
        print(f"{exp}")
        data[exp] = {}
        data_tmean[exp] = {}


        @cdod.cdo_mul(input2=f'data/frac_{exp}.nc') # Multiply by tile fractions.
        @cdod.cdo_mul(input2='data/cell_area.nc') # Multiply by cell area.
        @cdod.cdo_sellonlatbox('110','155','-45','-10') # Australia region.
        @cdod.cdo_vertsum() # Plant functional type sum.
        @cdod.cdo_fldsum # Australia sum.
        def load_global_sum(var, input:str)->np.ma.MaskedArray:
            """Load Australia sum of a carbon pool variable.
            """
            ncfile = cdo.copy(input=input, returnCdf=True, options='-L')
            return ncfile.variables[var][:]


        # Load the data.
        period = '0500-0700'
        # The ensemble members restart file is for year 510.
        if '-02' in exp:
            period = '0510-0700'
        if load_from_npy:
            for var in VARIABLES.keys():
                data[exp][var] = np.load(f'data/{var}_{exp}_australia_sum.npy') # [g(C)]
        else:
            for var in VARIABLES.keys():
                print(f"Loading {var}")
                data[exp][var] = load_global_sum(
                        var,
                        input=f'{ARCHIVE_DIR}/{exp}/{var}_{exp}_*.nc',
                        )
                np.save(f'data/{var}_{exp}_australia_sum.npy', data[exp][var].data) # [g(C)]

        # Calculate cLand.
        data[exp]['cLand'] = 0
        for var in VARIABLES.keys():
            if var=='tas': continue
            data[exp]['cLand'] += data[exp][var]

    # Find the difference for cLand.
    # 2nd ensemble members start at year 510
    if '-02' in nocrop_exp:
        cLand_diff[nocrop_exp] = data[nocrop_exp]['cLand'] - data[gwl_exp]['cLand'][12*10:]
    else:
        cLand_diff[nocrop_exp] = data[nocrop_exp]['cLand'] - data[gwl_exp]['cLand']

    ## Plot the difference time series of cLand.
    plt.figure(1)
    nmonths = cLand_diff[nocrop_exp].squeeze().shape[0]
    years = np.linspace(500, 500+nmonths*(1/12), nmonths)
    plt.plot(years, cLand_diff[nocrop_exp].squeeze()/G_IN_PG, color=COLORS[nocrop_exp], label=exp) # [Pg(C)]
    plt.xlabel('Time (years)')
    plt.ylabel('$\Delta$ cLand (Pg)')
    plt.xlim(left=500)
    plt.ylim(bottom=0, top=20)
    plt.legend(frameon=False)

    # Plot the difference of all the carbon pools.
    plt.figure(3)
    for v in VARIABLES.keys():
        nmonths = cLand_diff[nocrop_exp].squeeze().shape[0]
        years = np.linspace(500, 500+nmonths*(1/12), nmonths)
        if data[nocrop_exp][v].shape!=data[gwl_exp][v].shape:
            plot_this = data[nocrop_exp][v].squeeze() - data[gwl_exp][v][120:].squeeze()
        else:
            plot_this = data[nocrop_exp][v].squeeze() - data[gwl_exp][v].squeeze()
        plt.plot(
                years,
                plot_this/G_IN_PG,
                color=COLORS[nocrop_exp],
                label=v,
                )
    plt.xlim(left=500)
    plt.xlabel('Time (years)')
    plt.ylabel('$\Delta$ C (Pg)')
    plt.title("Carbon pools")

# Plot the ensemble mean of the experiments.
plt.figure(2)
ensembles = [
        ['GWL-NoCrops-B2030','GWL-NoCr-B2030-02'],
        ['GWL-NoCrops-B2035','GWL-NoCr-B2035-02'],
        ['GWL-NoCrops-B2040','GWL-NoCr-B2040-02'],
        ['GWL-NoCrops-B2045','GWL-NoCr-B2045-02'],
        ['GWL-NoCrops-B2050','GWL-NoCr-B2050-02'],
        ['GWL-NoCrops-B2055','GWL-NoCr-B2055-02'],
        ['GWL-NoCrops-B2060','GWL-NoCr-B2060-02'],
        ]
for ens in ensembles:
    total = np.zeros(cLand_diff[ens[0]].shape)
    i = 0
    for exp in ens:
        nmonths = cLand_diff[exp].squeeze().shape[0]
        years = np.linspace(500, 500+nmonths*(1/12), nmonths)
        plt.plot(years, cLand_diff[exp].squeeze()/G_IN_PG, color=COLORS[ens[0]], alpha=0.2)
        total = my_add(total, cLand_diff[exp])
        i += 1
    plot_this = total/i/G_IN_PG
    nmonths = plot_this.squeeze().shape[0]
    years = np.linspace(500, 500+nmonths*(1/12), nmonths)
    plt.plot(years, plot_this.squeeze(), color=COLORS[ens[0]], label=ens[0][-4:])
plt.legend()
plt.xlim(left=500)
plt.ylim(bottom=0)
plt.xlabel('Time (years)')
plt.ylabel('$\Delta$cLand Pg(C)')
plt.title('Australia sum total land carbon')


# Save figures.
plt.figure(1)
plt.savefig(f'plots/cLand_GWL_australia_sum.png', dpi=DPI)
plt.figure(2)
plt.savefig(f'plots/cLand_GWL_australia_sum_ensemble_mean.png', dpi=DPI)
plt.figure(3)
plt.savefig(f'plots/cpools_australia_sum.png', dpi=DPI)
plt.show()

