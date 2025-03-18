#!/usr/bin/env python3
# Plot the difference of no-crops and with crops for the global warmining level simulations.
import glob
import os

import scipy.stats as stats
import cartopy.crs as ccrs
import cdo_decorators as cdod
import ipdb
import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np
from cdo import Cdo

cdo = Cdo()
cdo.debug = False

ARCHIVE_DIR = '/g/data/p66/tfl561/archive_data'
COLORS = {
        'GWL-NoCrops-B2030':'#62EA00',
        #'GWL-NoCrops-B2030':'pink',
        'GWL-NoCr-B2030-02':'#62EA00',
        #'GWL-NoCr-B2030-02':'pink',
        'GWL-NoCrops-B2035':'#24CC00',
        'GWL-NoCrops-B2040':'#079F2A',
        'GWL-NoCrops-B2045':'#00786B',
        'GWL-NoCrops-B2050':'#055992',
        'GWL-NoCrops-B2055':'#1140AB',
        'GWL-NoCrops-B2060':'#1E31B6',
        #'GWL-NoCrops-B2060':'red',
        'GWL-NoCr-B2060-02':'#1E31B6',
        #'GWL-NoCr-B2060-02':'red',
        'GWL-EqFor-B2060':'deepskyblue',
        }
CROP_INDEX = 8
DPI = 200
EXAMPLE_FLIE = f'{ARCHIVE_DIR}/GWL-NoCrops-B2030/tas_GWL-NoCrops-B2030_0500-0700.nc'
EXPERIMENTS = {
        'PI-GWL-t6':'GWL-NoCrops-B2030',
        'PI-GWL-t6_duplicate':'GWL-NoCr-B2030-02',
        'PI-GWL-B2035':'GWL-NoCrops-B2035',
        'PI-GWL-B2040':'GWL-NoCrops-B2040',
        'PI-GWL-B2045':'GWL-NoCrops-B2045',
        'PI-GWL-B2050':'GWL-NoCrops-B2050',
        'PI-GWL-B2055':'GWL-NoCrops-B2055',
        'PI-GWL-B2060':'GWL-NoCrops-B2060',
        'PI-GWL-B2060-duplicate':'GWL-NoCr-B2060-02',
        #'PI-GWL-B2060_duplicate':'GWL-EqFor-B2060',
        }
LAND_FRAC_CODE = 'fld_s03i395'
LAST100 = str(-100*12)
LOAD_FROM_NP = True
RAW_CMIP_DIR = '/g/data/p73/archive/non-CMIP/ACCESS-ESM1-5'
RAW_NOCROP_DIR = '/g/data/p66/tfl561/ACCESS-ESM'
TILE_FRAC_CODE = 'fld_s03i317'
VARIABLES = {
        'tas':'fld_s03i860',
        }


def decorate_plot()->None:
    """Decorate the current figure with labels.
    """
    plt.xlabel('Time (years)')
    plt.ylabel('$\Delta$ TAS ($^{\circ}$C)')
    plt.xlim(left=500, right=700)
    plt.hlines(y=0, xmin=500, xmax=700, color='black')
    plt.ylim(bottom=-1.5, top=1)
    plt.legend(frameon=False)


@cdod.cdo_sellonlatbox('110','155','-45','-10') # Australia region.
@cdod.cdo_fldmean(weights='TRUE') # Australia mean.
def load_australia_mean(var:str, input:str)->np.ma.MaskedArray:
    """Load Australia data and spatially average.
    """
    ncfile = cdo.copy(input=input, returnCdf=True, options='-L')
    return ncfile.variables[var][:]


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
@cdod.cdo_fldmean(weights='TRUE') # Australia mean.
def load_forestation_mean(var:str, input:str)->np.ndarray:
    """Load data in Australian forestation regions and spatially average.
    """
    ncfile = cdo.copy(input=input, returnCdf=True, options='-L')
    return ncfile.variables[var][:]


def moving_average(x:np.ndarray, window:int)->np.ndarray:
    """Calculate a moving average using a window of size `window`.
    """
    return np.convolve(x, np.ones(window), mode='same')/window


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


def main()->None:
    """Run the main routine.
    """
    # Load lats and lons.
    ncfile = nc.Dataset(EXAMPLE_FLIE, 'r')
    lats = ncfile.variables['lat'][:]
    lons = ncfile.variables['lon'][:]

    # Create mask for forestation only regions.
    make_forestation_mask(input=f'{ARCHIVE_DIR}/PI-GWL-t6/frac_PI-GWL-t6_0500-0700.nc')

    data = {}
    data_forest_only = {}
    plt.figure()
    for gwl_exp,nocrop_exp in EXPERIMENTS.items():
        if 'duplicate' in gwl_exp: gwl_exp = gwl_exp[:-10]
        for exp in [gwl_exp, nocrop_exp]:
            print(f"{exp}")
            data[exp] = {}
            data_forest_only[exp] = {}

            # Load the data.
            period = '0500-0700'
            # The ensemble members restart file is for year 510.
            if '-02' in exp: period = '0510-0700'
            if LOAD_FROM_NP:
                for var in VARIABLES.keys():
                    data[exp][var] = np.load(f'data/{var}_{exp}_australia_mean.npy')
                    data_forest_only[exp][var] = np.load(
                            f'data/{var}_{exp}_aus_forestaiton_only_mean.npy',
                            )
            else:
                for var in VARIABLES.keys():
                    print(f"Loading {var}")
                    # Load the Australia mean.
                    data[exp][var] = load_australia_mean(var,
                            input=f'{ARCHIVE_DIR}/{exp}/{var}_{exp}_{period}.nc',
                            )
                    np.save(f'data/{var}_{exp}_australia_mean.npy', data[exp][var].data)
                    # Load the australia forestation region mean.
                    data_forest_only[exp][var] = load_forestation_mean(var,
                            input=f'{ARCHIVE_DIR}/{exp}/{var}_{exp}_{period}.nc',
                            )
                    np.save(
                            f'data/{var}_{exp}_aus_forestaiton_only_mean.npy',
                            data_forest_only[exp][var].data,
                            )

        # Find the difference for tas.
        # 2nd ensemble members start at year 510
        if '-02' in nocrop_exp:
            tas_diff = data[nocrop_exp]['tas'] - data[gwl_exp]['tas'][12*10:]
            tas_for_diff = data_forest_only[nocrop_exp]['tas'] - \
                    data_forest_only[gwl_exp]['tas'][12*10:]
        else:
            tas_diff = data[nocrop_exp]['tas'] - data[gwl_exp]['tas']
            tas_for_diff = data_forest_only[nocrop_exp]['tas'] - data_forest_only[gwl_exp]['tas']

        ## Plot the difference time series of tas.
        plt.figure(1)
        nmonths = tas_diff.squeeze().shape[0]
        years = np.linspace(500, 500+nmonths*(1/12), nmonths)
        plt.plot(yearly_mean_from_monthly(years), yearly_mean_from_monthly(tas_diff.squeeze()),
                color=COLORS[nocrop_exp],
                alpha=0.3,
                )
        if '-' in exp[-4:]:
            label = ''
        else:
            label = exp[-4:]
        plt.plot(
                yearly_mean_from_monthly(years),
                moving_average(yearly_mean_from_monthly(tas_diff.squeeze()), window=30),
                color=COLORS[nocrop_exp],
                label=label,
                )

        # Plot the difference for forestation regions only.
        plt.figure(2)
        nmonths = tas_diff.squeeze().shape[0]
        years = np.linspace(500, 500+nmonths*(1/12), nmonths)
        yearly_tas_for_diff = yearly_mean_from_monthly(tas_for_diff.squeeze())
        plt.plot(yearly_mean_from_monthly(years), yearly_tas_for_diff,
                color=COLORS[nocrop_exp],
                alpha=0.3,
                )
        plt.plot(yearly_mean_from_monthly(years), moving_average(yearly_tas_for_diff, window=30),
                color=COLORS[nocrop_exp],
                label=label,
                )

    """Significance testing the lines. I use here the Wilcoxon non-parametric test for the
    difference of paired means. There are 9 samples of the reference simulations and 9
    samples of the experiment simulations. This would test that the temperature change
    from forestation (under a range of CO2 fertilizations) is significantly different from the
    case where there is no forestation. Significance occurs where the null hypothesis that there
    is no effect on temperature is rejected at the 5% level.
    ...I need to rethink this.
    """
    NEXP = 9
    NTIMES = len(yearly_mean_from_monthly(data['PI-GWL-t6']['tas']))
    references = np.ones((NEXP,NTIMES))
    experiments = references.copy()
    i = 0
    for ref,exp in EXPERIMENTS.items():
        if 'duplicate' in ref: ref = ref[:-10]
        length = len(yearly_mean_from_monthly(data[exp]['tas'].squeeze()))
        references[i,:] = yearly_mean_from_monthly(data[ref]['tas'].squeeze())
        experiments[i,:length] = yearly_mean_from_monthly(data[exp]['tas'].squeeze())
        i += 1
    test_results = stats.wilcoxon(experiments, references, axis=0, nan_policy='omit')

    plt.figure(1)
    decorate_plot()
    plt.title('Australia mean')
    plt.savefig(f'plots/tas_GWL_australia.png', dpi=DPI)

    plt.figure(2)
    decorate_plot()
    plt.title('Australia forestation area mean')
    plt.savefig(f'plots/tas_GWL_australia_forestation_only.png', dpi=200)
    plt.show()

if __name__=='__main__':
    main()
    plt.show()

