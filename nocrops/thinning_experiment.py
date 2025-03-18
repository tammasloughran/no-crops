#!/usr/bin/env python3
""" Plot the results from the GWL-Thin-B2030.

This is complicated because:
- Natural carbon pools do not have tile or land fractions included in the variable.
- Wood product pools already have tile and land fractions included in the variable.
- Wood product pools are yearly and the natural pools are monthly.
"""
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import ipdb

READ_ONLY = 'r'
BASE_VARS = [
        'cLeaf',
        'cWood',
        'cRoot',
        'cMetabolic',
        'cStructural',
        'cCoarseWoodyDebris',
        'cMicrobial',
        'cSlow',
        'cPassive',
        ]
HARVEST_VARS = [
        'wood_harvest_1',
        'wood_harvest_2',
        'wood_harvest_3',
        ]
ALL_VARS = BASE_VARS + HARVEST_VARS
EXPERIMENTS = [
        'PI-GWL-t6',
        'GWL-NoCrops-B2030',
        'GWL-Thin-B2030',
        ]
BASE_DIR = '/g/data/p66/tfl561/archive_data'
AUS_LONS = slice(60, 85)
AUS_LATS = slice(35, 65)
LAND_FRAC_FILE = '/g/data/fs38/publications/CMIP6/CMIP/CSIRO/ACCESS-ESM1-5/historical/r1i1p1f1/'\
    'fx/sftlf/gn/latest/sftlf_fx_ACCESS-ESM1-5_historical_r1i1p1f1_gn.nc'
LAND_FRAC_VARNAME = 'sftlf'


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


def load_var(infile, varname):
    print(f'Loading {varname} from {infile}')
    ncfile = nc.Dataset(infile, READ_ONLY)
    if varname=='frac':
        return_data = ncfile.variables[varname][0,:,AUS_LATS,AUS_LONS].squeeze()
    else:
        return_data = ncfile.variables[varname][:,:,AUS_LATS,AUS_LONS].squeeze()
    ncfile.close()
    return return_data


# Load the land fraction and cell_area
ncfile = nc.Dataset(LAND_FRAC_FILE, READ_ONLY)
land_frac = ncfile.variables[LAND_FRAC_VARNAME][AUS_LATS,AUS_LONS].squeeze()/100.0
ncfile.close()
ncfile = nc.Dataset(f'{BASE_DIR}/GWL-Thin-B2030/cell_area.nc', READ_ONLY)
cell_area = ncfile.variables['cell_area'][AUS_LATS,AUS_LONS].squeeze()
ncfile.close()

# Load the experiments data
load_nc = False
if load_nc:
    data = {}
    for exper in EXPERIMENTS:
        tile_frac = load_var(f'{BASE_DIR}/{exper}/frac_{exper}_0500-0700.nc', 'frac')

        if 'Thin' not in exper:
            data[exper] = 0.0
            for var in BASE_VARS:
                data[exper] += load_var(f'{BASE_DIR}/{exper}/{var}_{exper}_0500-0700.nc', var)
            data[exper] = yearly_mean_from_monthly(data[exper])
            data[exper] = data[exper]*tile_frac*land_frac*cell_area
        else:
            data[exper] = 0.0
            for var in ALL_VARS:
                tmp = load_var(f'{BASE_DIR}/{exper}/{var}_{exper}_0500-0700.nc', var)
                if 'harvest' in var:
                    data[exper] += tmp[1:]*cell_area
                else:
                    data[exper] += yearly_mean_from_monthly(tmp*tile_frac*land_frac*cell_area)

    standard_no_crops = np.nansum(data['GWL-NoCrops-B2030'] - data['PI-GWL-t6'], axis=(1,2,3))
    with_thinning = np.nansum(data['GWL-Thin-B2030'] - data['PI-GWL-t6'], axis=(1,2,3))

    np.save('standard.npy', standard_no_crops.data)
    np.save('with_thinning.npy', with_thinning.data)
else:
    standard_no_crops = np.load('standard.npy')
    with_thinning = np.load('with_thinning.npy')

years = np.arange(500,701)
plt.plot(years, standard_no_crops/1E15, label='2030 GWL')
plt.plot(years, with_thinning/1E15, label='wood thinning')
plt.xlabel('Time (years)')
plt.ylabel('cLand [Pg(C)]')
plt.legend()
plt.show()
