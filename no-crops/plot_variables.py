#!/usr/bin/env python3
# Plot the test simulations of the no crop experiment.
import glob
import os
import re

import ipdb
import numpy as np
import cdo_decorators as cdod
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import netCDF4 as nc
from cdo import Cdo

cdo = Cdo()
cdo.debug = False

TILE_FRAC_VAR = 'fld_s03i317'
WOOD_VAR = 'fld_s03i853'
NOCROP_ARCHIVE_DIR = '/g/data/p66/tfl561/ACCESS-ESM'
PROCESSED_NOCROP_DIR = '/g/data/p66/tfl561/archive_data'
EXPERIMENTS = [
        'esm-esm-piNoCrops',
        'PI-EDC-01',
        ]
UM_DATA = 'history/atm/netCDF'
NOCROPS_VARIABLES = [
        'hfss', # Global mean time series, 20 year mean map.
        'hfls', # Global mean time series, 20 year mean map.
        'tas', # Global mean time series, 20 year mean map.
        'rss', # Global mean time series, 20 year mean map.
        'rls', # global mean time series, 20 year mean map.
        'gpp', # Multiply by land area, and global sum. 20 year mean
        'npp', # Multiply by land area, and global sum. 20 year mean
        'ra', # Multiply by land area, and global sum. 20 year mean
        'pr', # Global mean time series, 20 year mean map.
        'clt', # Global mean time series and 20 year mean map.
        'mrso', # Multiply by land area, global sum, 20 year mean.
        'rh',
        ]
CMAP = {
        'hfss':'seismic',
        'hfls':'PuOr',
        'tas':'seismic',
        'rss':'bone',
        'rls':'magma',
        'gpp':'BrBG',
        'npp':'BrBG',
        'ra':'PiYG',
        'rh':'PiYG',
        'pr':'RdBu',
        'clt':'bone',
        'mrso':'PuOr',
        }
TITLE = {
        'hfss':'Sensible heat flux',
        'hfls':'Latent heat flux',
        'tas':'Surface air temperature',
        'rss':'Net shortwave',
        'rls':'Net longwave',
        'gpp':'Gross primary production',
        'npp':'Net primary production',
        'ra':'Plant respiration',
        'rh':'Soil respiration',
        'pr':'Precipitation rate',
        'clt':'Cloud fraction',
        'mrso':'Total soil moisture',
        }
PLOT_VARS = NOCROPS_VARIABLES
PI_DIR = '/g/data/p73/archive/CMIP6/ACCESS-ESM1-5/PI-EDC-01/history/atm/netCDF'
LAND_FRAC = '/g/data/fs38/publications/CMIP6/CMIP/CSIRO/ACCESS-ESM1-5/esm-piControl/r1i1p1f1' \
        '/fx/sftlf/gn/latest/sftlf_fx_ACCESS-ESM1-5_esm-piControl_r1i1p1f1_gn.nc'


def pretty_units(units:str)->str:
    """Convert a units string into latex format.
    """
    units = f'${units}$'
    units = units.replace(' ', ' \cdot ')
    nums = re.findall('-*\d+', units)
    sups = ['^{'+s+'}' for s in nums]
    for num,sup in zip(nums, sups):
        units = units.replace(num, sup)
    return units


@cdod.cdo_mul(input1='cell_area.nc')
@cdod.cdo_divc('100')
@cdod.cdo_ifthen(input1=LAND_FRAC)
def cdo_get_land_areas(input:str, output:str)->None:
    cdo.copy(input=input, output=output, options='-L')


@cdod.cdo_mulc('86.4') # kg s-1 to tonnes day-1
@cdod.cdo_mul(input2='land_areas.nc')
@cdod.cdo_fldsum
def cdo_load_global_sum(input:str, varname:str)->np.ndarray:
    """Global sum using cdo and load.
    """
    ncfile = cdo.copy(input=input, returnCdf=True, options='-L')
    data = ncfile.variables[varname][:].squeeze()
    units = ncfile.variables[varname].units
    return data, units


@cdod.cdo_fldmean()
def cdo_load_global_mean(input:str, varname:str)->np.ndarray:
    ncfile = cdo.copy(input=input, returnCdf=True, options='-L')
    data = ncfile.variables[varname][:].squeeze()
    units = ncfile.variables[varname].units
    return data, units


@cdod.cdo_seltimestep('-240/-1')
@cdod.cdo_timmean
def cdo_load_temp_last(input:str, varname:str):
    ncfile = cdo.copy(input=input, returnCdf=True, options='-L')
    data = ncfile.variables[varname][:].squeeze()
    units = ncfile.variables[varname].units
    return data, units


@cdod.cdo_seltimestep('-240/-1')
@cdod.cdo_mulc('86.4') # kg s-1 to tonnes day-1
@cdod.cdo_mul(input2='land_areas.nc')
@cdod.cdo_timmean
def cdo_load_temp_last2(input:str, varname:str):
    ncfile = cdo.copy(input=input, returnCdf=True, options='-L')
    data = ncfile.variables[varname][:].squeeze()
    units = ncfile.variables[varname].units
    return data, units


def nc_selvar(filename:str, varname:str, outfile:str)->None:
    """Copy a variable from a netCDF file.
    """
    print(outfile)
    ncfile = nc.Dataset(filename, 'r')
    var = ncfile.variables[varname]
    ncout = nc.Dataset(outfile, 'w')
    for dim in var.dimensions:
        if 'time' in dim:
            ncout.createDimension(dim, size=None)
        else:
            ncout.createDimension(dim, size=ncfile.dimensions[dim].size)
        newdim = ncout.createVariable(dim, var.datatype, dim)
        newdim[:] = ncfile.variables[dim][:]
    outvar = ncout.createVariable(
            var.name,
            var.datatype,
            var.dimensions,
            fill_value=var.getncattr('_FillValue'),
            )
    attributes = ncfile.variables[varname].ncattrs()
    attributes.remove('_FillValue')
    for attr in attributes:
        outvar.setncattr(attr, var.getncattr(attr))
    outvar[:] = var[:]
    ncout.close()
    ncfile.close()


if __name__=='__main__':
    # Calculate tile areas.
    exp = 'esm-esm-piNoCrops'
    cdo.gridarea(input=LAND_FRAC, output='cell_area.nc')
    cdo_get_land_areas(input=LAND_FRAC, output='land_areas.nc')

    # Load the variables for the no-crop experiment.
    no_crops = {}
    var_units = {}
    for var in NOCROPS_VARIABLES:
        if var in ['gpp','npp','ra','rh']:
            no_crops[var], var_units[var] = cdo_load_global_sum(
                    input=f'{PROCESSED_NOCROP_DIR}/{exp}/{var}_{exp}.nc',
                    varname=var,
                    )
        else: # Load the other variables as globam mean.
            no_crops[var], var_units[var] = cdo_load_global_mean(
                    input=f'{PROCESSED_NOCROP_DIR}/{exp}/{var}_{exp}.nc',
                    varname=var,
                    )

    # PI-EDC-01.pa-052304_mon.nc
    # Recreate the land areas for the pre-industrial simulation.
    exp = 'PI-EDC-01'
    cdo_get_land_areas(input=LAND_FRAC, output='land_areas.nc')

    # Load the pre-industrial variables according to the relevant table.
    pi_data = {}
    for var in NOCROPS_VARIABLES:
        if var in ['gpp','npp','ra','rh']:
            pi_data[var], _ = cdo_load_global_sum(
                    input=f'{PROCESSED_NOCROP_DIR}/{exp}/{var}_{exp}.nc',
                    varname=var,
                    )
        else:
            pi_data[var], _ = cdo_load_global_mean(
                    input=f'{PROCESSED_NOCROP_DIR}/{exp}/{var}_{exp}.nc',
                    varname=var,
                    )

    # Units
    no_crops['pr'] *= 60*60*24
    pi_data['pr'] *= 60*60*24
    var_units['pr'] = 'kg m-2 day-1'
    var_units['gpp'] = 'tonnes day-1'
    var_units['npp'] = 'tonnes day-1'
    var_units['ra'] = 'tonnes day-1'
    var_units['rh'] = 'tonnes day-1'

    # Make units latex pretty.
    for key,val in var_units.items():
        var_units[key] = pretty_units(val)

    for var in NOCROPS_VARIABLES:
        plt.figure()
        difference = no_crops[var] - pi_data[var]
        yearly = difference.reshape((50,12)).mean(axis=1)
        plt.plot(np.linspace(0, 51, 600, endpoint=False), difference)
        plt.plot(range(1, 51), yearly, color='navy')
        plt.hlines(y=0, xmin=0, xmax=51, colors='black')
        plt.title(TITLE[var])
        plt.xlabel('Time (year)')
        if var in ['gpp','npp','ra','rh']:
            plt.ylabel(f'Tonnes day-1')
        else:
            plt.ylabel(var_units[var])
        plt.xlim(left=0, right=51)
        plt.savefig(f'plots/{var}_global_esm-piControl_esm-piNoCrops_difference.svg')
    #plt.show()

    no_crops_last = {}
    pi_last = {}
    for var in NOCROPS_VARIABLES:
        if var in ['gpp','npp','ra','rh']:
            exp = 'esm-esm-piNoCrops'
            no_crops_last[var], _ = cdo_load_temp_last2(
                    input=f'{PROCESSED_NOCROP_DIR}/{exp}/{var}_{exp}.nc',
                    varname=var,
                    )
            exp = 'PI-EDC-01'
            pi_last[var], _ = cdo_load_temp_last2(
                    input=f'[ {PROCESSED_NOCROP_DIR}/{exp}/{var}_{exp}.nc ]',
                    varname=var,
                    )
        else:
            exp = 'esm-esm-piNoCrops'
            no_crops_last[var], _ = cdo_load_temp_last(
                    input=f'{PROCESSED_NOCROP_DIR}/{exp}/{var}_{exp}.nc',
                    varname=var,
                    )
            exp = 'PI-EDC-01'
            pi_last[var], _ = cdo_load_temp_last(
                    input=f'[ {PROCESSED_NOCROP_DIR}/{exp}/{var}_{exp}.nc ]',
                    varname=var,
                    )

    ncin = nc.Dataset(f'{PROCESSED_NOCROP_DIR}/{exp}/tas_{exp}.nc')
    lats = ncin.variables['lat_v'][:]
    lons = ncin.variables['lon_u'][:]
    # pcolormesh expects +1 lon.
    llons = np.append(lons, [lons[-1] + lons[1] - lons[0]])
    # Shift plotting to centre grid, align with coastlines.
    llons = llons - (llons[1] - llons[0])/2
    llats = lats - (lats[1] - lats[0])/4

    for var in NOCROPS_VARIABLES:
        # Plot the map of the difference in temperature.
        plt.figure()
        ax = plt.axes(projection=ccrs.Robinson())
        data = no_crops_last[var] - pi_last[var]
        maxrange = max(np.max(data), np.abs(np.min(data)))
        plt.pcolormesh(
                llons,
                llats,
                no_crops_last[var] - pi_last[var],
                cmap=CMAP[var],
                vmin=-maxrange,
                vmax=maxrange,
                transform=ccrs.PlateCarree(),
                linewidth=0.2,
                edgecolors='face',
                )
        ax.coastlines()
        cbar = plt.colorbar(label=f'{var_units[var]}', orientation='horizontal', pad=0.05)
        cbar.solids.set_edgecolor('face')
        plt.title(f'{TITLE[var]} difference years 31-50')
        plt.tight_layout()
        plt.savefig(f'plots/{var}_map_esm-piControl_esm-piNoCrops_difference.svg')
    plt.show()

