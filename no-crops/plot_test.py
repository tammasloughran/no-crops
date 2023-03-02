#!/usr/bin/env python3
# Plot the test simulations of the no crop experiment.
import glob
import os

import ipdb
import numpy as np
import cdo_decorators as cdod
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import netCDF4 as nc
from cdo import Cdo

cdo = Cdo()
cdo.debug = True

TILE_FRAC_VAR = 'fld_s03i317'
WOOD_VAR = 'fld_s03i853'
NOCROP_ARCHIVE_DIR = '/g/data/p66/tfl561/ACCESS-ESM'
PROCESSED_NOCROP_DIR = '/g/data/p66/tfl561/archive_data'
EXPERIMENTS = [
        'esm-esm-piNoCrops',
        #'old_esm-piNoCrops',
        #'esm-esm-pre-industrial',
        ]
UM_DATA = 'history/atm/netCDF'
NOCROPS_VARIABLES = [
        'rh',
        'cCoarseWoodyDebris',
        'cLeaf',
        'cMetabolic',
        'cMicrobial',
        'cPassive',
        'cRoot',
        'cSlow',
        'cStructural',
        'cWood',
        'tas',
        ]
PI_L_VARIABLES = [
        'cLeaf',
        'cLitter',
        'cRoot',
        'cCwd',
        'rh',
        ]
PI_E_VARIABLES = [
        'cWood',
        'cSoil',
        ]
PI_A_VARIABLES = [
        'tas',
        ]
PLOT_VARS = ['rh','cLeaf','cWood','cLitter','cSoil','cRoot','cCwd','tas']
#PI_DIR = '/g/data/fs38/publications/CMIP6/CMIP/CSIRO/ACCESS-ESM1-5/esm-piControl/r1i1p1f1'
PI_DIR = '/g/data/p73/archive/CMIP6/ACCESS-ESM1-5/PI-EDC-01/history/atm/netCDF'
LAND_FRAC = '/g/data/fs38/publications/CMIP6/CMIP/CSIRO/ACCESS-ESM1-5/esm-piControl/r1i1p1f1' \
        '/fx/sftlf/gn/latest/sftlf_fx_ACCESS-ESM1-5_esm-piControl_r1i1p1f1_gn.nc'


@cdod.cdo_selvar(TILE_FRAC_VAR)
@cdod.cdo_mul(input1='cell_area.nc')
@cdod.cdo_mul(input2=LAND_FRAC)
@cdod.cdo_divc('100')
@cdod.cdo_ifthen(input1=LAND_FRAC)
def cdo_get_tile_areas(input:str, output:str)->None:
    cdo.copy(input=input, output=output, options='-L')


@cdod.cdo_cat(input2='')
@cdod.cdo_divc('1e15')
@cdod.cdo_mul(input2='tile_areas.nc')
@cdod.cdo_fldsum
#@cdod.cdo_vertsum # I've decided to do this in python, because I want the pools on PFTs.
def cdo_load_global_sum(input:str, varname:str)->np.ndarray:
    """Global sum using cdo and load.
    """
    return cdo.copy(input=input, returnCdf=True, options='-L').variables[varname][:].squeeze()


@cdod.cdo_cat(input2='')
@cdod.cdo_mul(input2=LAND_FRAC)
@cdod.cdo_divc('1e15')
@cdod.cdo_fldsum
def cdo_load_global_sum2(input:str, varname:str)->np.ndarray:
    """Global sum using cdo and load for rh.
    """
    return cdo.copy(input=input, returnCdf=True, options='-L').variables[varname][:].squeeze()


@cdod.cdo_cat(input2='')
@cdod.cdo_fldmean()
def cdo_load_global_mean(input:str, varname:str)->np.ndarray:
    return cdo.copy(input=input, returnCdf=True, options='-L').variables[varname][:].squeeze()


@cdod.cdo_cat(input2='')
@cdod.cdo_seltimestep('-20/-1')
@cdod.cdo_timmean
def cdo_load_temp_last(input:str, varname:str):
    return cdo.copy(input=input, returnCdf=True, options='-L').variables[varname][:].squeeze()


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
    example_file = f'{NOCROP_ARCHIVE_DIR}/{exp}/{UM_DATA}/{exp}.pa-010101_mon.nc'
    cdo.gridarea(input=example_file, output='cell_area.nc')
    cdo_get_tile_areas(input=example_file, output='tile_areas.nc')

    # Load the variables for the no-crop experiment.
    no_crops = {}
    for var in NOCROPS_VARIABLES:
        if not (var=='tas' or var=='rh'):
            no_crops[var] = cdo_load_global_sum(
                    input=f'[ {PROCESSED_NOCROP_DIR}/{exp}/{var}_{exp}.nc ]',
                    varname=var,
                    )
        elif var=='rh':
            no_crops[var] = cdo_load_global_sum2(
                    input=f'[ {PROCESSED_NOCROP_DIR}/{exp}/{var}_{exp}.nc ]',
                    varname=var,
                    )
        else: # Load the sureface temperature variable as a global mean.
            no_crops[var] = cdo_load_global_mean(
                    input=f'[ {PROCESSED_NOCROP_DIR}/{exp}/{var}_{exp}.nc ]',
                    varname=var,
                    )

    # Aggregate to the CMIP variables for cSoil and cLitter.
    #no_crops['cSoil'] = no_crops['cMicrobial'] + no_crops['cSlow'] + no_crops['cPassive']
    #no_crops['cLitter'] = no_crops['cMetabolic'] + no_crops['cStructural'] + \
    #        no_crops['cCoarseWoodyDebris']
    #no_crops['cCwd'] = no_crops['cCoarseWoodyDebris']


    # PI-EDC-01.pa-052304_mon.nc
    # Recreate the tile areas for the pre-industrial simulation.
    exp = 'PI-EDC-01'
    example_file = f'{PI_DIR}/PI-EDC-01.pa-010101_mon.nc'
    cdo_get_tile_areas(input=example_file, output='tile_areas.nc')

    # Load the pre-industrial variables according to the relevant table.
    pi_data = {}
    for var in NOCROPS_VARIABLES:
        if not (var=='tas' or var=='rh'):
            pi_data[var] = cdo_load_global_sum(
                    input=f'[ {PROCESSED_NOCROP_DIR}/{exp}/{var}_{exp}.nc ]',
                    varname=var,
                    )
        elif var=='rh':
            pi_data[var] = cdo_load_global_sum2(
                    input=f'[ {PROCESSED_NOCROP_DIR}/{exp}/{var}_{exp}.nc ]',
                    varname=var,
                    )
        else:
            pi_data[var] = cdo_load_global_mean(
                    input=f'[ {PROCESSED_NOCROP_DIR}/{exp}/{var}_{exp}.nc ]',
                    varname=var,
                    )

    for var in NOCROPS_VARIABLES:
        plt.figure()
        if not (var=='tas' or var=='rh'):
            plt.plot(no_crops[var].sum(axis=1), label='esm-piNoCrops')
            plt.plot(pi_data[var].sum(axis=1), label='esm-piControl')
        elif var=='rh':
            nocrop = no_crops[var]#.reshape(12,50).mean(axis=0)
            pidata = pi_data[var]#.reshape(12,50).mean(axis=0)
            plt.plot(nocrop, label='esm-piNoCrops')
            plt.plot(pidata, label='esm-piControl')
        else:
            plt.plot(no_crops[var], label='esm-piNoCrops')
            plt.plot(pi_data[var], label='esm-piControl')
        plt.title(var)
        plt.xlabel('Time (months)')
        if var=='tas':
            plt.ylabel('°C')
        else:
            plt.ylabel('Pg C')
        plt.legend()
    plt.show()

    pft_names = {
            0:'EvgNL',
            1:'EvgBL',
            2:'DecNL',
            3:'DecBL',
            4:'Shrub',
            5:'C3 grass',
            6:'C4 grass',
            7:'Tundra',
            8:'C3 Crop',
            }

    # Plot the difference for all PFTs.
    for var in NOCROPS_VARIABLES:
        if var=='tas' or var=='rh': continue
        plt.figure()
        for pft in range(9):
            plt.plot(no_crops[var][:,pft] - pi_data[var][:,pft], label=f'{pft_names[pft]}')
            plt.title(var)
            plt.xlabel('Time (months)')
            if var=='tas':
                plt.ylabel('°C')
            else:
                plt.ylabel('Pg C')
            plt.legend()
    plt.show()

    exp = 'esm-esm-piNoCrops'
    no_crops_temp = cdo_load_temp_last(input=f'[ {PROCESSED_NOCROP_DIR}/{exp}/tas_{exp}.nc ]',
            varname='tas')
    exp = 'PI-EDC-01'
    pi_temp = cdo_load_temp_last(input=f'[ {PROCESSED_NOCROP_DIR}/{exp}/tas_{exp}.nc ]',
            varname='tas')

    ncin = nc.Dataset(f'{PROCESSED_NOCROP_DIR}/{exp}/tas_{exp}.nc')
    lats = ncin.variables['lat_v'][:]
    lons = ncin.variables['lon_u'][:]
    # pcolormesh expects +1 lon.
    llons = np.append(lons, [lons[-1] + lons[1] - lons[0]])
    # Shift plotting to centre grid, align with coastlines.
    llons = llons - (llons[1] - llons[0])/2
    llats = lats - (lats[1] - lats[0])/4

    # Plot the map of the difference in temperature.
    plt.figure()
    ax = plt.axes(projection=ccrs.Robinson())
    plt.pcolormesh(
            llons,
            llats,
            no_crops_temp - pi_temp,
            cmap='seismic',
            vmin=-10,
            vmax=10,
            transform=ccrs.PlateCarree(),
            )
    ax.coastlines()
    plt.colorbar(label='$\Delta$ Temperature (°C)', orientation='horizontal')
    plt.title('Surface temperature esm-piNoCrops - esm-piControl')
    plt.show()

