#!/usr/bin/env python3
# Plot the test pre-industrial simulations of the no crop experiment.
import glob
import os

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cdo_decorators as cdod
import ipdb
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np
from cdo import Cdo

cdo = Cdo()
cdo.debug = False

TILE_FRAC_VAR = 'fld_s03i317'
WOOD_VAR = 'fld_s03i853'
NOCROP_ARCHIVE_DIR = '/g/data/p66/tfl561/ACCESS-ESM'
PROCESSED_NOCROP_DIR = '/g/data/p66/tfl561/archive_data'
EXPERIMENTS = [
        'esm-esm-piNoCrops',
        'esm-esm-piNoCrops-2',
        #'old_esm-piNoCrops',
        #'esm-esm-pre-industrial',
        ]
COLORS = {
        'esm-piNoCrops':'blue',
        'esm-piNoCrops-2':'lightblue',
        'esm-piControl':'gray',
        }
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

load_cdo = False

@cdod.cdo_selvar(TILE_FRAC_VAR)
@cdod.cdo_mul(input1='data/cell_area.nc')
@cdod.cdo_mul(input2=LAND_FRAC)
@cdod.cdo_divc('100')
@cdod.cdo_ifthen(input1=LAND_FRAC)
def cdo_get_tile_areas(input:str, output:str)->None:
    cdo.copy(input=input, output=output, options='-L')


@cdod.cdo_cat(input2='')
@cdod.cdo_divc('1e15')
@cdod.cdo_mul(input2='data/tile_areas.nc')
@cdod.cdo_fldsum
#@cdod.cdo_vertsum # I've decided to do this in python, because I want the pools on PFTs.
def cdo_load_global_sum(input:str, varname:str)->np.ndarray:
    """Global sum using cdo and load data into numpy arrays.
    """
    return cdo.copy(input=input, returnCdf=True, options='-L').variables[varname][:].squeeze()


@cdod.cdo_cat(input2='')
@cdod.cdo_mul(input2='data/cell_area.nc')
@cdod.cdo_mul(input2=LAND_FRAC) # From m-2 to -
@cdod.cdo_divc('100') # From % to frac
@cdod.cdo_mulc('86400')
@cdod.cdo_muldpm # From s-1 to mon-1
@cdod.cdo_divc('1e12') # From kg to Pg
@cdod.cdo_fldsum
def cdo_load_global_sum2(input:str, varname:str)->np.ndarray:
    """Global sum using cdo and load for carbon flux variables.
    """
    return cdo.copy(input=input, returnCdf=True, options='-L').variables[varname][:].squeeze()


@cdod.cdo_cat(input2='')
@cdod.cdo_fldmean()
def cdo_load_global_mean(input:str, varname:str)->np.ndarray:
    return cdo.copy(input=input, returnCdf=True, options='-L').variables[varname][:].squeeze()


@cdod.cdo_cat(input2='')
@cdod.cdo_seltimestep('-240/-1') # last 20 years (20*12)
@cdod.cdo_timmean
def cdo_load_temp_last(input:str, varname:str):
    return cdo.copy(input=input, returnCdf=True, options='-L').variables[varname][:].squeeze()


@cdod.cdo_seltimestep('-360/-1') # last 20 years (30*12)
@cdod.cdo_selseason('DJF')
@cdod.cdo_timmean
@cdod.cdo_mul(input2='data/tile_areas.nc')
@cdod.cdo_divc('1e15')
def cdo_load_djf_mean(input:str, varname:str)->np.ndarray:
    """Load a spatial map of DJF 30-year composites for carbon pools.
    """
    return cdo.copy(input=input, returnCdf=True, options='-L').variables[varname][:].squeeze()


@cdod.cdo_seltimestep('-360/-1') # last 20 years (30*12)
@cdod.cdo_selseason('DJF')
@cdod.cdo_timmean
def cdo_load_djf_mean2(input:str, varname:str)->np.ndarray:
    """Load a spatial map of DJF 30-year composites for temp.
    """
    return cdo.copy(input=input, returnCdf=True, options='-L').variables[varname][:].squeeze()


@cdod.cdo_mul(input2='data/cell_area.nc')
@cdod.cdo_mul(input2=LAND_FRAC) # From m-2 to -
@cdod.cdo_divc('100') # From % to frac
@cdod.cdo_mulc('86400')
@cdod.cdo_muldpm # From s-1 to mon-1
@cdod.cdo_divc('1e12') # From kg to Pg
@cdod.cdo_seltimestep('-360/-1') # last 20 years (30*12)
@cdod.cdo_selseason('DJF')
@cdod.cdo_timmean
def cdo_load_djf_mean3(input:str, varname:str)->np.ndarray:
    """Load a spatial map of DJF 30-year composites for flux.
    """
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
    # PI-EDC-01.pa-052304_mon.nc
    # Recreate the tile areas for the pre-industrial simulation.
    exp = 'PI-EDC-01'
    example_file = f'{PI_DIR}/PI-EDC-01.pa-010101_mon.nc'
    cdo_get_tile_areas(input=example_file, output='data/tile_areas.nc')

    # Load the pre-industrial variables according to the relevant table.
    pi_data = {}
    pi_djf_data = {}
    for var in NOCROPS_VARIABLES:
        if load_cdo==True:
            files = sorted(glob.glob(f'{PROCESSED_NOCROP_DIR}/{exp}/{var}_{exp}*.nc'))
            if not (var=='tas' or var=='rh'):
                pi_data[var] = cdo_load_global_sum(
                        input='[ '+' '.join(files)+' ]',
                        varname=var,
                        )
                pi_djf_data[var] = cdo_load_djf_mean(
                        input=files[-1],
                        varname=var,
                        )
            elif var=='rh':
                pi_data[var] = cdo_load_global_sum2(
                        input='[ '+' '.join(files)+' ]',
                        varname=var,
                        )
                pi_djf_data[var] = cdo_load_djf_mean3(
                        input=files[-1],
                        varname=var,
                        )
            else:
                pi_data[var] = cdo_load_global_mean(
                        input='[ '+' '.join(files)+' ]',
                        varname=var,
                        )
                pi_djf_data[var] = cdo_load_djf_mean2(
                        input=files[-1],
                        varname=var,
                        )
            np.save(f'data/{var}_{exp}_test.npy', pi_data[var].data)
            np.save(f'data/{var}_{exp}_djf_test.npy', pi_djf_data[var].data)
        else:
            pi_data[var] = np.load(f'data/{var}_{exp}_test.npy')
            pi_djf_data[var] = np.load(f'data/{var}_{exp}_djf_test.npy')

    # Calculate tile areas. They are the same for both pre-industrial ensemble members.
    exp = 'esm-esm-piNoCrops'
    example_file = f'{NOCROP_ARCHIVE_DIR}/{exp}/{UM_DATA}/{exp}.pa-010101_mon.nc'
    cdo.gridarea(input=example_file, output='data/cell_area.nc')
    cdo_get_tile_areas(input=example_file, output='data/tile_areas.nc')

    # Load the variables for the no-crop experiment.
    no_crops = {}
    no_crops_djf = {}
    for e in [1,2]:
        no_crops[e] = {}
        no_crops_djf[e] = {}
        if e==1:
            exp = 'esm-esm-piNoCrops'
        else:
            exp = 'esm-esm-piNoCrops-2'
        for var in NOCROPS_VARIABLES:
            if load_cdo==True:
                files = sorted(glob.glob(f'{PROCESSED_NOCROP_DIR}/{exp}/{var}_{exp}*.nc'))
                if not (var=='tas' or var=='rh'):
                    no_crops[e][var] = cdo_load_global_sum(
                            input='[ '+' '.join(files)+' ]',
                            varname=var,
                            )
                    no_crops_djf[e][var] = cdo_load_djf_mean(
                            input=files[-1],
                            varname=var,
                            )
                elif var=='rh':
                    no_crops[e][var] = cdo_load_global_sum2(
                            input='[ '+' '.join(files)+' ]',
                            varname=var,
                            )
                    no_crops_djf[e][var] = cdo_load_djf_mean3(
                            input=files[-1],
                            varname=var,
                            )
                else: # Load the sureface temperature variable as a global mean.
                    no_crops[e][var] = cdo_load_global_mean(
                            input='[ '+' '.join(files)+' ]',
                            varname=var,
                            )
                    no_crops_djf[e][var] = cdo_load_djf_mean2(
                            input=files[-1],
                            varname=var,
                            )
                np.save(f'data/{var}_{exp}_test.npy', no_crops[e][var].data)
                np.save(f'data/{var}_{exp}_djf_test.npy', no_crops_djf[e][var].data)
            else:
                no_crops[e][var] = np.load(f'data/{var}_{exp}_test.npy')
                no_crops_djf[e][var] = np.load(f'data/{var}_{exp}_djf_test.npy')


    # Aggregate to the CMIP variables for cSoil and cLitter.
    #no_crops['cSoil'] = no_crops['cMicrobial'] + no_crops['cSlow'] + no_crops['cPassive']
    #no_crops['cLitter'] = no_crops['cMetabolic'] + no_crops['cStructural'] + \
    #        no_crops['cCoarseWoodyDebris']
    #no_crops['cCwd'] = no_crops['cCoarseWoodyDebris']


    def yearly_mean(data:np.ndarray):
        """Calculate a yearly mean on a numpy array.
        The first dimension must be time and divisible by 12.
        """
        toshape = list(data.shape)
        toshape.pop(0)
        toshape.insert(0, 12)
        toshape.insert(0, int(data.shape[0]/12))
        return data.reshape(toshape).mean(axis=1)


#    # Plot the carbon pools.
#    for var in NOCROPS_VARIABLES:
#        plt.figure()
#        if var[0]=='c':
#            plt.plot(
#                    yearly_mean(no_crops[1][var].sum(axis=1)),
#                    label='esm-piNoCrops',
#                    color=COLORS['esm-piNoCrops'],
#                    )
#            plt.plot(
#                    yearly_mean(no_crops[2][var].sum(axis=1)),
#                    label='esm-piNoCrops-2',
#                    color=COLORS['esm-piNoCrops-2'],
#                    )
#            plt.plot(
#                    yearly_mean(pi_data[var].sum(axis=1)),
#                    label='esm-piControl',
#                    color=COLORS['esm-piControl'],
#                    )
#        else:
#            plt.plot(
#                    yearly_mean(no_crops[1][var]),
#                    label='esm-piNoCrops',
#                    color=COLORS['esm-piNoCrops'],
#                    )
#            plt.plot(
#                    yearly_mean(no_crops[2][var]),
#                    label='esm-esm-piNoCrops-2',
#                    color=COLORS['esm-piNoCrops-2'],
#                    )
#            plt.plot(
#                    yearly_mean(pi_data[var]),
#                    label='esm-piControl',
#                    color=COLORS['esm-piControl'],
#                    )
#        plt.title(var)
#        plt.xlabel('Time (months)')
#        if var=='tas':
#            plt.ylabel('°C')
#        else:
#            plt.ylabel('Pg C')
#        plt.legend()
#        plt.xlim(left=0)
#        plt.savefig(f'plots/{var}_timeseries.png', dpi=200)
#    plt.show()
#
#    # Plot as a % change of the pi control.
#    for var in NOCROPS_VARIABLES:
#        if var[0]!='c': continue
#        plt.figure()
#        no_crop1 = yearly_mean(no_crops[1][var].sum(axis=1))
#        no_crop2 = yearly_mean(no_crops[2][var].sum(axis=1))
#        picont1 = yearly_mean(pi_data[var][:len(no_crop1)*12].sum(axis=1))
#        picont2 = yearly_mean(pi_data[var][:len(no_crop2)*12].sum(axis=1))
#        plt.plot(
#                (no_crop1/picont1)*100,
#                label='esm-piNoCrops',
#                color=COLORS['esm-piNoCrops'],
#                )
#        plt.plot(
#                (no_crop2/picont2)*100,
#                label='esm-piNoCrops-2',
#                color=COLORS['esm-piNoCrops-2'],
#                )
#        plt.title(var)
#        plt.xlabel('Time (months)')
#        if var=='tas':
#            plt.ylabel('°C')
#        else:
#            plt.ylabel('Pg C')
#        plt.legend()
#        plt.xlim(left=0)
#        plt.savefig(f'plots/{var}_percent_change_timeseries.png', dpi=200)
#    plt.show()


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
    color_cycle = [k[4:] for k in  mcolors.TABLEAU_COLORS.keys()]
    colors = {pft:color_cycle[i] for i,pft in pft_names.items()}
    light_colors = [
            'lightblue',
            'wheat',
            'lightgreen',
            'tomato',
            'mediumorchid',
            'chocolate',
            'lightpink',
            'lightgray',
            'yellowgreen',
            ]

#    # Plot the difference for all PFTs.
#    for var in NOCROPS_VARIABLES:
#        if var=='tas' or var=='rh': continue
#        plt.figure()
#        handles = list()
#        for pft in range(9):
#            for e in [1,2]:
#                if e==1:
#                    color = colors[pft_names[pft]]
#                else:
#                    color = light_colors[pft]
#                lenexp = len(no_crops[e][var][:,pft])
#                handles.append(plt.plot(
#                        no_crops[e][var][:,pft] - pi_data[var][:lenexp,pft],
#                        label=f'{pft_names[pft]}',
#                        color=color
#                        )[0])
#                if e==2: handles.pop()
#        plt.title(var)
#        plt.xlabel('Time (months)')
#        plt.ylabel('Pg C')
#        plt.legend(handles, pft_names.values())
#        plt.savefig(f'plots/{var}_difference_pfts.png', dpi=200)
#    plt.show()

    exp = 'esm-esm-piNoCrops-2'
    files = sorted(glob.glob(f'{PROCESSED_NOCROP_DIR}/{exp}/tas_{exp}*.nc'))
    no_crops_temp = cdo_load_temp_last(input=files[-1], varname='tas')
    exp = 'PI-EDC-01'
    files = sorted(glob.glob(f'{PROCESSED_NOCROP_DIR}/{exp}/tas_{exp}*.nc'))
    pi_temp = cdo_load_temp_last(input=files[-1], varname='tas')

    #ncin = nc.Dataset(f'{PROCESSED_NOCROP_DIR}/{exp}/tas_{exp}_0101-0150.nc')
    ncin = nc.Dataset(f'{PROCESSED_NOCROP_DIR}/{exp}/tas_{exp}_0101-0150.nc')
    lats = ncin.variables['lat_v'][:]
    lons = ncin.variables['lon_u'][:]
    # pcolormesh expects +1 lon.
    llons = np.append(lons, [lons[-1] + lons[1] - lons[0]])
    # Shift plotting to centre grid, align with coastlines.
    llons = llons - (llons[1] - llons[0])/2
    llats = lats - (lats[1] - lats[0])/4

    # Plot the maps of djf.
    for var in NOCROPS_VARIABLES:
        plt.figure()
        ax = plt.axes(projection=ccrs.Robinson())
        if var[0]=='c':
            data = np.nansum(no_crops_djf[1][var], axis=0) - np.nansum(pi_djf_data[var], axis=0)
            data[data>1e15] = np.nan
        else:
            data = no_crops_djf[1][var] - pi_djf_data[var]
        mappable = ax.pcolormesh(
                llons,
                llats,
                data,
                cmap='seismic',
                transform=ccrs.PlateCarree(),
                vmin=-1.0*max(abs(data.min()), abs(data.max())),
                vmax=max(abs(data.min()), abs(data.max())),
                )
        ax.coastlines()
        plt.colorbar(mappable, label=f'$\Delta$ {var} (Pg)', orientation='horizontal')
        plt.savefig(f'plots/{var}_esm-esm-piNoCrops-1_djf.png', dpi=200)

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
    plt.savefig(f'plots/temperature_diff.png', dpi=200)
    plt.show()

