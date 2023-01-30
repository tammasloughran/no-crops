#!/usr/bin/env/python3
# Plot the test simulations of the no crop experiment.
import glob
import os

import cdo_decorators as cdod
import matplotlib.pyplot as plt
import netCDF4 as nc
from cdo import Cdo

cdo = Cdo()
cdo.debug = True

TILE_FRAC_VAR = 'fld_s03i317'
WOOD_VAR = 'fld_s03i853'
ARCHIVE_DIR = '/g/data/p66/tfl561/ACCESS-ESM'
EXPERIMENTS = [
        #'esm-esm-piFewCrops',
        'esm-esm-piNoCrops',
        #'esm-esm-pi',
        #'old_esm-piNoCrops',
        'esm-esm-pre-industrial',
        #'esm-esm-pi1850Luc',
        #'esm-historical',
        #'esm-historical-no-luc',
        ]
UM_DATA = 'history/atm/netCDF'


@cdod.cdo_selvar(TILE_FRAC_VAR)
@cdod.cdo_mul(input1='cell_area.nc')
@cdod.cdo_ifthen(input1='fractions.nc')
def cdo_get_tile_areas(input:str, output:str):
    cdo.copy(input=input, output=output, options='-L')


@cdod.cdo_divc('1e15')
@cdod.cdo_mul(input2='tile_areas.nc')
@cdod.cdo_sellevel('1','2','3','4')
@cdod.cdo_fldsum
@cdod.cdo_vertsum
def cdo_load_global_sum(input:str, varname:str):
    return cdo.copy(input=input, returnCdf=True, options='-L').variables[varname][:].squeeze()


def nc_selvar(filename, varname, outfile):
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
    # Create plot canvas
    plt.figure()

    for exp in EXPERIMENTS:
        # Calculate tile areas.
        example_file = f'{ARCHIVE_DIR}/{exp}/{UM_DATA}/{exp}.pa-010101_mon.nc'
        cdo.gridarea(input=example_file, output='cell_area.nc')
        cdo_get_tile_areas(input=example_file, output='tile_areas.nc')

        print("Create concatenated wood only file for", exp)
        ncfiles = sorted(glob.glob(f'{ARCHIVE_DIR}/{exp}/{UM_DATA}/*pa*.nc'))
        for ncfile in ncfiles:
            nc_selvar(ncfile, WOOD_VAR, f'wood_pa{ncfile[-14:]}')
        wood_files = sorted(glob.glob('./wood_*pa*.nc'))
        cdo.cat(input=' '.join(wood_files), output=f'wood_{exp}.nc', options='-L')
        for ncfile in wood_files: os.remove(ncfile)

        print("Calculate global sum of wood pool")
        wood_pool = cdo_load_global_sum(input=f'wood_{exp}.nc', varname=WOOD_VAR)

        # Plot the line for this experiment
        plt.plot(wood_pool, label=exp)

    plt.legend()
    plt.xlabel('Time (months)')
    plt.ylabel('Wood pool (PgC)')
    plt.show()

