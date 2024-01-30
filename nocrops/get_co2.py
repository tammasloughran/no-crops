#!/usr/env python3
# CO2 data from emissions driven simulations are output for (time,lev,lat,lon) which is too much
# data for long  monthly simulations. This script loads files one at a time using the python
# netcdf4 library and manually does global averaging with area and level thickness weghtings.
# Output files are .npy files in 100 year chunks.
# To run for each 100 year chunk:
# for i in 0 1 2 3 4 5 6; do {
#     nohup python get_co2.py ${i}
# } &
# done

import glob
import sys

import ipdb
import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np


def global_mean(data:np.ndarray, lats:np.ndarray)->np.ndarray:
    """Calculate an area weighted global mean. Weights are the cosine of lats.
    """
    coslats = np.cos(np.deg2rad(lats))[:,None]*np.ones(data.shape)
    return np.ma.average(data, axis=(-1,-2), weights=coslats)


# Variables
sect = int(sys.argv[-1])
CO2_FLD = 'fld_s00i252'
MOLMASS_AIR = 0.0289652 # kg/mol
MOLMASS_CO2 = 0.0440095 # kg/mol
KGKG_TO_MOLMOL = MOLMASS_AIR/MOLMASS_CO2 # Converion of kg/kg to mol/mol
MIL = 1000000
NYEARS = 12*100
#exp = 'PI-GWL-B2060'
#exp_dir = f'/g/data/p73/archive/non-CMIP/ACCESS-ESM1-5/{exp}/history/atm/netCDF'
exp = 'esm-esm-piNoCrops'
exp_dir = f'/g/data/p66/tfl561/ACCESS-ESM/{exp}/history/atm/netCDF'
files = sorted(glob.glob(f'{exp_dir}/{exp}.pa-*_mon.nc'))[:800*12]
files = files[0+sect*NYEARS:NYEARS+sect*NYEARS]

# Load
sample_nc = nc.Dataset(files[0], 'r')
lats = sample_nc.variables['lat'][:]

co2_data = np.empty((len(files),))
for i,f in enumerate(files):
    print('\r', i, '/', len(files), end='')
    ncfile = nc.Dataset(f, 'r')
    co2_data[i] = global_mean(ncfile.variables[CO2_FLD][:,0,...].squeeze(), lats)
    ncfile.close()

# Do vertical weighted average.
#layer_thickness = np.diff(sample_nc.variables['theta_level_height_bnds'][:]).squeeze()
#vertical_weights = (layer_thickness/layer_thickness.sum())[None,:]*np.ones(co2_data.shape)
#co2_concentration = np.average(co2_data, axis=1, weights=vertical_weights)*KGKG_TO_MOLMOL*MIL

# Save to file
np.save(f'data/co2_surface_{exp}_global_mean_{sect}.npy', co2_data.data)

# Ignore all this.

# Plot
#plt.figure()
#data = np.mean(co2_data, axis=1)
#years = np.arange(1, data.size+1)/12.0
#plt.plot(years, data)
#plt.xlim(left=0)
#plt.xlabel('Time (years)')
#plt.ylabel('CO$_2$ concentration (kg/kg)')
#plt.show()
#
##
#use_cdo = False
#
#if use_cdo:
#    import cdo_decorators as cdod
#    from cdo import Cdo
#
#    cdo = Cdo()
#    cdo.debug = False
#
#    @cdod.cdo_fldmean(weights='TRUE')
#    @cdod.cdo_vertmean(weights='TRUE')
#    def load_cdo_mean(input:str, var:str)->np.ndarray:
#        """Do a spatial and vertical mean and load the result into python as a numpy array.
#        """
#        return cdo.copy(input=input, returnCdf=True).variables[var][:].squeeze()
#co2var = sample_nc.variables[CO2_FLD]
#co2_attr = vars(co2var)
#if use_cdo:
#    co2_data = np.empty((len(files),))
#    for i,f in enumerate(files):
#        #print('\r', i, '/', len(files), end='')
#        co2_data[i] = load_cdo_mean(input=f, var=CO2_FLD)
#if not use_cdo:
#    # Save to netCDF.
#    ncout = nc.Dataset(f'co21D_{sect}.nc', 'w')
#    ncout.createDimension('time', size=None)
#    ncout.createDimension('model_theta_level_number', size=38)
#    ncout.createDimension('bnds', size=2)
#
#    timeout = ncout.createVariable('time', 'i', dimensions=('time',))
#    timeout[:] = np.arange(1+sect*1200, 1201+sect*1200)/12.0
#
#    thetavar = sample_nc.variables['model_theta_level_number']
#    thetaout = ncout.createVariable(
#            'model_theta_level_number',
#            'f',
#            dimensions=('model_theta_level_number',),
#            )
#    for att,val in vars(thetavar).items():
#       thetaout.setncattr(att, val)
#    thetaout[:] = thetavar[:]
#
#    co2out = ncout.createVariable(
#            'co2',
#            'f',
#            dimensions=('time','model_theta_level_number',),
#            fill_value=co2var._FillValue,
#            )
#    for att,val in vars(co2var).items():
#        if att=='_FillValue': continue
#        co2out.setncattr(att, val)
#    co2out[:] = co2_data
#    ncout.close()
