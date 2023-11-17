#!/usr/bin/env python3
import glob

import netCDF4 as nc
import numpy as np
import f90nml
import ipdb

from contrib.tools import yearly_mean_from_monthly, global_mean

MISSINGVAL = -32768.0
ARCHIVE = '/g/data/p66/tfl561/ACCESS-ESM/GWL-NoCrops-B2030/history/atm/netCDF'

example_file = f'{ARCHIVE}/GWL-NoCrops-B2030.pa-050101_mon.nc'
ncfile = nc.Dataset(example_file, 'r')
thickness = ncfile.variables['theta_level_height_bnds'][:]
ipdb.set_trace()

files = sorted(glob.glob(f'{ARCHIVE}/GWL-NoCrops-B2030.pa-*.nc'))

co2_data = np.ones((len(files)))
for f in files:
    co2_month = global_mean(nc.Dataset(f, 'r').variables['fld_s00i252'][:,:-1,...].squeeze())
    co2_val = np.average()
#nyears = co2_yearly.size
nyears = 100
start_year = 500
years = np.arange(start_year, start_year+nyears)
#missing_data = np.ones((nyears))*MISSINGVAL

years = [[y for y in range(400,600)] for i in range(11)]
nyears = [200 for i in range(11)]
nyears[5] = None
nyears[6] = None

# Create namelist
namelist_dict = {
        'CLMCHFCG':{
                'L_CLMCHFCG':True,
                'CLIM_FCG_NYEARS':nyears,
                'CLIM_FCG_YEARS':years,
                },
        }
namelist = f90nml.Namelist(namelist_dict)
namelist.write('foo.nml')

