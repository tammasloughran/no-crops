#!/usr/bin/env python3
import argparse
import sys

import netCDF4
import numpy as np
import umfile
from um_fileheaders import *
import ipdb
import netCDF4 as nc
import matplotlib.pyplot as plt

ORIGINAL_RESTART = '/g/data/p66/tfl561/restart_dump.astart'
NEW_RESTART = 'myrestart.astart'
VEGFRAC_CODE = 216
PREV_VEGFRAC_CODE = 835
CWOOD_CODE = 853
NTILES=17

original_restart_file = umfile.UMFile(ORIGINAL_RESTART, 'r')
new_restart_file = umfile.UMFile(NEW_RESTART, 'r')


def getvar(infile, code):
    var = []
    for k in range(infile.fixhd[FH_LookupSize2]):
        ilookup = infile.ilookup[k]
        lbegin = ilookup[LBEGIN]
        if lbegin==-99:
            break
        if ilookup[ITEM_CODE]==code:
            var.append(infile.readfld(k))
    assert len(var)==NTILES, 'Error - expected %d vegetation classes' % NTILES
    var = np.array(var)
    var[var==infile.missval_r] = np.nan
    return var


old_vegfrac = getvar(original_restart_file, VEGFRAC_CODE)
old_cwood = getvar(original_restart_file, CWOOD_CODE)
new_vegfrac = getvar(new_restart_file, VEGFRAC_CODE)
new_previous_year = getvar(new_restart_file, PREV_VEGFRAC_CODE)
new_cwood = getvar(new_restart_file, CWOOD_CODE)

area = nc.Dataset('gridarea.nc').variables['cell_area'][:]

def total_pool(density, fraction, cell_area):
    return np.nansum(density*fraction*cell_area)

ipdb.set_trace()
access_pa = nc.Dataset('/g/data/p66/tfl561/ACCESS-ESM/esm-esm-piNoCrops/history/atm/netCDF/esm-esm-piNoCrops.pa-010101_mon.nc', 'r')
access_density = access_pa.variables['fld_s03i853'][:].squeeze()
access_pool = access_density*new_vegfrac*area
old_pool = old_cwood*old_vegfrac*area
difference = old_cwood - access_density
plt.pcolormesh(difference[0]); plt.colorbar(); plt.show()

print(total_pool(old_cwood, old_vegfrac, area))

print(total_pool(new_cwood))

