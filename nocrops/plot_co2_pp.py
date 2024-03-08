#!/bin/env python3
import glob

import ipdb
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
from cdo import Cdo; cdo = Cdo()
import cdo_decorators as cdod
cdo.debug = True

@cdod.cdo_sellevidx('1/37')
@cdod.cdo_vertmean(weights='TRUE')
@cdod.cdo_fldmean(weights='TRUE')
def cdo_load(var:str, input:str):
    return cdo.copy(input=input, returnCdf=True, options='-L').variables[var][:].squeeze()


def load_co2(filename:str, var:str):
    """Load the surface co2 from a file.
    """
    return nc.Dataset(filename, 'r').variables[var][:,0,...].squeeze()


def yearly_mean(data):
    return np.mean(data.reshape((int(data.shape[-1]/12),12)), axis=1)


def global_mean(data:np.ndarray, lats:np.ndarray)->np.ndarray:
    """Calculate an area weighted global mean. Weights are the cosine of lats.
    """
    coslats = np.cos(np.deg2rad(lats))[:,None]*np.ones(data.shape)
    return np.ma.average(data, axis=(-1,-2), weights=coslats)


# Variables
CO2_FLD = 'fld_s00i252'
MOLMASS_AIR = 0.0289652 # kg/mol
MOLMASS_CO2 = 0.0440095 # kg/mol
KGKG_TO_MOLMOL = MOLMASS_AIR/MOLMASS_CO2 # Converion of kg/kg to mol/mol
MIL = 1000000
EXPERIMENTS = [
        #'esm-esm-piNoCrops',
        #'PI-GWL-t6',
        'PI-GWL-B2035',
        'PI-GWL-B2040',
        'PI-GWL-B2045',
        'PI-GWL-B2050',
        'PI-GWL-B2055',
        'PI-GWL-B2060',
        #'GWL-NoCrops-B2030',
        'GWL-NoCrops-B2035',
        'GWL-NoCrops-B2040',
        'GWL-NoCrops-B2045',
        'GWL-NoCrops-B2050',
        'GWL-NoCrops-B2055',
        'GWL-NoCrops-B2060',
        #'GWL-EqFor-B2060',
        ]
COLORS = {
        'esm-esm-piNoCrops':'black',
        'PI-GWL-t6':'#62EA00',
        'PI-GWL-B2035':'#24CC00',
        'PI-GWL-B2040':'#079F2A',
        'PI-GWL-B2045':'#00786B',
        'PI-GWL-B2050':'#055992',
        'PI-GWL-B2055':'#1140AB',
        'PI-GWL-B2060':'#1E31B6',
        'GWL-NoCrops-B2030':'#62EA00',
        'GWL-NoCrops-B2035':'#24CC00',
        'GWL-NoCrops-B2040':'#079F2A',
        'GWL-NoCrops-B2045':'#00786B',
        'GWL-NoCrops-B2050':'#055992',
        'GWL-NoCrops-B2055':'#1140AB',
        'GWL-NoCrops-B2060':'#1E31B6',
        'GWL-EqFor-B2060':'deepskyblue'
        }
YEARS = {
        'esm-esm-piNoCrops':np.arange(0*12, 100*12)/12 + 100,
        'PI-GWL-t6':np.arange(0*12, 700*12)/12,
        'PI-GWL-B2035':np.arange(400*12, 601*12)/12,
        'PI-GWL-B2040':np.arange(400*12, 601*12)/12,
        'PI-GWL-B2045':np.arange(400*12, 601*12)/12,
        'PI-GWL-B2050':np.arange(400*12, 601*12)/12,
        'PI-GWL-B2055':np.arange(400*12, 601*12)/12,
        'PI-GWL-B2060':np.arange(400*12, 601*12)/12,
        'GWL-NoCrops-B2030':np.arange(400*12, 601*12)/12,
        'GWL-NoCrops-B2035':np.arange(400*12, 601*12)/12,
        'GWL-NoCrops-B2040':np.arange(400*12, 601*12)/12,
        'GWL-NoCrops-B2045':np.arange(400*12, 601*12)/12,
        'GWL-NoCrops-B2050':np.arange(400*12, 601*12)/12,
        'GWL-NoCrops-B2055':np.arange(400*12, 601*12)/12,
        'GWL-NoCrops-B2060':np.arange(400*12, 601*12)/12,
        'GWL-EqFor-B2060':np.arange(400*12, 502*12)/12,
        }

load_nc = False

all_data = dict()
handles = list()
labels = list()
for exper in EXPERIMENTS:
    print(exper)
    if load_nc:
        filename = f'/g/data/p66/tfl561/archive_data/{exper}/co23D_{exper}_*.nc'
        all_data[exper] = cdo_load('co23D', input=filename)
        np.save(f'data/co2_mean_{exper}_all.npy', all_data[exper].data)
    else:
        all_data[exper] = np.load(f'data/co2_mean_{exper}_all.npy')

    x = yearly_mean(YEARS[exper])
    y = yearly_mean(all_data[exper])*KGKG_TO_MOLMOL*MIL

    lwidth = 1
    if 'PI' in exper:
        lwidth = 0.5
    handle = plt.plot(x, y,
            color=COLORS[exper],
            linewidth=lwidth,
            )
    if 'NoCrops' in exper or 'EqFor' in exper:
        handles.append(handle[0])
        labels.append(exper)

plt.legend(handles, labels, fontsize='small', frameon=False)
plt.xlabel('Time (years)')
plt.ylabel('CO$_2$ (ppm)')
plt.xlim(left=400, right=600)
plt.show()
#plt.savefig('plots/co2_surface_global_warming_level_no_crops.svg', dpi=200)

