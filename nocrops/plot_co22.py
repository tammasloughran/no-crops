#!/bin/env python3
import glob

import ipdb
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt


def load_co2(filename:str, var:str):
    """Load the surface co2 from a file.
    """
    return nc.Dataset(filename, 'r').variables[var][:,0,...].squeeze()


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
        'PI-GWL-t6',
        'PI-GWL-B2035',
        'PI-GWL-B2040',
        'PI-GWL-B2045',
        'PI-GWL-B2050',
        'PI-GWL-B2055',
        'PI-GWL-B2060',
        'GWL-NoCrops-B2030',
        'GWL-NoCrops-B2035',
        'GWL-NoCrops-B2040',
        'GWL-NoCrops-B2045',
        'GWL-NoCrops-B2050',
        'GWL-NoCrops-B2055',
        'GWL-NoCrops-B2060',
        'GWL-EqFor-B2060',
        ]
COLORS = {
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
        'PI-GWL-t6':np.arange(0*12, 700*12)/12,
        'PI-GWL-B2035':np.arange(0*12, 700*12)/12,
        'PI-GWL-B2040':np.arange(0*12, 700*12)/12,
        'PI-GWL-B2045':np.arange(0*12, 700*12)/12,
        'PI-GWL-B2050':np.arange(0*12, 700*12)/12,
        'PI-GWL-B2055':np.arange(0*12, 700*12)/12,
        'PI-GWL-B2060':np.arange(0*12, 700*12)/12,
        'GWL-NoCrops-B2030':np.arange(400*12, 601*12)/12,
        'GWL-NoCrops-B2035':np.arange(400*12, 502*12)/12,
        'GWL-NoCrops-B2040':np.arange(400*12, 502*12)/12,
        'GWL-NoCrops-B2045':np.arange(400*12, 502*12)/12,
        'GWL-NoCrops-B2050':np.arange(400*12, 502*12)/12,
        'GWL-NoCrops-B2055':np.arange(400*12, 502*12)/12,
        'GWL-NoCrops-B2060':np.arange(400*12, 502*12)/12,
        'GWL-EqFor-B2060':np.arange(400*12, 502*12)/12,
        }

load_nc = True

all_data = dict()
handles = list()
labels = list()
for exp in EXPERIMENTS:
    if load_nc:
        if 'NoCrops' in exp or 'EqFor' in exp:
            exp_dir = f'/g/data/p66/tfl561/ACCESS-ESM/{exp}/history/atm/netCDF'
        else:
            exp_dir = f'/g/data/p73/archive/non-CMIP/ACCESS-ESM1-5/{exp}/history/atm/netCDF'
        files = sorted(glob.glob(f'{exp_dir}/{exp}.pa-*_mon.nc'))[:700*12]
        sample_nc = nc.Dataset(files[0], 'r')
        lats = sample_nc.variables['lat'][:]
        all_data[exp] = np.ones((len(files),))*np.nan
        for i,f in enumerate(files):
            print(i, '/', len(files), '\r', end='')
            all_data[exp][i] = global_mean(load_co2(f, CO2_FLD), lats)*KGKG_TO_MOLMOL*MIL
        np.save(f'data/co2_surface_{exp}.npy', all_data[exp])
    else:
        all_data[exp] = np.load(f'data/co2_surface_{exp}.npy')

    if 'NoCrops' in exp or 'EqFor' in exp:
        lstyle = 'solid'
    else:
        lstyle = 'dotted'
    handle = plt.plot(
            YEARS[exp],
            all_data[exp],
            color=COLORS[exp],
            linestyle=lstyle,
            )
    if 'NoCrops' in exp or 'EqFor' in exp:
        handles.append(handle[0])
        labels.append(exp)

plt.legend(handles, labels, frameon=False)
plt.xlabel('Time (years)')
plt.ylabel('CO$_2$ (ppm)')
plt.xlim(left=0, right=700)
plt.savefig('plots/co2_surface_global_warming_level_no_crops.svg', dpi=200)
plt.show()

