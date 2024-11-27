#!/bin/env python3
"""The .npy files should be created with get_co2 and join_co2
"""
import glob

import ipdb
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt


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
        'PI-GWL-t6',
        'GWL-10pct-B2030',
        'GWL-25pct-B2030',
        #'GWL-50pct-B2030',
        #'GWL-EGNL-B2030',
        'GWL-EGBL-B2030',
        'GWL-DCBL-B2030',
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
        'GWL-10pct-B2030':'green',
        'GWL-25pct-B2030':'blue',
        'GWL-50pct-B2030':'yellow',
        'GWL-EGBL-B2030':'pink',
        'GWL-EGNL-B2030':'red',
        'GWL-DCBL-B2030':'purple',
        'GWL-NoCrops-B2030':'#62EA00',
        'GWL-NoCrops-B2035':'#24CC00',
        'GWL-NoCrops-B2040':'#079F2A',
        'GWL-NoCrops-B2045':'#00786B',
        'GWL-NoCrops-B2050':'#055992',
        'GWL-NoCrops-B2055':'#1140AB',
        'GWL-NoCrops-B2060':'#1E31B6',
        'GWL-EqFor-B2060':'deepskyblue'
        }
LABELS = {
        'GWL-10pct-B2030':'10% deployment',
        'GWL-25pct-B2030':'25% deployment',
        'GWL-50pct-B2030':'50% deployment',
        'GWL-EGBL-B2030':'Evergreen broad leaf',
        'GWL-EGNL-B2030':'Evergreen needle leaf',
        'GWL-DCBL-B2030':'Deciduous broad leaf',
        'GWL-NoCrops-B2030':'2030',
        'GWL-NoCrops-B2035':'2035',
        'GWL-NoCrops-B2040':'2040',
        'GWL-NoCrops-B2045':'2045',
        'GWL-NoCrops-B2050':'2050',
        'GWL-NoCrops-B2055':'2055',
        'GWL-NoCrops-B2060':'2060',
        'GWL-EqFor-B2060':'2030'
        }
YEARS = {
        'esm-esm-piNoCrops':np.arange(0*12, 100*12)/12 + 100,
        'PI-GWL-t6':np.arange(0*12, 700*12)/12,
        'PI-GWL-B2035':np.arange(0*12, 700*12)/12,
        'PI-GWL-B2040':np.arange(0*12, 700*12)/12,
        'PI-GWL-B2045':np.arange(0*12, 700*12)/12,
        'PI-GWL-B2050':np.arange(0*12, 700*12)/12,
        'PI-GWL-B2055':np.arange(0*12, 700*12)/12,
        'PI-GWL-B2060':np.arange(0*12, 700*12)/12,
        'GWL-10pct-B2030':np.arange(400*12, 601*12)/12,
        'GWL-25pct-B2030':np.arange(400*12, 601*12)/12,
        'GWL-50pct-B2030':np.arange(400*12, 601*12)/12,
        'GWL-EGNL-B2030':np.arange(400*12, 601*12)/12,
        'GWL-EGBL-B2030':np.arange(400*12, 601*12)/12,
        'GWL-DCBL-B2030':np.arange(400*12, 601*12)/12,
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
    if load_nc:
        if 'NoCrops' in exper or 'EqFor' in exper:
            exp_dir = f'/g/data/p66/tfl561/ACCESS-ESM/{exper}/history/atm/netCDF'
        else:
            exp_dir = f'/g/data/p73/archive/non-CMIP/ACCESS-ESM1-5/{exper}/'+\
                    'history/atm/netCDF'
        files = sorted(glob.glob(f'{exp_dir}/{exper}.pa-*_mon.nc'))[:700*12]
        sample_nc = nc.Dataset(files[0], 'r')
        lats = sample_nc.variables['lat'][:]
        all_data[exper] = np.ones((len(files),))*np.nan
        for i,f in enumerate(files):
            print(i, '/', len(files), '\r', end='')
            all_data[exper][i] = global_mean(load_co2(f, CO2_FLD), lats)
        np.save(f'data/co2_surface_{exper}_all.npy', all_data[exper])
    else:
        all_data[exper] = np.load(f'data/co2_surface_{exper}_all.npy')

    if 'NoCrops' in exper or 'EqFor' in exper:
        lwidth = 1.3
    else:
        lwidth = 0.7
    print(exper)
    x = yearly_mean(YEARS[exper])
    y = yearly_mean(all_data[exper])*KGKG_TO_MOLMOL*MIL
    handle = plt.plot(x, y,
            color=COLORS[exper],
            linewidth=lwidth,
            )
    if 'GWL' in exper[:4] in exper:
        handles.append(handle[0])
        labels.append(LABELS[exper])

plt.legend(handles, labels, fontsize='small', frameon=False, ncols=2)
plt.xlabel('Time (years)')
plt.ylabel('CO$_2$ (ppm)')
plt.xlim(left=0, right=700)
plt.savefig('plots/co2_surface_global_warming_level_no_crops.png', dpi=200)


def pretty_print(data):
    x = 0
    for i in range(data.size):
        if x==5:
            if 'int' in str(type(data[i])):
                print(f"{data[i]}, ", end='\n')
            else:
                print(f"{data[i]:e}, ", end='\n')
            x = 0
        else:
            if 'int' in str(type(data[i])):
                print(f"{data[i]}, ", end='')
            else:
                print(f"{data[i]:e}, ", end='')
            x += 1
    print("")


print("Namelist for GWL-NoCrops-B2030 co2")
print()

exper = 'GWL-NoCrops-B2030'
MISSINGVAL = -32768
nyears = int(YEARS[exper].size/12)
print("L_CLIM=.TRUE.,", end='\n')
print(f"CLIM_FCG_NYEARS={nyears},", end='\n')
print("CLIM_FCG_YEARS(:,1)=", end='')
pretty_print(yearly_mean(YEARS[exper]).astype(int))
print(f"CLIM_FCG_LEVELS(:,1)=", end='')
pretty_print(yearly_mean(all_data[exper]))
print(f"CLIM_FCG_RATES(:,1)=", end='')
pretty_print(np.ones(yearly_mean(YEARS[exper]).size, dtype=int)*MISSINGVAL)
