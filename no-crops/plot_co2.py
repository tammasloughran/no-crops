#!/usr/env python3
import glob

import numpy as np
import matplotlib.pyplot as plt

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
        }

all_data = dict()
handles = list()
labels = list()
for experiment in EXPERIMENTS:
    files = sorted(glob.glob(f'data/co2_{experiment}_global_mean_?.npy'))

    all_data[experiment] = np.empty((0))
    for f in files:
        all_data[experiment] = np.concatenate((all_data[experiment],np.load(f)))

    np.save(f'data/co2_{experiment}_all.npy', all_data[experiment])

    if 'NoCrops' in experiment:
        lstyle = 'dotted'
    else:
        lstyle = 'solid'
    handle = plt.plot(
            YEARS[experiment],
            all_data[experiment],
            color=COLORS[experiment],
            linestyle=lstyle,
            )
    if 'NoCrops' in experiment:
        handles.append(handle[0])
        labels.append(experiment)

plt.legend(handles, labels, frameon=False)
plt.xlabel('Time (years)')
plt.ylabel('CO$_2$ (ppm)')
plt.xlim(left=0, right=700)
plt.savefig('co2_global_warming_level_no_crops.svg', dpi=200)
plt.show()

