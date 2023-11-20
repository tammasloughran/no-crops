#!/usr/bin/env python3
# Plot the carbon pools from the single pft no-crops forestation experiments.
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc
import glob
import sys
[sys.path.append(i) for i in ['.', '..']]
from contrib.tools import global_mean

EXPERIMENTS = [
        'PI-GWL-t6',
        'GWL-NoCrops-B2030',
        'GWL-EGNL-B2030',
        'GWL-EGBL-B2030',
        'GWL-DCBL-B2030',
        ]
ARCHIVE_DIR = '/g/data/p66/tfl561/archive_data'
WOODFIG = 1

example_file = f'{ARCHIVE_DIR}/GWL-NoCrops-B2030/cWood_GWL-NoCrops-B2030_0500-0700.nc'
lats = nc.Dataset(example_file, 'r').variables['lat'][:]

plt.figure(WOODFIG)

data = {}
for exper in EXPERIMENTS:
    data[exper] = {}
    var = 'cWood'
    files = glob.glob(f'{ARCHIVE_DIR}/{exper}/{var}_{exper}_*')
    print("Loading: ", files[0])
    ncfile = nc.Dataset(files[0], 'r')
    data[exper][var] = global_mean(np.sum(ncfile.variables[var][:,0:3,...], axis=1), lats)
    np.save(f'data/{var}_{exper}_foo.npy', data[exper][var].data)

    # Plot
    plt.plot(range(len(data[exper][var])), data[exper][var], label=exper)

plt.legend()
plt.show()
