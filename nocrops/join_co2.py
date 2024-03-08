#!/usr/bin/env python3
"""Join the .npy files created by get_co2.py into a single file.
"""
import glob
import numpy as np
import ipdb


EXPERIMENTS = [
        'GWL-10pct-B2030',
        'GWL-25pct-B2030',
        'GWL-50pct-B2030',
        'GWL-EGNL-B2030',
        'GWL-EGBL-B2030',
        'GWL-DCBL-B2030',
        'GWL-NoCrops-B2030',
        'GWL-NoCrops-B2035',
        'GWL-NoCrops-B2040',
        'GWL-NoCrops-B2045',
        'GWL-NoCrops-B2050',
        'GWL-NoCrops-B2055',
        'GWL-NoCrops-B2060',
        'GWL-EqFor-B2060',
        'PI-GWL-t6',
        'PI-GWL-B2035',
        'PI-GWL-B2040',
        'PI-GWL-B2045',
        'PI-GWL-B2050',
        'PI-GWL-B2055',
        'PI-GWL-B2060',
        ]

for exper in EXPERIMENTS:
    print(exper)
    files = sorted(glob.glob(f'data/co2_surface_{exper}_global_mean_*.npy'))
    data = np.array([])
    for f in files:
        data = np.concatenate([data, np.load(f)])
    np.save(f'data/co2_surface_{exper}_all.npy', data)
