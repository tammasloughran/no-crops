#!/usr/bin/env python3
import matplotlib as mpl
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import netCDF4 as nc
import numpy as np
from cdo import Cdo
cdo = Cdo()
import cdo_decorators as cdod

PI_GWL_DIR = '/g/data/p73/archive/non-CMIP/ACCESS-ESM1-5'
GWL_EXPERIMENTS = [
        'PI-GWL-t6',
        'PI-GWL-B2035',
        'PI-GWL-B2040',
        'PI-GWL-B2045',
        'PI-GWL-B2050',
        'PI-GWL-B2055',
        'PI-GWL-B2060',
        'GWL-NoCrops-B2060',
        ]
DATA_DIR = 'history/atm/netCDF'


# Create list of PI-GWL files
files = [f'{PI_GWL_DIR}/{expe}/{DATA_DIR}/{expe}.pa-040101_mon.nc' for expe in GWL_EXPERIMENTS]
files[-1] = f'/g/data/p66/tfl561/ACCESS-ESM/GWL-NoCrops-B2060/history/atm/netCDF/'\
        'GWL-NoCrops-B2060.pa-050001_mon.nc'

# Create land area map
land_area = cdo.gridarea(input=files[0], returnCdf='data/land_area.nc').variables['cell_area'][:]

# Load PI-GWL maps
data = {}
trees = {}
for expe,f in zip(GWL_EXPERIMENTS, files):
    data[expe] = nc.Dataset(f, 'r').variables['fld_s03i317'][:].squeeze()
    trees[expe] = (data[expe][0:4,...]*land_area).sum(axis=0)

lons = nc.Dataset(f, 'r').variables['lon'][:]
lats = nc.Dataset(f, 'r').variables['lat'][:]

for expe in GWL_EXPERIMENTS:
    print(expe, trees[expe].sum()*10**-12)

# Plot data
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1, projection=ccrs.Robinson())
diff = trees['GWL-NoCrops-B2060'] - trees['PI-GWL-B2060']
diff /= 1e10 # Mha
bnds = np.arange(0, 3.2, 0.2)
bnds[0] = 0.01
discrete_bins = mpl.colors.BoundaryNorm(
        boundaries=bnds,
        ncolors=256,
        )
cmap = plt.cm.get_cmap('YlGn')
cmap.set_under('white')
colors = ax.pcolormesh(lons, lats, diff,
        cmap=cmap,
        edgecolors='face',
        linewidth=0.2,
        norm=discrete_bins,
        transform=ccrs.PlateCarree(),
        )
ax.coastlines()
plt.colorbar(colors,
        ticks=np.arange(0.1, 3.1, 0.2),
        label='Area (Mha)',
        orientation='horizontal',
        pad=0.05,
        )
plt.title('Tree area difference')
plt.savefig(f'plots/tree_area_diff.png', dpi=200)
plt.show()

def plot_australia(data:np.ndarray, title:str)->None:
    """Plot a map of Asutralia anomalies.
    """
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.set_extent([110,160,-45,-10], crs=ccrs.PlateCarree())
    bnds = np.arange(0, 3.2, 0.2)
    bnds[0] = 0.01
    discrete_bins = mpl.colors.BoundaryNorm(
            boundaries=bnds,
            ncolors=256
            )
    cmap = plt.cm.get_cmap('YlGn')
    cmap.set_under('white')
    shading = ax.pcolormesh(lons, lats, data,
            cmap=cmap,
            edgecolors='face',
            linewidth=0.2,
            norm=discrete_bins,
            transform=ccrs.PlateCarree(),
            )
    cbar = plt.colorbar(shading,
            ticks=np.arange(0.1, 3.1, 0.2),
            label='Area (Mha)',
            orientation='horizontal',
            pad=0.05,
            )
    ax.coastlines()
    plt.title(title)

fig = plt.figure()
plot_australia(diff, "Tree area difference")
plt.tight_layout()
plt.savefig('plots/tree_area_diff_australia.png', dpi=200)
plt.show()
