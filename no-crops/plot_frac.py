import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import netCDF4 as nc
import numpy as np
from cdo import Cdo; cdo = Cdo()
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
files[-1] = f'/g/data/p66/tfl561/ACCESS-ESM/GWL-NoCrops-B2060/history/atm/netCDF/GWL-NoCrops-B2060.pa-050001_mon.nc'

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

# define Nocrop versions

# Load NoCrops versions

# Plot data
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1, projection=ccrs.Robinson())
diff = trees['GWL-NoCrops-B2060'] - trees['PI-GWL-B2060']
diff /= 1e12
colors = ax.pcolormesh(lons, lats, diff,
        vmin=0,
        vmax=0.03,
        cmap='YlGn',
        transform=ccrs.PlateCarree(),
        )
ax.coastlines()
plt.colorbar(colors,
        label='Area Mkm$^2$',
        orientation='horizontal',
        pad=0.05,
        )
plt.title('Tree area difference')
plt.savefig(f'plots/tree_area_diff.png', dpi=200)
plt.show()
