#!/usr/bin/env python3
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import cartopy
import cartopy.crs as ccrs
import netCDF4 as nc
import numpy as np
from cdo import Cdo
cdo = Cdo()
import cdo_decorators as cdod
import ipdb

PI_GWL_FILE = '/g/data/p73/archive/non-CMIP/ACCESS-ESM1-5/PI-GWL-B2060/history/atm/netCDF/PI-GWL-B2060.pa-010101_mon.nc'
READ_ONLY = 'r'

# Create land area map
land_area = cdo.gridarea(
        input=PI_GWL_FILE,
        returnCdf='data/land_area.nc',
        ).variables['cell_area'][:]

# Load
frac_GWL = nc.Dataset(PI_GWL_FILE, READ_ONLY).variables['fld_s03i317'][:].squeeze()

# Calculate tree areas.
trees_GWL = frac_GWL[8,...]*land_area

ncfile = nc.Dataset(PI_GWL_FILE, READ_ONLY)
lons = ncfile.variables['lon'][:]
lats = ncfile.variables['lat'][:]
ncfile.close()

# Plot data
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1, projection=ccrs.Robinson())
diff = trees_GWL
diff /= 1e10 # Mha
bnds = np.arange(0, 3.2, 0.2)
bnds[0] = 0.01
discrete_bins = mpl.colors.BoundaryNorm(
        boundaries=bnds,
        ncolors=256,
        )
cmap = plt.cm.get_cmap('YlGn')
cmap.set_under('white')
colors = ax.pcolormesh(
        lons,
        lats,
        diff,
        cmap=cmap,
        norm=discrete_bins,
        transform=ccrs.PlateCarree(),
        )
ax.coastlines()
plt.colorbar(
        colors,
        ticks=np.arange(0.1, 3.1, 0.2),
        label='Area (Mha)',
        orientation='horizontal',
        pad=0.05,
        )
plt.title('Tree area difference')
axins = inset_axes(
        ax,
        width="60%",
        height="60%",
        loc='lower left',
        bbox_to_anchor=(0.55,-0.05,0.9,0.9),
        bbox_transform=ax.transAxes,
        axes_class=cartopy.mpl.geoaxes.GeoAxes,
        axes_kwargs=dict(map_projection=ccrs.PlateCarree()),
        )
axins.add_feature(cartopy.feature.COASTLINE.with_scale('50m'))
axins.set_extent([110,155,-42,-10])
axins.pcolormesh(
        lons,
        lats,
        diff,
        cmap=cmap,
        norm=discrete_bins,
        transform=ccrs.PlateCarree(),
        )
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

#fig = plt.figure()
#plot_australia(diff, "Tree area difference")
#plt.tight_layout()
#plt.savefig('plots/tree_area_diff_australia.png', dpi=200)
plt.show()
