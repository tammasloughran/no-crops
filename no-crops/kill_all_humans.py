#!/usr/bin/env python3
# kill_all_humans.py removes all crops from CABLE land-use data and replaces them with trees.
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt


def euclidean(x1, y1, x2, y2):
    """Euclidean distance between points (x1, y1) and (x2, y2).
    """
    return np.sqrt((x1 - x2)**2 + (y1 - y2)**2)


class MyDict(dict):
    """A dictionary with a method `get_key()` to search it's values and return it's key.
    """


    def get_key(self, value):
        """Get a dictionary key for a given value.
        """
        return list(self.keys())[list(self.values()).index(value)]


# Constants
CABLE_TILES = MyDict({
        1:'Evergreen needleleaf',
        2:'Evergreen broadleaf',
        3:'Deciduous needleleaf',
        4:'Deciduous broadleaf',
        5:'Shrub',
        6:'C3 grass',
        7:'C4 grass',
        8:'Tundra',
        9:'C3 crop',
        10:'C4 crop',
        11:'Wetland',
        12:'Not used',
        13:'Not used',
        14:'barren',
        15:'urban',
        16:'lakes',
        17:'ice',
        })
CROP_INDEX = CABLE_TILES.get_key('C3 crop') - 1

# Load land-use map.
year = 1850
a = 'txz'; b = '599'; uname = a + b # This is to avoid putting a complete username in code.
ncin = nc.Dataset(f'/g/data/p66/{uname}/data/luc_hist_maps/cableCMIP6_LC_{year}.nc', 'r')
lats = ncin.variables['latitude'][:]
lons = ncin.variables['longitude'][:]
xx, yy = np.meshgrid(lons, lats)
fraction = ncin.variables['fraction'][:].squeeze()
crops = fraction[CROP_INDEX]

# Create map of all 0's on land gridpoints.
zeros = crops.copy()
zeros[crops<2] = 0

# Create fractions map without crops.
no_humans = fraction.copy()
no_humans[CROP_INDEX,...] = zeros

# Add tree plant functional types to fill in the missing crops, preserving the proportions of
# tree fractions.
TREE_INDEX = [CABLE_TILES.get_key(val) - 1 for val in [
        'Evergreen needleleaf',
        'Evergreen broadleaf',
        'Deciduous needleleaf',
        'Deciduous broadleaf',
        ]]
tree_proportions = fraction[TREE_INDEX,...]/fraction[TREE_INDEX].sum(axis=0)
# Rmove missing values from above divide by 0's.
no_trees = fraction[TREE_INDEX].sum(axis=0)==0
for t in TREE_INDEX:
    tree_proportions[t,no_trees] = 0
no_humans[TREE_INDEX,...] = fraction[TREE_INDEX] + crops*tree_proportions

# Some gridpoints have crops but no trees.
has_trees = fraction[TREE_INDEX].sum(axis=0)>0
has_crops_but_no_trees = np.logical_and(
        fraction[CROP_INDEX]>0,
        no_trees,
        )
# Put trees there 1/3 evergreen needle, 1/3 evergreen broad and 1/3 deciduous broadleaf.
# Deciduous needle stands are rare. This is the quick and dirty way. See a better way below.
#COMMON_TREES = [CABLE_TILES.get_key(val) - 1 for val in [
#        'Evergreen needleleaf',
#        'Evergreen broadleaf',
#        'Deciduous broadleaf',
#        ]]
#for t in COMMON_TREES:
#    no_humans[t,has_crops_but_no_trees] = crops[has_crops_but_no_trees]/3.0

# Or even better: find the nearest gridpoint with trees and use the tree proportions from there.
for j in range(len(lats)):
    for i in range(len(lons)):
        if has_crops_but_no_trees[j,i]:
            distance = euclidean(xx, yy, lons[i], lats[j])
            distance[np.logical_not(has_trees)] = np.nan
            neary, nearx = np.unravel_index(np.nanargmin(distance), distance.shape)
            no_humans[TREE_INDEX,j,i] = fraction[TREE_INDEX,j,i] + \
                    crops[j,i]*tree_proportions[:,neary,nearx]

# Make sure there are no fractions less than 1.0E-8
for t in TREE_INDEX:
    invalid_values = (no_humans[t,...]<1.0E-8)&(no_humans[t,...]>0.0)
    if invalid_values.any():
        print("There are tiles less than the threshold in ", CABLE_TILES[t+1])
        #plt.figure()
        #plt.pcolormesh(invalid_values)
        #plt.show()
        yy, xx = np.where(invalid_values==True)
        for i in range(len(xx)):
            dominant_tree = no_humans[0:4,yy[i],xx[i]].argmax()
            no_humans[dominant_tree,yy[i],xx[i]] = fraction[CROP_INDEX,yy[i],xx[i]]
            for j in TREE_INDEX:
                if j!=dominant_tree: no_humans[j,yy[i],xx[i]] = 0.0
            print(xx[i], yy[i], "altered")
        # Debug info
        print("Crops")
        for i in range(len(xx)):
            print(xx[i], yy[i], fraction[CROP_INDEX, yy[i], xx[i]])
        print("trees")
        for tree in TREE_INDEX:
            for i in range(len(xx)):
                print(xx[i], yy[i], fraction[tree, yy[i], xx[i]])
        print("after transition")
        print("Crops")
        for i in range(len(xx)):
            print(xx[i], yy[i], no_humans[CROP_INDEX, yy[i], xx[i]])
        print("trees")
        for tree in TREE_INDEX:
            for i in range(len(xx)):
                print(xx[i], yy[i], no_humans[tree, yy[i], xx[i]])

# Check that tiles add to 1.
plt.figure()
plt.pcolormesh(no_humans.sum(axis=0))
plt.colorbar()
plt.title('Sum of tiles')

for t in TREE_INDEX:
    plt.figure()
    plt.pcolormesh(no_humans[t,...] - fraction[t,...])
    plt.colorbar()
    plt.title(f'{CABLE_TILES[t+1]} difference')

plt.show()

# Save to netcdf file.
ncout = nc.Dataset(f'{year}_no_humans_CABLE_fraction.nc', 'w')
for d in ncin.dimensions:
    ncout.createDimension(d, ncin.dimensions[d].size)
outvars = {}
for var in ncin.variables:
    fill = None
    try:
        if ncin.variables[var]._FillValue: fill = ncin.variables[var]._FillValue
    except:
        pass
    outvars[var] = ncout.createVariable(
            var,
            ncin.variables[var].dtype,
            ncin.variables[var].dimensions,
            fill_value=fill,
            )
outvars['longitude'].units = 'degrees_east'
outvars['longitude'].long_name = 'longitude'
outvars['latitude'].units = 'degrees_north'
outvars['latitude'].long_name = 'latitude'
outvars['vegtype'].units = '\n'.join(list(CABLE_TILES.values()))
#for att in ncin.variables['fraction'].ncattrs():
#    if att=='_FillValue': continue
#    setattr(outvars['fraction'], att, ncin.variables['fraction'].getncattr(att))
for att in ncin.variables['time'].ncattrs():
    setattr(outvars['time'], att, ncin.variables['time'].getncattr(att))
for var in ['longitude','latitude','vegtype','time']:
    outvars[var][:] = ncin.variables[var][:]
outvars['fraction'][:] = no_humans
ncout.note = "Created by /home/561/tfl561/sources/sensitivity_lu_map/kill_all_humans.py"
ncout.close()
ncin.close()

