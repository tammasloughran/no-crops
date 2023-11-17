import numpy as _np
import matplotlib.pyplot as _plt
import matplotlib as _mpl
import cartopy.crs as _ccrs
import netCDF4 as _nc
import datetime as _dt


def yearly_mean_from_monthly(data:_np.ndarray)->_np.ndarray:
    """Calculate a yearly mean on a numpy array of monthly data.
    The 0th dimension must be time and divisible by 12.
    """
    if _np.mod(data.shape[0], 12)!=0:
        raise ValueError("Not enough months in 0th dimension.")
    toshape = list(data.shape)
    toshape.pop(0)
    toshape.insert(0, 12)
    toshape.insert(0, int(data.shape[0]/12))
    fraction_of_year = _np.array([31,28,31,30,31,30,31,31,30,31,30,31])/365.0
    return _np.average(data.reshape(toshape), axis=1, weights=fraction_of_year)


def global_mean(data:_np.ndarray, lats:_np.ndarray)->_np.ndarray:
    """Calculate an area weighted global mean. Weights are the cosine of lats.
    """
    if data.ndim<2: raise ValueError("No lats or lons.")
    if len(lats)!=data.shape[-2]:
        raise ValueError("size of lats does not match the provided data.")
    coslats = _np.cos(_np.deg2rad(lats))[:,None]*_np.ones(data.shape)
    return _np.ma.average(data, axis=(-1,-2), weights=coslats)


def num2date(times:_np.ndarray, units:str, calendar:str)->list:
    """Convert numbers to a list of datetime objects.
    """
    dates = _nc.num2date(times, units, calendar=calendar)
    return [_dt.datetime(d.year, d.month, d.day, d.hour) for d in dates]


def plot_globe_diverging_discrete(
            lats:_np.ndarray,
            lons:_np.ndarray,
            data:_np.ndarray,
            title:str,
            filename:str,
            cmap:str='bwr',
            )->None:
    """Plot a global map of some data using the Robinson projection and a
    discrete colormap.
    """
    _plt.figure()
    ax = _plt.axes(projection=_ccrs.Robinson())
    absmax = max(abs(_np.nanmin(data)), _np.nanmax(data))
    discrete_bins = _mpl.colors.BoundaryNorm(
            boundaries=_np.arange(-absmax, absmax+1, 1),
            ncolors=256,
            )
    _plt.pcolormesh(lons, lats, data,
            norm=discrete_bins,
            cmap=cmap,
            linewidth=0.2,
            edgecolors='face',
            transform=_ccrs.PlateCarree(),
            )
    ax.coastlines()
    _plt.colorbar(
            ticks=_np.arange(-absmax, absmax, 1),
            orientation='horizontal',
            pad=0.05,
            )
    _plt.title(title)
    _plt.tight_layout()
    _plt.savefig(filename, dpi=200)

