import numpy
from pandas import Series, DatetimeIndex
from xarray import DataArray


def show_versions():
    """Method to print the version of dependencies."""
    from pyet import __version__ as ps_version
    from pandas import __version__ as pd_version
    from numpy import __version__ as np_version
    from sys import version as os_version
    from xarray import __version__ as xr_version

    msg = (
        f"Python version: {os_version}\n"
        f"Numpy version: {np_version}\n"
        f"Pandas version: {pd_version}\n"
        f"xarray version: {xr_version}\n"
        f"Pyet version: {ps_version}"
    )
    return print(msg)


def deg_to_rad(lat):
    """Method to convert latitude in degrees to radians.

    Parameters
    ----------
    lat: float or xarray.DataArray
        The site latitude [deg].

    Returns
    -------
    float or pandas.Series or xarray.DataArray containing the calculated latitude in
    radians [rad].
    """
    return lat * numpy.pi / 180


def check_rad(rad):
    """Method to check if radiation was probably provided in MJ/m2d."""
    if rad is not None:
        if numpy.nanmax(rad) < 100:
            return rad
        else:
            raise Exception(
                "The radiation input provided is greater than 100 MJ/m2d, "
                "which is not realistic. Please convert the radiation input"
                " to MJ/m2d."
            )


def check_rh(rh):
    """Method to check if relative humidity is provided in percentage."""
    if rh is not None:
        if numpy.nanmax(rh) > 1.0:
            return rh
        else:
            raise Exception(
                "The maximum value of relative humidity provided is smaller "
                "than 1 [%], which is not realistic. Please convert the "
                "relative humidity to [%]."
            )
    else:
        pass


def check_lat(lat, shape=None):
    """Method to check if latitude was (most likely) given in radians."""
    if not isinstance(lat, (float, int, DataArray)):
        raise TypeError("lat must be a float, int or DataArray")
    if isinstance(lat, (float, int)) and (shape is None or len(shape) == 1):
        lat1 = lat
    else:
        if isinstance(lat, (float, int)) and (len(shape) > 1):
            raise ValueError(f"lat must be a shaped as 2D DataArray")
        lat1 = lat.values
    if not (-1.6 < numpy.mean(lat1) < 1.6):
        raise Exception(
            "Latitude must be provided in radians! Use pyet.deg_to_rad()"
            "to convert from degrees to radians."
        )
    return lat1


def clip_zeros(s, clip_zero):
    """Method to replace negative values with 0 for Pandas.Series and xarray.DataArray."""
    if clip_zero:
        s = s.where((s >= 0) | s.isnull(), 0)
    return s


def pet_out(tmean, pet, name):
    """Method to create pandas.Series or xarray.DataArray from numpy.ndarray"""
    if isinstance(tmean, (Series, DataArray)):
        return pet.rename(name)
    else:
        raise TypeError("Input must be either pandas.Series or xarray.DataArray!")


def get_index(df):
    """Method to return the index of the input data."""
    try:
        index = DatetimeIndex(df.index)
    except AttributeError:
        index = DatetimeIndex(df.time)
    return index


def vectorize(*arrays):
    """Vectorize pandas.Series or xarray.DataArray inputs."""
    vec_arrays = []
    for arr in arrays:
        if arr is None:
            vec_arr = None
        elif isinstance(arr, (int, float)):
            vec_arr = arr
        elif isinstance(arr, Series):
            vec_arr = arr.copy().values
        elif isinstance(arr, DataArray):
            vec_arr = arr.copy().values
        else:
            raise TypeError(
                f"Input must be a pandas.Series or xarray.DataArray, "
                f"but got {type(arr)}"
            )
        vec_arrays.append(vec_arr)
    return vec_arrays
