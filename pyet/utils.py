import numpy
import pandas


def show_versions():
    """Method to print the version of dependencies.

    """
    from pyet import __version__ as ps_version
    from pandas import __version__ as pd_version
    from numpy import __version__ as np_version
    from sys import version as os_version

    msg = (
        f"Python version: {os_version}\n"
        f"Numpy version: {np_version}\n"
        f"Pandas version: {pd_version}\n"
        f"Pyet version: {ps_version}"
    )
    return print(msg)


def deg_to_rad(lat):
    """Method to convert latitude in degrees to radians.
    lat: float/xarray.DataArray
        the site latitude [deg]

    Returns
    -------
    float/pandas.Series/xarray.DataArray containing the calculated
            latitude in radians [rad].
    """
    return lat * numpy.pi / 180


def check_rad(rad):
    """Method to check if radiation was provided in MJ/m2d."""
    if rad is not None:
        if numpy.max(rad) < 100:
            return rad
        else:
            raise Exception(
                "The radiation input provided is greater than 100 MJ/m2d, "
                "which is not realistic. Please convert the radiation input"
                " to MJ/m2d.")


def check_rh(rh):
    """Method to check if relative humidity is provided in percentage."""
    if rh is not None:
        if numpy.max(rh) > 1.0:
            return rh
        else:
            raise Exception(
                "The maximum value of relative humidity provided is smaller "
                "than 1 [%], which is not realistic. Please convert the "
                "relative humidity to [%].")
    else:
        pass


def check_lat(lat):
    """Method to check if latitude was (most likely) given in radians."""
    if lat is None:
        pass
    else:
        if -1.6 < numpy.mean(lat) < 1.6:
            return lat
        else:
            raise Exception(
                "Latitude must be provided in radians! Use pyet.deg_to_rad()"
                "to convert from degrees to radians.")


def clip_zeros(s, clip_zero):
    """Method to replace negative values with 0 for Pandas.Series and
        xarray.DataArray.

    """
    if clip_zero:
        return s.where((s > 0) | (s.isnull()), 0)
    else:
        return s


def get_index(df):
    """Method to return the index of the input data.

    """
    try:
        index = pandas.DatetimeIndex(df.index)
    except AttributeError:
        index = pandas.DatetimeIndex(df.time)
    return index
