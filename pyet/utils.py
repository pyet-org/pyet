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


def check_lat(lat):
    """Method to check if latitude was (most likely) given in radians."""
    if lat is None:
        pass
    else:
        if -1.6 < numpy.mean(lat) < 1.6:
            return lat
        else:
            print("Latitude must be provided in radians!")
            return lat


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
