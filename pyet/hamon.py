import numpy as np
import pandas as pd


def hamon(temperature, latitude):
    """Returns evapotranspiration calculated with the Hamon (1961) method.
    From Prudhomme(hess, 2013)
    Parameters
    ----------
    temperature: Series
        mean day temperature [°C]
    latitude: float/int
        the site latitude [rad]
    Returns
    -------
        Series containing the calculated evapotranspiration
    """
    j = day_of_year(temperature.index)
    sol_dec = solar_declination(j)
    omega = sunset_angle(latitude, sol_dec)

    dl = 24 / np.pi * omega  # maximum possible daylight length

    return (dl / 12) ** 2 * np.exp(temperature / 16)


def day_of_year(meteoindex):
    """
    Return day of the year (1-365) based on pandas.series.index

    Parameters
    ----------
    meteoindex: pandas.Series.index

    Returns
    -------
        array of with ints specifyng day of year.

    """
    return pd.to_numeric(meteoindex.strftime('%j'))


def relative_distance(j):
    """
    Inverse relative distance between earth and sun from day of the year.

    Based on equation 23 in Allen et al (1998).

    Parameters
    ----------
    j: array.py
        day of the year (1-365)
    Returns
    -------
        array.py specifyng day of year.
    """
    return 1 + 0.033 * np.cos(2 * np.pi / 365 * j)


def sunset_angle(sol_dec, lat):
    """
    Sunset hour angle from latitude and solar declination [rad].

    Based on equations 25 in ALLen et al (1998).
    Parameters
    ----------
    sol_dec: pandas.Series
        solar declination [rad]
    lat: float/int
        the site latitude [rad]
    Returns
    -------
        pandas.Series of sunset hour angle [rad].

    """
    return np.arccos(-np.tan(lat) * np.tan(sol_dec))


def solar_declination(j):
    """
    Solar declination [rad] from day of year [rad].

    Based on equations 24 in ALLen et al (1998).
    Parameters
    ----------
    j: array.py
        day of the year (1-365)
    Returns
    -------
        array.py of solar declination [rad].

    """
    return 0.4093 * np.sin(2 * np.pi / 365 * j - 1.39)
