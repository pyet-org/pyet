import numpy as np
import pandas as pd


def hamon(temperature, latitude):
    """Returns evapotranspiration calculated with the Hamon (1961) method.
    Parameters
    ----------
    temperature: pandas.Series
        mean day temperature [°C]
    lat: float/int
        the site latitude [rad]
    Returns
    -------
        pandas.Series containing the calculated evapotranspiration
    Examples
    --------
    >>> ham = et.hamon(tmean, lat)
    """
    j = day_of_year(temperature.index)
    sol_dec = solar_declination(j)
    omega = sunset_angle(latitude, sol_dec)

    dl = 24 / np.pi * omega

    return (dl / 12) ** 2 * np.exp(temperature / 16)


def day_of_year(meteoindex):
    return pd.to_numeric(meteoindex.strftime('%j'))


def sunset_angle(lat, sol_dec):
    """
    Sunset hour angle [rad] (omega)
    From FAO (1990), ANNEX V, eq. 20
    """
    return np.arccos(-np.tan(lat) * np.tan(sol_dec))


def solar_declination(j):
    """
    Solar declination [rad] (sol_dec)
    From FAO (1990), ANNEX V, eq. 22
    """
    return 0.4093 * np.sin(2 * np.pi / 365 * j - 1.39)
