import numpy as np
import pandas as pd


def hargreaves(tmax, tmin, lat):
    """Returns evapotranspiration calculated with the Hargreaves and
    Samani (1975) method.
    Parameters
    ----------
    tmax: pandas.Series
        maximum day temperature [°C]
    tmin: pandas.Series
        minimum day temperature [°C]
    lat: float/int
        the site latitude [rad]
    Returns
    -------
        pandas.Series containing the calculated evapotranspiration
    Examples
    --------
    >>> har = et.har(tmax, tmin, lat)
    """
    ta = (tmax + tmin) / 2
    lambd = lambda_calc(ta)
    ra = ra_calc(ta.index, lat)
    return 0.0023 * ra * (ta + 17.8) * np.sqrt(tmax - tmin) / lambd


def lambda_calc(temperature):
    """
    From FAO (1990), ANNEX V, eq. 1
    """
    return 2.501 - 0.002361 * temperature


def ra_calc(meteoindex, lat):
    """
    Extraterrestrial Radiation (Ra)
    From FAO (1990), ANNEX V, eq. 18
    """

    j = day_of_year(meteoindex)
    dr = relative_distance(j)
    sol_dec = solar_declination(j)

    omega = sunset_angle(lat, sol_dec)
    gsc = 0.082 * 24 * 60
    # gsc = 1360
    return gsc / 3.141592654 * dr * (omega * np.sin(sol_dec) * np.sin(lat) +
                                     np.cos(sol_dec) * np.cos(lat) * np.sin(
                omega))


def relative_distance(j):
    """
    Relative distance Earth - Sun
    From FAO (1990), ANNEX V, eq. 21
    """
    return 1 + 0.033 * np.cos(2 * np.pi / 365 * j)


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
