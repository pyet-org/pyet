from numpy import tan, pi, sin, exp, arccos
from pandas import to_numeric


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

    dl = 24 / pi * omega

    return (dl / 12) ** 2 * exp(temperature / 16)


def day_of_year(meteoindex):
    return to_numeric(meteoindex.strftime('%j'))


def sunset_angle(lat, sol_dec):
    """
    Sunset hour angle [rad] (omega)
    From FAO (1990), ANNEX V, eq. 20
    """
    return arccos(-tan(lat) * tan(sol_dec))


def solar_declination(j):
    """
    Solar declination [rad] (sol_dec)
    From FAO (1990), ANNEX V, eq. 22
    """
    return 0.4093 * sin(2 * pi / 365 * j - 1.39)
