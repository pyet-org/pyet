from numpy import pi, exp

from .utils import day_of_year, sunset_angle, solar_declination


def hamon(temperature, lat):
    """Evapotranspiration calculated with the Hamon (1961) method.

    Parameters
    ----------
    temperature: pandas.Series
        mean day temperature [°C]
    lat: float
        the site latitude [rad]

    Returns
    -------
    pandas.Series
        Series containing the calculated evapotranspiration

    Notes
    -----
    From Prudhomme(hess, 2013)

    """
    j = day_of_year(tindex=temperature.index)
    sol_dec = solar_declination(j=j)
    omega = sunset_angle(sol_dec=sol_dec, lat=lat)

    dl = 24 / pi * omega  # maximum possible daylight length

    return (dl / 12) ** 2 * exp(temperature / 16)
