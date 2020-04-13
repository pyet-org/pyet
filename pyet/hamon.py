import numpy as np

from .utils import day_of_year, sunset_angle, solar_declination


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
