import numpy as np

from .utils import extraterrestrial_r


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
    >>> har_et = hargreaves(tmax, tmin, lat)
    """
    ta = (tmax + tmin) / 2
    ra = extraterrestrial_r(tmax.index, lat)
    return 0.408 * 0.0023 * ra * (ta + 17.8) * np.sqrt(tmax - tmin)
