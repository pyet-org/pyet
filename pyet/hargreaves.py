from numpy import sqrt

from .utils import extraterrestrial_r


def hargreaves(tmax, tmin, lat):
    """Evapotranspiration calculated with Hargreaves and Samani (1975) method.

    Parameters
    ----------
    tmax: pandas.Series
        maximum day temperature [°C]
    tmin: pandas.Series
        minimum day temperature [°C]
    lat: float
        the site latitude [rad]

    Returns
    -------
    pandas.Series
        Series containing the calculated evapotranspiration.

    Examples
    --------
    >>> har_et = hargreaves(tmax, tmin, lat)

    """
    ta = (tmax + tmin) / 2
    ra = extraterrestrial_r(tindex=tmax.index, lat=lat)
    return 0.408 * 0.0023 * ra * (ta + 17.8) * sqrt(tmax - tmin)
