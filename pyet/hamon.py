from numpy import pi, exp

from .utils import day_of_year, sunset_angle, solar_declination


def hamon(tindex, tmean, lat):
    """Evaporation calculated according to [hamon_1961]_.

    Parameters
    ----------
    tindex: pandas.Series.index
    tmean: pandas.Series, optional
        average day temperature [°C]
    lat: float, optional
        the site latitude [rad]

    Returns
    -------
    pandas.Series containing the calculated evaporation.

    Examples
    --------
    >>> et_ham = hamon(tmean, lat)

    .. math::
    -----
        E = (\\frac{DL}{12})^2 exp(\\frac{T_a}{16})

    References
    -----
    .. [hamon_1961] Hamon, W. R. (1963). Estimating potential
       evapotranspiration. Transactions of the American Society of Civil
       Engineers, 128(1), 324-338.
    .. [oudin_2005] Oudin, L., Hervieu, F., Michel, C., Perrin, C.,
       Andréassian, V., Anctil, F., & Loumagne, C. (2005). Which potential
       evapotranspiration input for a lumped rainfall–runoff model?:
       Part 2—Towards a simple and efficient potential evapotranspiration model
       for rainfall–runoff modelling. Journal of hydrology, 303(1-4), 290-306.

    """
    j = day_of_year(tindex=tindex)
    sol_dec = solar_declination(j=j)
    omega = sunset_angle(sol_dec=sol_dec, lat=lat)

    dl = 24 / pi * omega  # maximum possible daylight length

    return (dl / 12) ** 2 * exp(tmean / 16)
