"""The temeprature module contains functions of temeprature PET methods

"""

from numpy import exp

from .meteo_utils import daylight_hours, calc_ea, calc_es


def blaney_criddle(tmean, p, k=0.85):
    """Evaporation calculated according to [blaney_1952]_.

    Parameters
    ----------
    tmean: pandas.Series, optional
        average day temperature [°C]
    p: pandas.Series/float, optional
        bright sunshine (hour day-1)
    k: float, optional
        calibration coefficient [-]

    Returns
    -------
    pandas.Series containing the calculated evaporation.

    Examples
    --------
    >>> et_blaney_criddle = blaney_criddle(tmean)

    Notes
    -----
    Based on equation 6 in [xu_2001]_.

    .. math:: PE=kp(0.46 * T_a + 8.13)

    References
    ----------
    .. [blaney_1952] Blaney, H. F. (1952). Determining water requirements in
       irrigated areas from climatological and irrigation data.
    .. [xu_2001] Xu, C. Y., & Singh, V. P. (2001). Evaluation and
       generalization of temperature‐based methods for calculating evaporation.
       Hydrological processes, 15(2), 305-319.
    """
    et = k * p * (0.46 * tmean + 8.13)
    return et


def hamon(tmean, lat):
    """Evaporation calculated according to [hamon_1961]_.

    Parameters
    ----------
    tmean: pandas.Series, optional
        average day temperature [°C]
    lat: float, optional
        the site latitude [rad]

    Returns
    -------
    pandas.Series containing the calculated evaporation.

    Examples
    --------
    >>> et_hamon = hamon(tmean, lat)

    Notes
    -----
    Following [hamon_1961]_ and [oudin_2005]_.

    .. math:: PE = (\\frac{DL}{12})^2 exp(\\frac{T_a}{16})

    References
    ----------
    .. [hamon_1961] Hamon, W. R. (1963). Estimating potential
       evapotranspiration. Transactions of the American Society of Civil
       Engineers, 128(1), 324-338.
    .. [oudin_2005] Oudin, L., Hervieu, F., Michel, C., Perrin, C.,
       Andréassian, V., Anctil, F., & Loumagne, C. (2005). Which potential
       evapotranspiration input for a lumped rainfall–runoff model?:
       Part 2—Towards a simple and efficient potential evapotranspiration model
       for rainfall–runoff modelling. Journal of hydrology, 303(1-4), 290-306.

    """

    dl = daylight_hours(tmean.index, lat)

    return (dl / 12) ** 2 * exp(tmean / 16)


def romanenko(tmean, rh, k=4.5):
    """Evaporation calculated according to [romanenko_1961]_.

    Parameters
    ----------
    tmean: pandas.Series, optional
        average day temperature [°C]
    rh: pandas.Series, optional
        mean daily relative humidity [%]
    k: float, optional
        calibration coefficient [-]

    Returns
    -------
    pandas.Series containing the calculated evaporation.

    Examples
    --------
    >>> et_romanenko = romanenko(tmean, rh)

    Notes
    -----
    Based on equation 11 in [xu_2001]_.

    .. math:: PE=4.5(1 + (\\frac{T_a}{25})^2 (1  \\frac{e_a}{e_s})

    References
    ----------
    .. [romanenko_1961] Romanenko, V. A. (1961). Computation of the autumn soil
       moisture using a universal relationship for a large area. Proc. of
       Ukrainian Hydrometeorological Research Institute, 3, 12-25.
    """
    ea = calc_ea(tmean=tmean, rh=rh)
    es = calc_es(tmean=tmean)

    return k * (1 + tmean / 25) ** 2 * (1 - ea / es)


def linacre(tmean, elevation, lat, tdew=None, tmax=None, tmin=None):
    """Evaporation calculated according to [linacre_1977]_.

    Parameters
    ----------
    tmean: pandas.Series, optional
        average day temperature [°C]
    elevation: float, optional
        the site elevation [m]
    lat: float, optional
        the site latitude [°]
    tdew: pandas.Series, optional
        mean dew-point temperature [°C]
    tmax: pandas.Series, optional
        maximum day temperature [°C]
    tmin: pandas.Series, optional
        minimum day temperature [°C]

    Returns
    -------
    pandas.Series containing the calculated evaporation.

    Examples
    --------
    >>> et_linacre = linacre(tmean, elevation, lat)

    Notes
    -----
    Based on equation 5 in [xu_2001]_.

    .. math:: PE = \\frac{\\frac{500 T_m}{(100-A)}+15 (T_a-T_d)}{80-T_a}

    References
    -----
    .. [linacre_1977] Linacre, E. T. (1977). A simple formula for estimating
       evaporation rates in various climates, using temperature data alone.
       Agricultural meteorology, 18(6), 409-424.
    """
    if tdew is None:
        tdew = 0.52 * tmin + 0.6 * tmax - 0.009 * tmax ** 2 - 2
    tm = tmean + 0.006 * elevation
    et = (500 * tm / (100 - lat) + 15 * (tmean - tdew)) / (80 - tmean)
    return et
