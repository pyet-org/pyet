"""The temeprature module contains functions of temeprature PET methods

"""

from numpy import exp, broadcast_to

from .meteo_utils import daylight_hours, calc_ea, calc_es, calc_e0

from .utils import get_index_shape


def blaney_criddle(tmean, lat, k=0.65):
    """Evaporation calculated according to [blaney_1952]_.

    Parameters
    ----------
    tmean: pandas.Series, optional
        average day temperature [°C]
    lat: float, optional
        the site latitude [rad]
    k: float, optional
        calibration coefficient [-]

    Returns
    -------
    pandas.Series containing the calculated evaporation.

    Examples
    --------
    >>> et_blaney_criddle = blaney_criddle(tmean, lat)

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
    index, shape = get_index_shape(tmean)
    dl = broadcast_to(daylight_hours(index, lat, shape), shape)
    et = k * dl / (365 * 12) * 100 * (0.46 * tmean + 8.13)
    return et


def haude(tmean, rh, k=1):
    """Evaporation calculated according to [haude]_.

    Parameters
    ----------
    tmean: pandas.Series, optional
        temperature at 2pm or maximum dailty temperature [°C]
    rh: float, optional
        average relative humidity at 2pm [%]
    k: float, optional
        calibration coefficient [-]

    Returns
    -------
    pandas.Series containing the calculated evaporation.

    Examples
    --------
    >>> et_haude = haude(tmean, rh)

    Notes
    -----
    Following [haude_1955]_ and [schiff_1975]_.

    .. math:: PE = f * (e_s-e_a)

    References
    ----------
    .. [haude_1955] Haude, W. (1955). Determination of evapotranspiration by
        an approach as simple as possible. Mitt Dt Wetterdienst, 2(11).
    .. [schiff_1975] Schiff, H. (1975). Berechnung der potentiellen Verdunstung
        und deren Vergleich mit aktuellen Verdunstungswerten von Lysimetern.
        Archiv für Meteorologie, Geophysik und Bioklimatologie, Serie B, 23(4),
        331-342.
    """
    e0 = calc_e0(tmean)
    ea = rh * e0 / 100
    # Haude coefficients from [schiff_1975]_
    fk = [0.27, 0.27, 0.28, 0.39, 0.39, 0.37, 0.35, 0.33, 0.31, 0.29, 0.27,
          0.27]
    index, shape = get_index_shape(tmean)
    f = [fk[x - 1] for x in index.month]
    return k * f * (e0 - ea) * 10  # kPa to hPa


def hamon(tmean, lat, k=1, c=13.97, cc=218.527, method=1):
    """Evaporation calculated according to [hamon_1961]_.

    Parameters
    ----------
    tmean: pandas.Series, optional
        average day temperature [°C]
    lat: float, optional
        the site latitude [rad]
    k: float, optional
        calibration coefficient if method = 0 [-]
    c: float, optional
        c is a constant for calculation in mm per day if method = 1.
    cc: float, optional
        calibration coefficient if method = 2 [-].
    method: float, optional
        0 => Hamon after [oudin_2005]_
        1 => Hamon after equation 7 in [ansorge_2019]_
        2 => Hamon after equation 12 in [ansorge_2019]_.

    Returns
    -------
    pandas.Series containing the calculated evaporation.

    Examples
    --------
    >>> et_hamon = hamon_1(tmean, lat)

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
    .. [ansorge_2019] Ansorge, L., & Beran, A. (2019). Performance of simple
       temperature-based evaporation methods compared with a time series of pan
       evaporation measures from a standard 20 m 2 tank. Journal of Water and
       Land Development.
    """
    index, shape = get_index_shape(tmean)
    dl = broadcast_to(daylight_hours(index, lat, shape), shape)
    if method == 0:
        et = k * (dl / 12) ** 2 * exp(tmean / 16)
        return et[:]
    if method == 1:
        pt = 4.95 * exp(
            0.062 * tmean) / 100  # saturated water content after Xu and Singh (2001)
        et = c * (dl / 12) ** 2 * pt
        et = et.where(tmean > 0, 0)
        return et[:]
    if method == 2:
        et = cc * (dl / 12) * 1 / (tmean + 273.3) * exp(
            (17.26939 * tmean) / (tmean + 273.3))
        et = et.where(tmean > 0, 0)
        return et


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
