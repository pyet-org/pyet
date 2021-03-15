from .penman import calc_lambda, calc_ea, calc_es

from .utils import extraterrestrial_r


def oudin(tmean, lat, k1=100, k2=5):
    """Evaporation calculated according to [oudin_2005]_.

    Parameters
    ----------
    tmean: pandas.Series, optional
        average day temperature [°C]
    lat: float, optional
        the site latitude [rad]
    k1: float, optional
        calibration coefficient [-]
    k2: float, optional
        calibration coefficient [-]

    Returns
    -------
    pandas.Series containing the calculated evaporation.

    Examples
    --------
    >>> et_oudin = oudin(tmean, lat)

    Notes
    -----
        Based on equation 3 in [oudin_2005]_.

    .. math::
    -----
        E = \\frac{R_a (T_a +5)}{\\lambda \\rho 100};   if T_a+5>0
        P = 0,                                          otherwise

    References
    -----
    .. [oudin_2005] Oudin, L., Hervieu, F., Michel, C., Perrin, C.,
       Andréassian, V., Anctil, F., & Loumagne, C. (2005). Which potential
       evapotranspiration input for a lumped rainfall–runoff model?:
       Part 2—Towards a simple and efficient potential evapotranspiration model
       for rainfall–runoff modelling. Journal of hydrology, 303(1-4), 290-306.
    """
    lambd = calc_lambda(tmean)
    ra = extraterrestrial_r(tmean.index, lat)
    et = ra * (tmean + k2) / lambd / k1
    et[(tmean + k2) < 0] = 0
    return et


def abtew(tmean, rs, k=0.53):
    """Evaporation calculated according to [abtew_1996]_.

    Parameters
    ----------
    tmean: pandas.Series, optional
        average day temperature [°C]
    rs: pandas.Series, optional
        incoming solar radiation [MJ m-2 d-1]
    k: float, optional
        calibration coefficient [-]

    Returns
    -------
    pandas.Series containing the calculated evaporation.

    Examples
    --------
    >>> et_abtew = abtew(tmean, rs)

    Notes
    -----
        Based on equation 14 in [xu_2000]_.

    .. math::
    -----
        E = \\frac{k Rs}{\\lambda}

    References
    -----
    .. [abtew_1996] Abtew, W. (1996). Evapotranspiration measurements and
       modeling for three wetland systems in South Florida 1. JAWRA Journal of
       the American Water Resources Association, 32(3), 465-473.
    .. [xu_2000] Xu, C. Y., & Singh, V. P. (2000). Evaluation and
       generalization of radiation‐based methods for calculating evaporation.
       Hydrological processes, 14(2), 339-349.
    """
    lambd = calc_lambda(tmean)
    et = k * rs / lambd
    return et


def turc(tmean, rs, rh, k=0.0133):
    """Evaporation calculated according to [turc_1961]_.

    Parameters
    ----------
    tmean: pandas.Series
        average day temperature [°C]
    rs: pandas.Series
        incoming solar radiation [MJ m-2 d-1]
    rh: pandas.Series
        mean daily relative humidity [%]
    k: float, optional
        calibration coefficient [-]

    Returns
    -------
    pandas.Series containing the calculated evaporation.

    Examples
    --------
    >>> et_turc = turc(tmean, rs, rh)

    Notes
    -----
        Based on equation 2 and 3 in [xu_2000]_.

    .. math::
    -----
        E = (k * T_a / (T_a + 15) * (R_s/4.184 + 50)) / 4.184

    References
    -----
    .. [turc_1961] Xu, C‐Y., and V. P. Singh. "Evaluation and generalization of
       radiation‐based methods for calculating evaporation." Hydrological
       processes 14.2 (2000): 339-349.
    .. [xu_2000] Xu, C. Y., & Singh, V. P. (2000). Evaluation and
       generalization of radiation‐based methods for calculating evaporation.
       Hydrological processes, 14(2), 339-349.
    """
    et = k * tmean / (tmean + 15) * (rs / 4.184 + 50)
    et[rh < 50] = k * tmean / (tmean + 15) * (rs / 4.184 + 50) * (
            1 + (50 - rh) / 70)
    return et * 4.184


def mcguinness_bordne(tmean, lat, k=0.0147):
    """Evaporation calculated according to [mcguinness_bordne_1972]_.

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
    >>> et_mcguinness_bordne = mcguinness_bordne(tmean, lat)

    Notes
    -----
        Based on equation 13 in [xu_2000]_.

    .. math::
    -----
        E = \\frac{0.0147 R_a (T_a + 5)}{\\lambda}

    References
    -----
    .. [mcguinness_bordne_1972] McGuinness, J. L., & Bordne, E. F. (1972).
       A comparison of lysimeter derived potential evapotranspiration with
       computed values, Tech. Bull., 1452. Agric. Res. Serv., US Dep. of
       Agric., Washington, DC.
    .. [xu_2000] Xu, C. Y., & Singh, V. P. (2000). Evaluation and
       generalization of radiation‐based methods for calculating evaporation.
       Hydrological processes, 14(2), 339-349.
    """
    lambd = calc_lambda(tmean)
    ra = extraterrestrial_r(tmean.index, lat)
    et = k * ra * (tmean + 5) / lambd
    return et


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

    .. math::
    -----
        E = \\frac{500 T_m / (100-A)+15 (T_a-T_d)}{(80-T_a)}

    References
    -----
    .. [linacre_1977] Linacre, E. T. (1977). A simple formula for estimating
       evaporation rates in various climates, using temperature data alone.
       Agricultural meteorology, 18(6), 409-424.
    .. [linacre_1992] Linacre, E. (1992). Climate data and resources: a
       reference and guide. Psychology Press.
    .. [xu_2001] Xu, C. Y., & Singh, V. P. (2001). Evaluation and
       generalization of temperature‐based methods for calculating evaporation.
       Hydrological processes, 15(2), 305-319.
    """
    if tdew is None:
        tdew = 0.52 * tmin + 0.6 * tmax - 0.009 * tmax ** 2 - 2
    tm = tmean + 0.006 * elevation
    et = (500 * tm / (100 - lat) + 15 * (tmean - tdew)) / (80 - tmean)
    return et


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

    .. math::
    -----
        E =  kp(0.46 * T_a + 8.13)

    References
    -----
    .. [blaney_1952] Blaney, H. F. (1952). Determining water requirements in
       irrigated areas from climatological and irrigation data.
    .. [xu_2001] Xu, C. Y., & Singh, V. P. (2001). Evaluation and
       generalization of temperature‐based methods for calculating evaporation.
       Hydrological processes, 15(2), 305-319.
    """
    et = k * p * (0.46 * tmean + 8.13)
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

    .. math::
    -----
        E =  4.5 (1 + T_a/25)^2 (1 - \\frac{ea}{es})

    References
    -----
    .. [romanenko_1961] Romanenko, V. A. (1961). Computation of the autumn soil
       moisture using a universal relationship for a large area. Proc. of
       Ukrainian Hydrometeorological Research Institute, 3, 12-25.
    .. [xu_2001] Xu, C. Y., & Singh, V. P. (2001). Evaluation and
       generalization of temperature‐based methods for calculating evaporation.
       Hydrological processes, 15(2), 305-319.
    """
    ea = calc_ea(tmean=tmean, rh=rh)
    es = calc_es(tmean=tmean)

    return k * (1 + tmean / 25) ** 2 * (1 - ea / es)
