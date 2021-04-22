"""The radiation module contains functions of radiation PET methods

"""

from numpy import sqrt

from .combination import calc_lambda

from .meteo_utils import extraterrestrial_r, calc_press, calc_psy, calc_vpc


def turc(tmean, rs, rh, k=0.31):
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

    .. math:: PE=k(\\frac{T_a}{T_a+15})(R_s/4.184 + 50)*4.184; for RH>50

    References
    ----------
    .. [turc_1961] Xu, C‐Y., and V. P. Singh. "Evaluation and generalization of
       radiation‐based methods for calculating evaporation." Hydrological
       processes 14.2 (2000): 339-349.
    .. [xu_2000] Xu, C. Y., & Singh, V. P. (2000). Evaluation and
       generalization of radiation‐based methods for calculating evaporation.
       Hydrological processes, 14(2), 339-349.
    """
    c = tmean / tmean
    c[rh] = 1 - (50 - rh) / 70
    et = k * c * tmean / (tmean + 15) * (rs + 2.094)
    return et


def jensen_haise(tmean, rs=None, cr=0.025, tx=-3, lat=None, method=1):
    """Evaporation calculated accordinf to [jensen_haise_1963]_.

    Parameters
    ----------
    tmean: pandas.Series, optional
        average day temperature [°C]
    rs: pandas.Series, optional
        incoming solar radiation [MJ m-2 d-1]
    cr: float, optional
        temperature coefficient [-]
    tx: float, optional
        intercept of the temperature axis [°C]
    lat: float
        the site latitude [rad]
    method: float, optional
        1 => after [jensen_allen_2016]
        2 => after [oudin_2005]

    Returns
    -------
    pandas.Series containing the calculated evaporation.

    Examples
    --------
    >>> jh_et = jensen_haise(tmean, rs)

    Notes
    -----
    Based on equation (K-15) in [jensen_allen_2016]_.

    .. math:: PE = \\frac{C_r(T-T_x)R_s}{\\lambda}

    References
    ----------
    .. [jensen_haise_1963] Jensen, M.E., Haise, H.R., 1963. Estimating
       evapotranspiration from solar radiation. Journal of Irrigation and
       Drainage Division, ASCE 89 (LR4), 15–41.
    .. [jensen_allen_2016] Task Committee on Revision of Manual 70. (2016).
       Evaporation, evapotranspiration, and irrigation water requirements.
       American Society of Civil Engineers.
    """
    lambd = calc_lambda(tmean)
    if method == 1:
        return rs / lambd * cr * (tmean - tx)
    else:
        ra = extraterrestrial_r(tmean.index, lat)
        return ra * (tmean + 5) / 68 / lambd


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

    .. math:: PE = \\frac{0.0147 R_a (T_a + 5)}{\\lambda}

    References
    ----------
    .. [mcguinness_bordne_1972] McGuinness, J. L., & Bordne, E. F. (1972).
       A comparison of lysimeter derived potential evapotranspiration with
       computed values, Tech. Bull., 1452. Agric. Res. Serv., US Dep. of
       Agric., Washington, DC.
    """
    lambd = calc_lambda(tmean)
    ra = extraterrestrial_r(tmean.index, lat)
    et = k * ra * (tmean + 5) / lambd
    return et


def hargreaves(tmean, tmax, tmin, lat):
    """Evaporation calculated according to [hargreaves_samani_1982]_.

        Parameters
        ----------
        tmean: pandas.Series
            average day temperature [°C]
        tmax: pandas.Series, optional
            maximum day temperature [°C]
        tmin: pandas.Series, optional
            minimum day temperature [°C]
        lat: float, optional
            the site latitude [rad]

        Returns
        -------
        pandas.Series containing the calculated evaporation.

        Examples
        --------
        >>> et_har = hargreaves(tmean, tmax, tmin, lat)

        Notes
        -----
        Based on equation (8-16) in [jensen_allen_2016]_.

        .. math:: PE = 0.0023 \\frac{R_a (T_a+17.8)\\sqrt{(T_{max}-T_{min})}}\
            {\\lambda}

        References
        ----------
        .. [hargreaves_samani_1982] Hargreaves, G. H., & Samani, Z. A. (1982).
           Estimating potential evapotranspiration. Journal of the irrigation
           and Drainage Division, 108(3), 225-230.
        """
    lambd = calc_lambda(tmean)
    ra = extraterrestrial_r(tindex=tmean.index, lat=lat)
    return 0.0023 * (tmean + 17.8) * sqrt(tmax - tmin) * ra / lambd


def fao_24(tmean, wind, rs, rh, pressure=None, elevation=None, albedo=0.23):
    """Evaporation calculated according to [jensen_1990]_.

    Parameters
    ----------
    tmean: pandas.Series
        average day temperature [°C]
    wind: pandas.Series
        mean day wind speed [m/s]
    rs: pandas.Series
        incoming solar radiation [MJ m-2 d-1]
    rh: pandas.Series
        mean daily relative humidity [%]
    pressure: float, optional
        atmospheric pressure [kPa]
    elevation: float, optional
        the site elevation [m]
    albedo: float, optional
        surface albedo [-]

    Returns
    -------
        pandas.Series containing the calculated evaporation

    Examples
    --------
    >>> et_fao24 = fao_24(tmean, wind, rs, rh, pressure=pressure)

    .. math:: PE = \\frac{- 0.3 \\Delta + R_s (1-\\alpha) w}\
        {\\lambda(\\Delta +\\gamma)}
    .. math:: w = 1.066-0.13*\\frac{rh}{100}+0.045*u_2-0.02*\\frac{rh}{100}\
        *u_2-3.15*(\\frac{rh}{100})^2-0.0011*u_2$

    References
    ----------
    .. [jensen_1990] Jensen, M. E., Burman, R. D., & Allen, R. G. (1990).
       Evapotranspiration and irrigation water requirements. ASCE.

    """
    if pressure is None:
        pressure = calc_press(elevation)
    gamma = calc_psy(pressure)
    dlt = calc_vpc(tmean)
    lambd = calc_lambda(tmean)

    w = 1.066 - 0.13 * rh / 100 + 0.045 * wind - 0.02 * rh / 100 * wind - \
        0.315 * (rh / 100) ** 2 - 0.0011 * wind

    return -0.3 + dlt / (dlt + gamma) * rs * (1 - albedo) * w / lambd


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

    .. math:: PE = \\frac{k R_s}{\\lambda}

    References
    -----
    .. [abtew_1996] Abtew, W. (1996). Evapotranspiration measurements and
       modeling for three wetland systems in South Florida 1. JAWRA Journal of
       the American Water Resources Association, 32(3), 465-473.
    """
    lambd = calc_lambda(tmean)
    et = k * rs / lambd
    return et


def makkink(tmean, rs, pressure=None, elevation=None):
    """"Evaporation calculated according to [makkink1957]_.

    Parameters
    ----------
    tmean: pandas.Series
        average day temperature [°C]
    rs: pandas.Series, optional
        incoming solar radiation [MJ m-2 d-1]
    pressure: float, optional
        atmospheric pressure [kPa]
    elevation: float, optional
        the site elevation [m]

    Returns
    -------
        pandas.Series containing the calculated evaporation

    Examples
    --------
    >>> mak = makkink(tmean, rs)

    Notes
    -----

    .. math:: PE = \\frac{0.65 \\Delta (R_s)}{\\lambda(\\Delta +\\gamma)}

    References
    ----------
    .. [makkink1957] Makkink, G.F., 1957. Testing the Penman formula by means
        of lysimeters. Journal of the Institution of Water Engineers 11,
        277–288.
    """
    if pressure is None:
        pressure = calc_press(elevation)
    gamma = calc_psy(pressure)
    dlt = calc_vpc(tmean)
    lambd = calc_lambda(tmean)

    return 0.65 * dlt / (dlt + gamma) * rs / lambd


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

    .. math:: PE = \\frac{R_a (T_a +5)}{\\lambda 100}; if T_a+5>0
        else: P = 0

    """
    lambd = calc_lambda(tmean)
    ra = extraterrestrial_r(tmean.index, lat)
    et = ra * (tmean + k2) / lambd / k1
    et[(tmean + k2) < 0] = 0
    return et
