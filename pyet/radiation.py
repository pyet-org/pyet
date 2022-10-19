"""The radiation module contains functions of radiation PET methods

"""

from numpy import sqrt

from .combination import calc_lambda

from .meteo_utils import extraterrestrial_r, calc_press, calc_psy, calc_vpc

from .utils import *


def turc(tmean, rs, rh, k=0.31, clip_zero=True):
    """Evaporation calculated according to [turc_1961]_.

    Parameters
    ----------
    tmean: float/pandas.Series/xarray.DataArray
        average day temperature [°C]
    rs: float/pandas.Series/xarray.DataArray
        incoming solar radiation [MJ m-2 d-1]
    rh: float/pandas.Series/xarray.DataArray
        mean daily relative humidity [%]
    k: float, optional
        calibration coefficient [-]
    clip_zero: bool, optional
        if True, replace all negative values with 0.

    Returns
    -------
    float/pandas.Series/xarray.DataArray containing the calculated
            potential evaporation [mm d-1].

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
    c.where(rh > 50, 1 + (50 - rh) / 70)
    pe = k * c * tmean / (tmean + 15) * (rs + 2.094)
    pe = clip_zeros(pe, clip_zero)
    return pe.rename("Turc")


def jensen_haise(tmean, rs=None, cr=0.025, tx=-3, lat=None, method=1,
                 clip_zero=True):
    """Potential evaporation calculated accordinf to [jensen_haise_1963]_.

    Parameters
    ----------
    tmean: andas.Series/xarray.DataArray
        average day temperature [°C]
    rs: float/pandas.Series/xarray.DataArray, optional
        incoming solar radiation [MJ m-2 d-1]
    cr: float, optional
        temperature coefficient [-]
    tx: float, optional
        intercept of the temperature axis [°C]
    lat: float/xarray.DataArray
        the site latitude [rad]
    method: float, optional
        1 => after [jensen_allen_2016]
        2 => after [oudin_2005]
    clip_zero: bool, optional
        if True, replace all negative values with 0.

    Returns
    -------
    float/pandas.Series/xarray.DataArray containing the calculated
            potential evaporation [mm d-1].

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
        pe = rs / lambd * cr * (tmean - tx)
    elif method == 2:
        index = get_index(tmean)
        ra = extraterrestrial_r(index, lat)
        pe = ra * (tmean + 5) / 68 / lambd
    pe = clip_zeros(pe, clip_zero)
    return pe.rename("Jensen_Haise")


def mcguinness_bordne(tmean, lat, k=0.0147, clip_zero=True):
    """Potential evaporation calculated according to [mcguinness_bordne_1972]_.

    Parameters
    ----------
    tmean: pandas.Series/xarray.DataArray
        average day temperature [°C]
    lat: float/xarray.DataArray, optional
        the site latitude [rad]
    k: float, optional
        calibration coefficient [-]
    clip_zero: bool, optional
        if True, replace all negative values with 0.

    Returns
    -------
    float/pandas.Series/xarray.DataArray containing the calculated
            potential evaporation [mm d-1].

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
    index = get_index(tmean)
    ra = extraterrestrial_r(index, lat)
    pe = k * ra * (tmean + 5) / lambd
    pe = clip_zeros(pe, clip_zero)
    return pe.rename("Mcguinness_Bordne")


def hargreaves(tmean, tmax, tmin, lat, k=0.0135, method=0, clip_zero=True):
    """Potential evaporation calculated according to [hargreaves_samani_1982]_.

    Parameters
    ----------
    tmean: pandas.Series/xarray.DataArray
        average day temperature [°C]
    tmax: float/pandas.Series/xarray.DataArray
        maximum day temperature [°C]
    tmin: float/pandas.Series/xarray.DataArray
        minimum day temperature [°C]
    lat: float/xarray.DataArray
        the site latitude [rad]
    k: float, optional
        calirbation coefficient [-]
    method: float, optional
        0 => after [jensen_allen_2016]
        1 => after [mcmahon_2013]
    clip_zero: bool, optional
        if True, replace all negative values with 0.

    Returns
    -------
    float/pandas.Series/xarray.DataArray containing the calculated
            potential evaporation [mm d-1].

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
    index = get_index(tmean)

    ra = extraterrestrial_r(index, lat)
    chs = 0.00185 * (tmax - tmin) ** 2 - 0.0433 * (tmax - tmin) + 0.4023
    if method == 0:
        pe = k / 0.0135 * 0.0023 * (tmean + 17.8) * sqrt(
            tmax - tmin) * ra / lambd
    elif method == 1:
        pe = k * chs * sqrt(tmax - tmin) * ra / lambd * (tmean + 17.8)
    pe = clip_zeros(pe, clip_zero)
    return pe.rename("Hargreaves")


def fao_24(tmean, wind, rs, rh, pressure=None, elevation=None, albedo=0.23,
           clip_zero=True):
    """Potential evaporation calculated according to [jensen_1990]_.

    Parameters
    ----------
    tmean: float/pandas.Series/xarray.DataArray
        average day temperature [°C]
    wind: float/pandas.Series/xarray.DataArray
        mean day wind speed [m/s]
    rs: float/pandas.Series/xarray.DataArray
        incoming solar radiation [MJ m-2 d-1]
    rh: float/pandas.Series/xarray.DataArray
        mean daily relative humidity [%]
    pressure: float/pandas.Series/xarray.DataArray, optional
        atmospheric pressure [kPa]
    elevation: float/xarray.DataArray, optional
        the site elevation [m]
    albedo: float/xarray.DataArray, optional
        surface albedo [-]
    clip_zero: bool, optional
        if True, replace all negative values with 0.

    Returns
    -------
    float/pandas.Series/xarray.DataArray containing the calculated
            potential evaporation [mm d-1].

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
    pressure = calc_press(elevation, pressure)
    gamma = calc_psy(pressure)
    dlt = calc_vpc(tmean)
    lambd = calc_lambda(tmean)

    w = 1.066 - 0.13 * rh / 100 + 0.045 * wind - 0.02 * rh / 100 * wind - \
        0.315 * (rh / 100) ** 2 - 0.0011 * wind
    pe = -0.3 + dlt / (dlt + gamma) * rs * (1 - albedo) * w / lambd
    pe = clip_zeros(pe, clip_zero)
    return pe.rename("FAO_24")


def abtew(tmean, rs, k=0.53, clip_zero=True):
    """Potential evaporation calculated according to [abtew_1996]_.

    Parameters
    ----------
    tmean: float/pandas.Series/xarray.DataArray
        average day temperature [°C]
    rs: float/pandas.Series/xarray.DataArray
        incoming solar radiation [MJ m-2 d-1]
    k: float, optional
        calibration coefficient [-]
    clip_zero: bool, optional
        if True, replace all negative values with 0.

    Returns
    -------
    float/pandas.Series/xarray.DataArray containing the calculated
            potential evaporation [mm d-1].

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
    pe = k * rs / lambd
    pe = clip_zeros(pe, clip_zero)
    return pe.rename("Abtew")


def makkink(tmean, rs, pressure=None, elevation=None, k=0.65, clip_zero=True):
    """"Potential evaporation calculated according to [makkink1957]_.

    Parameters
    ----------
    tmean: float/pandas.Series/xarray.DataArray
        average day temperature [°C]
    rs: float/pandas.Series/xarray.DataArray
        incoming solar radiation [MJ m-2 d-1]
    pressure: float/pandas.Series/xarray.DataArray, optional
        atmospheric pressure [kPa]
    elevation: float/xarray.DataArray, optional
        the site elevation [m]
    k: float, optional
        calirbation coefficient [-]
    clip_zero: bool, optional
        if True, replace all negative values with 0.

    Returns
    -------
    float/pandas.Series/xarray.DataArray containing the calculated
            potential evaporation [mm d-1].

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
    pressure = calc_press(elevation, pressure)
    gamma = calc_psy(pressure)
    dlt = calc_vpc(tmean)
    lambd = calc_lambda(tmean)
    pe = k * dlt / (dlt + gamma) * rs / lambd
    pe = clip_zeros(pe, clip_zero)
    return pe.rename("Makkink")


def oudin(tmean, lat, k1=100, k2=5, clip_zero=True):
    """Potential evaporation calculated according to [oudin_2005]_.

    Parameters
    ----------
    tmean: pandas.Series/xarray.DataArray
        average day temperature [°C]
    lat: float/xarray.DataArray
        the site latitude [rad]
    k1: float, optional
        calibration coefficient [-]
    k2: float, optional
        calibration coefficient [-]
    clip_zero: bool, optional
        if True, replace all negative values with 0.

    Returns
    -------
    float/pandas.Series/xarray.DataArray containing the calculated
            potential evaporation [mm d-1].
    clip_zero: bool, optional
        if True, replace all negative values with 0.

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
    index = get_index(tmean)
    # Add transpose to be able to work with lat in float or xarray.DataArray
    ra = extraterrestrial_r(index, lat)
    pe = ra * (tmean + k2) / lambd / k1
    pe = pe.where(((tmean + k2) > 0) | (pe.isnull()), 0)
    pe = clip_zeros(pe, clip_zero)
    return pe.rename("Oudin")
