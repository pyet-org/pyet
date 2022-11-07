"""The radiation module contains functions of radiation PET methods

"""

from numpy import sqrt

from .combination import calc_lambda

from .meteo_utils import extraterrestrial_r, calc_press, calc_psy, calc_vpc

from pyet.utils import get_index, clip_zeros, check_lat


def turc(tmean, rs, rh, k=0.31, clip_zero=True):
    """Evaporation calculated according to :cite:t:`turc_estimation_1961`.

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
    Based on equation 2 and 3 in :cite:t:`xu_evaluation_2000`.

    .. math:: PE=k(\\frac{T_a}{T_a+15})(R_s/4.184 + 50)*4.184; for RH>50

    """
    c = tmean / tmean
    c.where(rh > 50, 1 + (50 - rh) / 70)
    pe = k * c * tmean / (tmean + 15) * (rs + 2.094)
    pe = clip_zeros(pe, clip_zero)
    return pe.rename("Turc")


def jensen_haise(tmean, rs=None, cr=0.025, tx=-3, lat=None, method=0,
                 clip_zero=True):
    """Potential evaporation calculated accordinf to
    :cite:t:`jensen_estimating_1963`.

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
        0 => after :cite:t:`jensen_evaporation_2016`
        1 => after :cite:t:`oudin_which_2005`
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
    Based on equation (K-15) in :cite:t:`jensen_evaporation_2016`.

    .. math:: PE = \\frac{C_r(T-T_x)R_s}{\\lambda}

    """
    lambd = calc_lambda(tmean)
    if method == 0:
        pe = rs / lambd * cr * (tmean - tx)
    elif method == 1:
        index = get_index(tmean)
        ra = extraterrestrial_r(index, check_lat(lat))
        pe = ra * (tmean + 5) / 68 / lambd
    pe = clip_zeros(pe, clip_zero)
    return pe.rename("Jensen_Haise")


def mcguinness_bordne(tmean, lat, k=0.0147, clip_zero=True):
    """Potential evaporation calculated according to
    :cite:t:`mcguinness_comparison_1972`.

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
    Based on equation 13 in :cite:t:`xu_evaluation_2000`.

    .. math:: PE = \\frac{0.0147 R_a (T_a + 5)}{\\lambda}

    """
    lambd = calc_lambda(tmean)
    index = get_index(tmean)
    ra = extraterrestrial_r(index, check_lat(lat))
    pe = k * ra * (tmean + 5) / lambd
    pe = clip_zeros(pe, clip_zero)
    return pe.rename("Mcguinness_Bordne")


def hargreaves(tmean, tmax, tmin, lat, k=0.0135, method=0, clip_zero=True):
    """Potential evaporation calculated according to
    :cite:t:`hargreaves_estimating_1982`.

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
        0 => after :cite:t:`jensen_evaporation_2016`
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
    Based on equation (8-16) in :cite:t:`jensen_evaporation_2016`.

    .. math:: PE = 0.0023 \\frac{R_a (T_a+17.8)\\sqrt{(T_{max}-T_{min})}}\
        {\\lambda}

    """
    lambd = calc_lambda(tmean)
    index = get_index(tmean)

    ra = extraterrestrial_r(index, check_lat(lat))
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
    """Potential evaporation calculated according to
    :cite:t:`jensen_evapotranspiration_1990`.

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
    """Potential evaporation calculated according to
    :cite:t:`abtew_evapotranspiration_1996`.

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
    Based on equation 14 in :cite:t:`xu_evaluation_2000`.

    .. math:: PE = \\frac{k R_s}{\\lambda}

    """
    lambd = calc_lambda(tmean)
    pe = k * rs / lambd
    pe = clip_zeros(pe, clip_zero)
    return pe.rename("Abtew")


def makkink(tmean, rs, pressure=None, elevation=None, k=0.65, clip_zero=True):
    """Potential evaporation calculated according to
    :cite:t:`makkink_testing_1957`.

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

    """
    pressure = calc_press(elevation, pressure)
    gamma = calc_psy(pressure)
    dlt = calc_vpc(tmean)
    lambd = calc_lambda(tmean)
    pe = k * dlt / (dlt + gamma) * rs / lambd
    pe = clip_zeros(pe, clip_zero)
    return pe.rename("Makkink")


def oudin(tmean, lat, k1=100, k2=5, clip_zero=True):
    """Potential evaporation calculated according to :cite:t:`oudin_which_2005`.

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
    Based on equation 3 in :cite:t:`oudin_which_2005`.

    .. math:: PE = \\frac{R_a (T_a +5)}{\\lambda 100}; if T_a+5>0
        else: P = 0

    """
    lambd = calc_lambda(tmean)
    index = get_index(tmean)
    # Add transpose to be able to work with lat in float or xarray.DataArray
    ra = extraterrestrial_r(index, check_lat(lat))
    pe = ra * (tmean + k2) / lambd / k1
    pe = pe.where(((tmean + k2) > 0) | (pe.isnull()), 0)
    pe = clip_zeros(pe, clip_zero)
    return pe.rename("Oudin")
