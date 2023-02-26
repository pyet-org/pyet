"""The radiation module contains functions of radiation PET methods

"""

from numpy import sqrt, log, ones

from .combination import calc_lambda

from .meteo_utils import extraterrestrial_r, calc_press, calc_psy, calc_vpc

from pyet.utils import *


def turc(tmean, rs, rh, k=0.31, clip_zero=True):
    """Potential evapotranspiration calculated according to
    :cite:t:`turc_estimation_1961`.

    Parameters
    ----------
    tmean: pandas.Series or xarray.DataArray
        average day temperature [°C]
    rs: pandas.Series or xarray.DataArray
        incoming solar radiation [MJ m-2 d-1]
    rh: pandas.Series or xarray.DataArray
        mean daily relative humidity [%]
    k: float, optional
        calibration coefficient [-]
    clip_zero: bool, optional
        if True, replace all negative values with 0.

    Returns
    -------
    pandas.Series or xarray.DataArray containing the calculated
            Potential evapotranspiration [mm d-1].

    Examples
    --------
    >>> pet_turc = turc(tmean, rs, rh)

    Notes
    -----
    Based on equation 2 and 3 in :cite:t:`xu_evaluation_2000`.

    .. math:: PET=k(\\frac{T_{mean}}{T_{mean}+15})(\\frac{R_s}{4.184}+50)4.184;
        for RH>50

    .. math:: PET=k(1+\\frac{50-RH}{70})(\\frac{T_{mean}}{T_{mean}+15})
        (\\frac{R_s}{4.184}+50)4.184; for RH<50

    """
    vtmean, vrs, vrh = vectorize(tmean, rs, rh)
    c = ones(vtmean.shape)
    mask = check_rh(vrh) < 50
    c[mask] = (1 + (50 - vrh[mask]) / 70)
    vtmean[vtmean == 15] = 14.99
    epsilon = 1e-8  # small constant to avoid division by zero
    pet = k * c * vtmean / (vtmean + 15 + epsilon) * (check_rad(vrs) + 2.094)
    pet[vtmean < 0] = 0
    pet = clip_zeros(pet, clip_zero)
    return pet_out(tmean, pet, "Turc")


def jensen_haise(tmean, rs=None, cr=0.025, tx=-3, lat=None, method=0,
                 clip_zero=True):
    """Potential evapotranspiration calculated accordinf to
    :cite:t:`jensen_estimating_1963`.

    Parameters
    ----------
    tmean: pandas.Series orxarray.DataArray
        average day temperature [°C]
    rs: pandas.Series or xarray.DataArray, optional
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
    pandas.Series or xarray.DataArray containing the calculated
            Potential evapotranspiration [mm d-1].

    Examples
    --------
    >>> pet_jh = jensen_haise(tmean, lat)

    Notes
    -----
    Method = 0: Based on :cite:t:`jensen_evaporation_2016`.

    .. math:: PET = \\frac{C_r(T_{mean}-T_x)R_s}{\\lambda}

    Method = 1: Based on :cite:t:`oudin_which_2005`.

    .. math:: PET = \\frac{R_a(T_{mean}+5)}{68\\lambda}

    """
    vtmean, vrs = vectorize(tmean, rs)
    lambd = calc_lambda(vtmean)
    if method == 0:
        if vrs is None:
            raise Exception("If you choose method == 0, provide rs!")
        pet = check_rad(vrs) / lambd * cr * (vtmean - tx)
    elif method == 1:
        if lat is None:
            raise Exception("If you choose method == 1, provide lat!")
        index = get_index(tmean)
        ra = extraterrestrial_r(index, lat, vtmean.shape)
        pet = ra * (vtmean + 5) / 68 / lambd
    else:
        raise Exception("Method can be either 0 or 1.")
    pet = clip_zeros(pet, clip_zero)
    return pet_out(tmean, pet, "Jensen_Haise")


def mcguinness_bordne(tmean, lat, k=0.0147, clip_zero=True):
    """Potential evapotranspiration calculated according to
    :cite:t:`mcguinness_comparison_1972`.

    Parameters
    ----------
    tmean: pandas.Series or xarray.DataArray
        average day temperature [°C]
    lat: float/xarray.DataArray, optional
        the site latitude [rad]
    k: float, optional
        calibration coefficient [-]
    clip_zero: bool, optional
        if True, replace all negative values with 0.

    Returns
    -------
    pandas.Series or xarray.DataArray containing the calculated
            Potential evapotranspiration [mm d-1].

    Examples
    --------
    >>> pet_mcguinness_bordne = mcguinness_bordne(tmean, lat)

    Notes
    -----
    Based on equation 13 in :cite:t:`xu_evaluation_2000`.

    .. math:: PET = k\\frac{R_a (T_{mean} + 5)}{\\lambda}

    """
    vtmean, = vectorize(tmean)
    lambd = calc_lambda(vtmean)
    index = get_index(tmean)
    ra = extraterrestrial_r(index, lat, vtmean.shape)
    pet = k * ra * (vtmean + 5) / lambd
    pet = clip_zeros(pet, clip_zero)
    return pet_out(tmean, pet, "Mcguinness_Bordne")


def hargreaves(tmean, tmax, tmin, lat, k=0.0135, method=0, clip_zero=True):
    """Potential evapotranspiration calculated according to
    :cite:t:`hargreaves_estimating_1982`.

    Parameters
    ----------
    tmean: pandas.Series or xarray.DataArray
        average day temperature [°C]
    tmax: pandas.Series or xarray.DataArray
        maximum day temperature [°C]
    tmin: pandas.Series or xarray.DataArray
        minimum day temperature [°C]
    lat: float/xarray.DataArray
        the site latitude [rad]
    k: float, optional
        calirbation coefficient [-]
    method: float, optional
        0 => after :cite:t:`jensen_evaporation_2016`
        1 => after :cite:t:`mcmahon_estimating_2013`
    clip_zero: bool, optional
        if True, replace all negative values with 0.

    Returns
    -------
    pandas.Series or xarray.DataArray containing the calculated
            Potential evapotranspiration [mm d-1].

    Examples
    --------
    >>> pet_har = hargreaves(tmean, tmax, tmin, lat)

    Notes
    -----
    Method = 0; Based on equation (8-16) in :cite:t:`jensen_evaporation_2016`.

    .. math:: PET = k \\frac{R_a (T_{mean}+17.8)\\sqrt{(T_{max}-T_{min})}}\
        {\\lambda}

    Method = 1; Based on :cite:t:`mcmahon_estimating_2013`.

    .. math:: PET = chs k \\frac{R_a (T_{mean}+17.8)\\sqrt{(T_{max}-T_{min})}}\
        {\\lambda}

    , where

    .. math:: chs=0.00185*(T_{max}-T_{min})^2-0.0433*(T_{max}-T_{min})+0.4023

    """
    vtmean, vtmax, vtmin = vectorize(tmean, tmax, tmin)
    lambd = calc_lambda(vtmean)
    index = get_index(tmean)

    ra = extraterrestrial_r(index, lat, vtmean.shape)
    chs = 0.00185 * (vtmax - vtmin) ** 2 - 0.0433 * (vtmax - vtmin) + 0.4023
    if method == 0:
        pet = k / 0.0135 * 0.0023 * (vtmean + 17.8) * sqrt(
            vtmax - vtmin) * ra / lambd
    elif method == 1:
        pet = k * chs * sqrt(vtmax - vtmin) * ra / lambd * (vtmean + 17.8)
    else:
        raise Exception("Method can be either 0 or 1.")
    pet = clip_zeros(pet, clip_zero)
    return pet_out(tmean, pet, "Hargreaves")


def fao_24(tmean, wind, rs, rh, pressure=None, elevation=None, albedo=0.23,
           clip_zero=True):
    """Potential evapotranspiration calculated according to
    :cite:t:`jensen_evapotranspiration_1990`.

    Parameters
    ----------
    tmean: pandas.Series or xarray.DataArray
        average day temperature [°C]
    wind: pandas.Series xarray.DataArray
        mean day wind speed [m/s]
    rs: pandas.Series xarray.DataArray
        incoming solar radiation [MJ m-2 d-1]
    rh: pandas.Series or xarray.DataArray
        mean daily relative humidity [%]
    pressure: pandas.Series or xarray.DataArray, optional
        atmospheric pressure [kPa]
    elevation: float/xarray.DataArray, optional
        the site elevation [m]
    albedo: float/xarray.DataArray, optional
        surface albedo [-]
    clip_zero: bool, optional
        if True, replace all negative values with 0.

    Returns
    -------
    pandas.Series or xarray.DataArray containing the calculated
            Potential evapotranspiration [mm d-1].

    Examples
    --------
    >>> pet_fao24 = fao_24(tmean, wind, rs, rh, elevation=elevation)

    .. math:: PE = \\frac{- 0.3 \\Delta + R_s (1-\\alpha) w}\
        {\\lambda(\\Delta +\\gamma)}
    .. math:: w = 1.066-0.13*\\frac{rh}{100}+0.045*u_2-0.02*\\frac{rh}{100}\
        *u_2-3.15*(\\frac{rh}{100})^2-0.0011*u_2

    """
    vtmean, vwind, vrs, vrh, vpressure, velevation = vectorize(tmean, wind, rs,
                                                               rh, pressure,
                                                               elevation)
    vpressure = calc_press(velevation, vpressure)
    gamma = calc_psy(vpressure)
    dlt = calc_vpc(vtmean)
    lambd = calc_lambda(vtmean)

    w = 1.066 - 0.13 * check_rh(
        vrh) / 100 + 0.045 * vwind - 0.02 * vrh / 100 * vwind - \
        0.315 * (vrh / 100) ** 2 - 0.0011 * vwind
    pet = -0.3 + dlt / (dlt + gamma) * check_rad(vrs) * (
            1 - albedo) * w / lambd
    pet = clip_zeros(pet, clip_zero)
    return pet_out(tmean, pet, "FAO_24")


def abtew(tmean, rs, k=0.53, clip_zero=True):
    """Potential evapotranspiration calculated according to
    :cite:t:`abtew_evapotranspiration_1996`.

    Parameters
    ----------
    tmean: pandas.Series or xarray.DataArray
        average day temperature [°C]
    rs: pandas.Series or xarray.DataArray
        incoming solar radiation [MJ m-2 d-1]
    k: float, optional
        calibration coefficient [-]
    clip_zero: bool, optional
        if True, replace all negative values with 0.

    Returns
    -------
    pandas.Series or xarray.DataArray containing the calculated
            Potential evapotranspiration [mm d-1].

    Examples
    --------
    >>> pet_abtew = abtew(tmean, rs)

    Notes
    -----
    Based on equation 14 in :cite:t:`xu_evaluation_2000`.

    .. math:: PE = \\frac{k R_s}{\\lambda}

    """
    vtmean, vrs = vectorize(tmean, rs)
    lambd = calc_lambda(vtmean)
    pet = k * check_rad(vrs) / lambd
    pet = clip_zeros(pet, clip_zero)
    return pet_out(tmean, pet, "Abtew")


def makkink(tmean, rs, pressure=None, elevation=None, k=0.65, clip_zero=True):
    """Potential evapotranspiration calculated according to
    :cite:t:`makkink_testing_1957`.

    Parameters
    ----------
    tmean: pandas.Series or xarray.DataArray
        average day temperature [°C]
    rs: pandas.Series or xarray.DataArray
        incoming solar radiation [MJ m-2 d-1]
    pressure: pandas.Series or xarray.DataArray, optional
        atmospheric pressure [kPa]
    elevation: float/xarray.DataArray, optional
        the site elevation [m]
    k: float, optional
        calirbation coefficient [-]
    clip_zero: bool, optional
        if True, replace all negative values with 0.

    Returns
    -------
    pandas.Series or xarray.DataArray containing the calculated
            Potential evapotranspiration [mm d-1].

    Examples
    --------
    >>> pet_mak = makkink(tmean, rs, elevation=elevation)

    Notes
    -----

    .. math:: PET = \\frac{0.65 \\Delta (R_s)}{\\lambda(\\Delta +\\gamma)}

    """
    vtmean, vrs, vpressure, velevation = vectorize(tmean, rs, pressure,
                                                   elevation)
    vpressure = calc_press(velevation, vpressure)
    gamma = calc_psy(vpressure)
    dlt = calc_vpc(vtmean)
    lambd = calc_lambda(vtmean)
    pet = k * dlt / (dlt + gamma) * check_rad(vrs) / lambd
    pet = clip_zeros(pet, clip_zero)
    return pet_out(tmean, pet, "Makkink")


def makkink_knmi(tmean, rs, clip_zero=True):
    """Potential evapotranspiration calculated according to The Royal
    Netherlands Meteorological Institute (KNMI)

    Parameters
    ----------
    tmean : pandas.Series or xarray.DataArray
        average day temperature [°C]
    rs : pandas.Series or xarray.DataArray
        incoming solar radiation [MJ m-2 d-1]
    clip_zero: bool, optional
        if True, replace all negative values with 0.

    Returns
    -------
    pandas.Series or xarray.DataArray containing the calculated
            Potential evapotranspiration [mm d-1].

    Examples
    --------
    >>> pet_mak = makkink_knmi(tmean, rs)

    Notes
    ----
    This method is only applicable to the Netherlands (~sea level) due to some
    emperical values. Calculating the Makkink evaporation with the original
    formula is more suitable for general purposes. To obtain the same value as
    EV24 round the value up to one decimal.
    """
    vtmean, vrs = vectorize(tmean, rs)
    pet = (
            650
            * (
                    1
                    - (0.646 + 0.0006 * vtmean)
                    / (
                            7.5
                            * log(10)
                            * 6.107
                            * 10 ** (7.5 * (1 - 1 / (1 + vtmean / 237.3)))
                            / (237.3 * (1 + vtmean / 237.3) * (
                            1 + vtmean / 237.3))
                            + 0.646
                            + 0.0006 * vtmean
                    )
            )
            / (2501 - 2.38 * vtmean)
            * vrs
    )
    pet = clip_zeros(pet, clip_zero)
    return pet_out(tmean, pet, "Makkink_KNMI")


def oudin(tmean, lat, k1=100, k2=5, clip_zero=True):
    """Potential evapotranspiration calculated according to
     :cite:t:`oudin_which_2005`.

    Parameters
    ----------
    tmean: pandas.Series or xarray.DataArray
        average day temperature [°C]
    lat: float or xarray.DataArray
        the site latitude [rad]
    k1: float, optional
        calibration coefficient [-]
    k2: float, optional
        calibration coefficient [-]
    clip_zero: bool, optional
        if True, replace all negative values with 0.

    Returns
    -------
    pandas.Series or xarray.DataArray containing the calculated
            Potential evapotranspiration [mm d-1].
    clip_zero: bool, optional
        if True, replace all negative values with 0.

    Examples
    --------
    >>> pet_oudin = oudin(tmean, lat)

    Notes
    -----
    Based on equation 3 in :cite:t:`oudin_which_2005`.

    .. math:: PET = \\frac{R_a (T_{mean} +5)}{\\lambda 100}; if T_{mean}+5>0
        else: PET = 0

    """
    vtmean, = vectorize(tmean)
    lambd = calc_lambda(vtmean)
    index = get_index(tmean)
    ra = extraterrestrial_r(index, lat, vtmean.shape)
    pet = ra * (vtmean + k2) / lambd / k1
    pet[(tmean + k2) < 0] = 0
    pet = clip_zeros(pet, clip_zero)
    return pet_out(tmean, pet, "Oudin")
