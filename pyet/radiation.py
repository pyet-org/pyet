"""The radiation module contains functions of radiation PET methods

"""

from numpy import sqrt

from .combination import calc_lambda

from .meteo_utils import extraterrestrial_r, calc_press, calc_psy, calc_vpc

from pyet.utils import get_index, clip_zeros, check_lat


def turc(tmean, rs, rh, k=0.31, clip_zero=True):
    """Potential evapotranspiration calculated according to
    :cite:t:`turc_estimation_1961`.

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
            Potential evapotranspiration [mm d-1].

    Examples
    --------
    >>> iet_turc = turc(tmean, rs, rh)

    Notes
    -----
    Based on equation 2 and 3 in :cite:t:`xu_evaluation_2000`.

    .. math:: PET=k(\\frac{T_{mean}}{T_{mean}+15})(\\frac{R_s}{4.184}+50)4.184;
        for RH>50

    .. math:: PET=k(1+\\frac{50-RH}{70})(\\frac{T_{mean}}{T_{mean}+15})
        (\\frac{R_s}{4.184}+50)4.184; for RH<50

    """
    c = tmean / tmean
    c.where(rh > 50, 1 + (50 - rh) / 70)
    pet = k * c * tmean / (tmean + 15) * (rs + 2.094)
    pet = clip_zeros(pet, clip_zero)
    return pet.rename("Turc")


def jensen_haise(tmean, rs=None, cr=0.025, tx=-3, lat=None, method=0,
                 clip_zero=True):
    """Potential evapotranspiration calculated accordinf to
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
    lambd = calc_lambda(tmean)
    if method == 0:
        pet = rs / lambd * cr * (tmean - tx)
    elif method == 1:
        index = get_index(tmean)
        ra = extraterrestrial_r(index, check_lat(lat))
        pet = ra * (tmean + 5) / 68 / lambd
    pet = clip_zeros(pet, clip_zero)
    return pet.rename("Jensen_Haise")


def mcguinness_bordne(tmean, lat, k=0.0147, clip_zero=True):
    """Potential evapotranspiration calculated according to
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
            Potential evapotranspiration [mm d-1].

    Examples
    --------
    >>> pet_mcguinness_bordne = mcguinness_bordne(tmean, lat)

    Notes
    -----
    Based on equation 13 in :cite:t:`xu_evaluation_2000`.

    .. math:: PET = k\\frac{R_a (T_{mean} + 5)}{\\lambda}

    """
    lambd = calc_lambda(tmean)
    index = get_index(tmean)
    ra = extraterrestrial_r(index, check_lat(lat))
    pet = k * ra * (tmean + 5) / lambd
    pet = clip_zeros(pet, clip_zero)
    return pet.rename("Mcguinness_Bordne")


def hargreaves(tmean, tmax, tmin, lat, k=0.0135, method=0, clip_zero=True):
    """Potential evapotranspiration calculated according to
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
        1 => after :cite:t:`mcmahon_estimating_2013`
    clip_zero: bool, optional
        if True, replace all negative values with 0.

    Returns
    -------
    float/pandas.Series/xarray.DataArray containing the calculated
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
    lambd = calc_lambda(tmean)
    index = get_index(tmean)

    ra = extraterrestrial_r(index, check_lat(lat))
    chs = 0.00185 * (tmax - tmin) ** 2 - 0.0433 * (tmax - tmin) + 0.4023
    if method == 0:
        pet = k / 0.0135 * 0.0023 * (tmean + 17.8) * sqrt(
            tmax - tmin) * ra / lambd
    elif method == 1:
        pet = k * chs * sqrt(tmax - tmin) * ra / lambd * (tmean + 17.8)
    pet = clip_zeros(pet, clip_zero)
    return pet.rename("Hargreaves")


def fao_24(tmean, wind, rs, rh, pressure=None, elevation=None, albedo=0.23,
           clip_zero=True):
    """Potential evapotranspiration calculated according to
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
            Potential evapotranspiration [mm d-1].

    Examples
    --------
    >>> pet_fao24 = fao_24(tmean, wind, rs, rh, elevation=elevation)

    .. math:: PE = \\frac{- 0.3 \\Delta + R_s (1-\\alpha) w}\
        {\\lambda(\\Delta +\\gamma)}
    .. math:: w = 1.066-0.13*\\frac{rh}{100}+0.045*u_2-0.02*\\frac{rh}{100}\
        *u_2-3.15*(\\frac{rh}{100})^2-0.0011*u_2

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
    """Potential evapotranspiration calculated according to
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
            Potential evapotranspiration [mm d-1].

    Examples
    --------
    >>> pet_abtew = abtew(tmean, rs)

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
    """Potential evapotranspiration calculated according to
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
            Potential evapotranspiration [mm d-1].

    Examples
    --------
    >>> pet_mak = makkink(tmean, rs, elevation=elevation)

    Notes
    -----

    .. math:: PET = \\frac{0.65 \\Delta (R_s)}{\\lambda(\\Delta +\\gamma)}

    """
    pressure = calc_press(elevation, pressure)
    gamma = calc_psy(pressure)
    dlt = calc_vpc(tmean)
    lambd = calc_lambda(tmean)
    pet = k * dlt / (dlt + gamma) * rs / lambd
    pet = clip_zeros(pet, clip_zero)
    return pet.rename("Makkink")


def oudin(tmean, lat, k1=100, k2=5, clip_zero=True):
    """Potential evapotranspiration calculated according to
     :cite:t:`oudin_which_2005`.

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
    lambd = calc_lambda(tmean)
    index = get_index(tmean)
    ra = extraterrestrial_r(index, check_lat(lat))
    pet = ra * (tmean + k2) / lambd / k1
    pet = pet.where(((tmean + k2) > 0) | (pet.isnull()), 0)
    pet = clip_zeros(pet, clip_zero)
    return pet.rename("Oudin")
