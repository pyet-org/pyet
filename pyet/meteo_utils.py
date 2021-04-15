"""The meteo_utils module contains utility functions for meteorological data

"""

from numpy import tan, cos, pi, sin, arccos, clip, maximum, exp, log
from pandas import to_numeric

# Specific heat of air [MJ kg-1 °C-1]
CP = 1.013 * 10 ** -3


def calc_psy(pressure, tmean=None):
    """Psychrometric constant [kPa °C-1].

    Parameters
    ----------
    pressure: float
        atmospheric pressure [kPa].
    tmean: float, optional
        average day temperature [°C].

    Returns
    -------
        pandas.Series containing the Psychrometric constant [kPa °C-1].

    Examples
    --------
    >>> psy = calc_psy(pressure, tmean)

    Notes
    -----
    if tmean is none:
        Based on equation 8 in [allen_1998]_.
    elif rh is None:
        From FAO (1990), ANNEX V, eq. 4.

    References
    ----------
    .. [allen_1998] Allen, R. G., Pereira, L. S., Raes, D., & Smith, M. (1998).
       Crop evapotranspiration-Guidelines for computing crop water
       requirements-FAO Irrigation and drainage paper 56. Fao, Rome, 300.
       (http://www.fao.org/3/x0490e/x0490e06.htm#TopOfPage).
    """
    if tmean is None:
        return 0.000665 * pressure
    else:
        lambd = calc_lambda(tmean)  # MJ kg-1
        return CP * pressure / (0.622 * lambd)


def calc_vpc(tmean):
    """Slope of saturation vapour pressure curve at air Temperature [kPa °C-1].

    Parameters
    ----------
    tmean: pandas.Series, optional
        average day temperature [°C]

    Returns
    -------
        pandas.Series containing the calculated Saturation vapour pressure
        [kPa °C-1].

    Examples
    --------
    >>> vpc = calc_vpc(tmean)

    Notes
    -----
    Based on equation 13. in [allen_1998]_.
    """
    es = calc_e0(tmean)
    return 4098 * es / (tmean + 237.3) ** 2


def calc_lambda(tmean):
    """ Latent Heat of Vaporization [MJ kg-1].

    Parameters
    ----------
    tmean: pandas.Series/float, optional
        average day temperature [°C]

    Returns
    -------
    pandas.Series containing the calculated Latent Heat of Vaporization
        [MJ kg-1].

    Examples
    --------
    >>> lambd = calc_lambda(tmean)

    Notes
    -----
    Based on equation (3-1) in [allen_1998]_.
    """
    return 2.501 - 0.002361 * tmean


def calc_press(elevation):
    """Atmospheric pressure [kPa].

    Parameters
    ----------
    elevation: float, optional
        the site elevation [m]

    Returns
    -------
    pandas.Series containing the calculated atmospheric pressure [kPa].

    Examples
    --------
    >>> pressure = calc_press(elevation)

    Notes
    -----
    Based on equation 7 in [allen_1998]_.
    """
    return 101.3 * ((293 - 0.0065 * elevation) / 293) ** 5.26


def calc_rho(pressure, tmean, ea):
    """atmospheric air density calculated according to [allen_1998]_..

    Parameters
    ----------
    pressure: pandas.Series/float
        atmospheric pressure [kPa]
    tmean: pandas.Series/float, optional
        average day temperature [°C]
    ea: pandas.Series/float, optional
        actual vapour pressure [kPa]

    Returns
    -------
    pandas.Series containing the calculated mean air density

    Examples
    --------
    >>> rho = calc_rho(pressure, tmean, ea)

    Notes
    -----
    Based on equation (3-5) in [allen_1998]_.

    .. math:: rho = 3.486 \\frac{P}{T_{KV}}
    """
    tkv = (273.16 + tmean) * (1 - 0.378 * ea / pressure) ** -1
    return 3.486 * pressure / tkv


def calc_e0(tmean):
    """ Saturation vapor pressure at the air temperature T [kPa].

    Parameters
    ----------
    tmean: pandas.Series, optional
        average day temperature [°C]

    Returns
    -------
    pandas.Series containing the calculated saturation vapor pressure at the
        air temperature tmean [kPa].

    Examples
    --------
    >>> e0 = calc_e0(tmean)

    Notes
    -----
    Based on equation 11 in [allen_1998]_.
    """
    return 0.6108 * exp(17.27 * tmean / (tmean + 237.3))


def calc_es(tmean=None, tmax=None, tmin=None):
    """ Saturation vapor pressure [kPa].

    Parameters
    ----------
    tmean: pandas.Series, optional
        average day temperature [°C]
    tmax: pandas.Series, optional
        maximum day temperature [°C]
    tmin: pandas.Series, optional
        minimum day temperature [°C]

    Returns
    -------
    pandas.Series containing the calculated saturation vapor pressure [kPa].

    Examples
    --------
    >>> es = calc_es(tmean)

    Notes
    -----
    Based on equation 11, 12 in [allen_1998]_.
    """
    if tmax is not None:
        eamax = calc_e0(tmax)
        eamin = calc_e0(tmin)
        return (eamax + eamin) / 2
    else:
        return calc_e0(tmean)


def calc_ea(tmean=None, tmax=None, tmin=None, rhmax=None, rhmin=None, rh=None):
    """ Actual vapor pressure [kPa].

    Parameters
    ----------
    tmean: pandas.Series, optional
        average day temperature [°C]
    tmax: pandas.Series, optional
        maximum day temperature [°C]
    tmin: pandas.Series, optional
        minimum day temperature [°C]
    rhmax: pandas.Series, optional
        maximum daily relative humidity [%]
    rhmin: pandas.Series, optional
        mainimum daily relative humidity [%]
    rh: pandas.Series, optional
        mean daily relative humidity [%]

    Returns
    -------
    pandas.Series containing the calculated actual vapor pressure [kPa].

    Examples
    --------
    >>> ea = calc_ea(tmean, rh)

    Notes
    -----
    Based on equation 17, 19 in [allen_1998]_.
    """
    if rhmax is not None:  # eq. 11
        esmax = calc_e0(tmax)
        esmin = calc_e0(tmin)
        return (esmin * rhmax / 200) + (esmax * rhmin / 200)
    else:  # eq. 14
        if tmax is not None:
            es = calc_es(tmax=tmax, tmin=tmin)
        else:
            es = calc_e0(tmean)
        return rh / 100 * es


def day_of_year(tindex):
    """Day of the year (1-365) based on pandas.Index

    Parameters
    ----------
    tindex: pandas.Index

    Returns
    -------
    array of with ints specifying day of year.

    """
    return to_numeric(tindex.strftime('%j'))


def daylight_hours(tindex, lat):
    """Daylight hours [hour].

    Parameters
    ----------
    tindex: pandas.Index
    lat: float
        the site latitude [rad]

    Returns
    -------
    pandas.Series containing the calculated daylight hours [hour]

    Notes
    -----
    Based on equation 34 in [allen_1998]_.
    """
    j = day_of_year(tindex)
    sol_dec = solar_declination(j)
    sangle = sunset_angle(sol_dec, lat)
    return 24 / pi * sangle


def sunset_angle(sol_dec, lat):
    """Sunset hour angle from latitude and solar declination - daily [rad].

    Parameters
    ----------
    sol_dec: pandas.Series
        solar declination [rad]
    lat: float
        the site latitude [rad]

    Returns
    -------
    pandas.Series containing the calculated sunset hour angle - daily [rad]

    Notes
    -----
    Based on equations 25 in [allen_1998]_.
    """
    return arccos(-tan(sol_dec) * tan(lat))


def sunset_angle_hour(tindex, sol_dec, lat, lz, lm):
    """Sunset hour angle from latitude and solar declination - hourly [rad].

    Parameters
    ----------
    tindex: pandas.Index
    sol_dec: pandas.Series
        solar declination [rad]
    lat: float
        the site latitude [rad]
    lz: float
        longitude of the local time zone [°]
    lm: float
        longitude of the measurement site [°]

    Returns
    -------
    pandas.Series containing the calculated sunset hour angle - hourly [rad]

    Notes
    -----
    Based on equations 29, 30, 31, 32, 33 in [allen_1998]_.
    """
    j = day_of_year(tindex)
    b = 2 * pi * (j - 81) / 364
    sc = 0.1645 * sin(2 * b) - 0.1255 * cos(b) - 0.025 * sin(b)
    t = tindex.hour + 0.5
    sol_t = t + 0.06667 * (lz - lm) + sc - 12  # equation 31
    omega = pi / 12 * sol_t

    omega1 = omega - pi / 24
    omega2 = omega + pi / 24

    omegas = arccos(-tan(lat) * tan(sol_dec))

    omega1 = clip(omega1, -omegas, omegas)
    omega2 = clip(omega2, -omegas, omegas)
    omega1 = maximum(omega1, omega1, )
    omega1 = clip(omega1, -100000000, omega2)

    return omega2, omega1


def solar_declination(j):
    """Solar declination from day of year [rad].

    Parameters
    ----------
    j: array.py
        day of the year (1-365)
    Returns
    -------
    array.py of solar declination [rad].

    Notes
    -------
    Based on equations 24 in [allen_1998]_.
    """
    return 0.409 * sin(2. * pi / 365. * j - 1.39)


def relative_distance(j):
    """Inverse relative distance between earth and sun from day of the year.

    Parameters
    ----------
    j: array.py
        day of the year (1-365)
    Returns
    -------
    array.py specifyng relative distance between earth and sun.

    Notes
    -------
    Based on equations 23 in [allen_1998]_.
    """
    return 1 + 0.033 * cos(2. * pi / 365. * j)


def extraterrestrial_r(tindex, lat):
    """Extraterrestrial daily radiation [MJ m-2 d-1].

    Parameters
    ----------
    tindex: pandas.Index
    lat: float
        the site latitude [rad]

    Returns
    -------
    pandas.Series containing the calculated extraterrestrial radiation

    Notes
    -----
    Based on equation 21 in [allen_1998]_.
    """
    j = day_of_year(tindex)
    dr = relative_distance(j)
    sol_dec = solar_declination(j)

    omega = sunset_angle(sol_dec, lat)
    xx = sin(sol_dec) * sin(lat)
    yy = cos(sol_dec) * cos(lat)
    return 118.08 / 3.141592654 * dr * (omega * xx + yy * sin(omega))


def extraterrestrial_r_hour(tindex, lat, lz=0, lm=0):
    """Extraterrestrial hourly radiation [MJ m-2 h-1].

    Parameters
    ----------
    tindex: pandas.Index
    lat: float
        the site latitude [rad]
    lz: float, optional
        longitude of the centre of the local time zone (0° for Greenwich) [°]
    lm: float, optional
        longitude of the measurement site [degrees west of Greenwich] [°]

    Returns
    -------
    pandas.Series containing the calculated extraterrestrial radiation

    Notes
    -----
    Based on equation 28 in [allen_1998]_.

    """
    j = day_of_year(tindex)
    dr = relative_distance(j)
    sol_dec = solar_declination(j)

    omega2, omega1 = sunset_angle_hour(tindex, lz=lz, lm=lm, lat=lat,
                                       sol_dec=sol_dec)
    xx = sin(sol_dec) * sin(lat)
    yy = cos(sol_dec) * cos(lat)
    gsc = 4.92
    return 12 / pi * gsc * dr * ((omega2 - omega1) * xx + yy *
                                 (sin(omega2) - sin(omega1)))


def calc_res_surf(lai=None, r_s=70, r_l=100, lai_eff=0, srs=None, co2=None):
    """Surface resistance [s m-1].

    Parameters
    ----------
    lai: pandas.Series/float, optional
        leaf area index [-]
    r_s: pandas.series/float, optional
        surface resistance [s m-1]
    r_l: float, optional
        bulk stomatal resistance [s m-1]
    lai_eff: float, optional
        1 => LAI_eff = 0.5 * LAI
        2 => LAI_eff = lai / (0.3 * lai + 1.2)
        3 => LAI_eff = 0.5 * LAI; (LAI>4=4)
        4 => see [zhang_2008]_.
    srs: float, optional
        Relative sensitivity of rl to Δ[CO2] [yang_2019]_
    co2: float
        CO2 concentration [ppm]

    Returns
    -------
    pandas.Series containing the calculated surface resistance

    References
    -----
    .. [zhang_2008] Zhang, B., Kang, S., Li, F., & Zhang, L. (2008). Comparison
       of three evapotranspiration models to Bowen ratio-energy balance method
       for a vineyard in an arid desert region of northwest China. Agricultural
       and Forest Meteorology, 148(10), 1629-1640.
    .. [yang_2019] Yang, Y., Roderick, M. L., Zhang, S., McVicar, T. R., &
       Donohue, R. J. (2019). Hydrologic implications of vegetation response to
       elevated CO 2 in climate projections. Nature Climate Change, 9, 44-48.

    """
    if lai is None:
        return r_s
    else:
        fco2 = (1 + srs * (co2 - 300))
        return fco2 * r_l / calc_laieff(lai=lai, lai_eff=lai_eff)


def calc_laieff(lai=None, lai_eff=0):
    """Effective leaf area index [-].

    Parameters
    ----------
    lai: pandas.Series/float, optional
        leaf area index [-]
    lai_eff: float, optional
        0 => LAI_eff = 0.5 * LAI
        1 => LAI_eff = lai / (0.3 * lai + 1.2)
        2 => LAI_eff = 0.5 * LAI; (LAI>4=4)
        3 => see [zhang_2008]_.

    Returns
    -------
    pandas.Series containing the calculated effective leaf area index
    """
    if lai_eff == 0:
        return 0.5 * lai
    if lai_eff == 1:
        return lai / (0.3 * lai + 1.2)
    if lai_eff == 2:
        laie = lai.copy()
        laie[(lai > 2) & (lai < 4)] = 2
        laie[lai > 4] = 0.5 * lai
        return laie
    if lai_eff == 3:
        laie = lai.copy()
        laie[lai > 4] = 4
        return laie * 0.5


def calc_res_aero(wind, croph=None, zw=2, zh=2, ra_method=1):
    """Aerodynamic resistance [s m-1].

    Parameters
    ----------
    wind: pandas.Series
        mean day wind speed [m/s]
    croph: pandas.series/float, optional
        crop height [m]
    zw: float, optional
        height of wind measurement [m]
    zh: float, optional
         height of humidity and or air temperature measurement [m]
    ra_method: float, optional
        1 => ra = 208/wind
        2 => ra is calculated based on equation 36 in FAO (1990), ANNEX V.
    Returns
    -------
    pandas.Series containing the calculated aerodynamic resistance
    """
    if ra_method == 1:
        return 208 / wind
    else:
        d = 0.667 * croph
        zom = 0.123 * croph
        zoh = 0.0123 * croph
        return (log((zw - d) / zom)) * (log((zh - d) / zoh) /
                                        (0.41 ** 2) / wind)
