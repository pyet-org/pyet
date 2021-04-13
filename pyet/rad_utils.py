"""The rad_utils module contains utility funtions for radiation data

"""

from numpy import sqrt, log, clip

from .meteo_utils import calc_ea, extraterrestrial_r_hour, \
    extraterrestrial_r, daylight_hours

# Stefan Boltzmann constant - hourly [MJm-2K-4h-1]
STEFAN_BOLTZMANN_HOUR = 2.042 * 10 ** -10
# Stefan Boltzmann constant - daily [MJm-2K-4d-1]
STEFAN_BOLTZMANN_DAY = 4.903 * 10 ** -9


def calc_rad_long(rs, tmean=None, tmax=None, tmin=None, rhmax=None,
                  rhmin=None, rh=None, elevation=None, lat=None, rso=None,
                  a=1.35, b=-0.35, ea=None, freq="D"):
    """Net longwave radiation [MJ m-2 d-1].

    Parameters
    ----------
    rs: pandas.Series, optional
        incoming solar radiation [MJ m-2 d-1]
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
    elevation: float, optional
        the site elevation [m]
    lat: float, optional
        the site latitude [rad]
    rso: pandas.Series/float, optional
        clear-sky solar radiation [MJ m-2 day-1]
    a: float, optional
        empirical coefficient for Net Long-Wave radiation [-]
    b: float, optional
        empirical coefficient for Net Long-Wave radiation [-]
    ea: pandas.Series, optional
        actual vapor pressure [kPa]
    freq: string, optional
        "D" => daily estimation
        "H" => hourly estimation

    Returns
    -------
        pandas.Series containing the calculated net longwave radiation

    Notes
    -----
        Based on equation 39 in [allen_1998]_.

    References
    -----
    .. [allen_1998] Allen, R. G., Pereira, L. S., Raes, D., & Smith, M. (1998).
       Crop evapotranspiration-Guidelines for computing crop water
       requirements-FAO Irrigation and drainage paper 56. Fao, Rome, 300.
       (http://www.fao.org/3/x0490e/x0490e06.htm#TopOfPage)

    """
    if ea is None:
        ea = calc_ea(tmean=tmean, tmax=tmax, tmin=tmin, rhmax=rhmax,
                     rhmin=rhmin, rh=rh)

    if freq == "H":
        if rso is None:
            ra = extraterrestrial_r_hour(tindex=rs.index, lat=lat)
            rso = calc_rso(ra=ra, elevation=elevation)
        solar_rat = clip(rs / rso, 0.3, 1)
        tmp1 = STEFAN_BOLTZMANN_HOUR * (tmean + 273.2) ** 4
    else:
        if rso is None:
            ra = extraterrestrial_r(tindex=rs.index, lat=lat)
            rso = calc_rso(ra=ra, elevation=elevation)
        solar_rat = clip(rs / rso, 0.3, 1)
        if tmax is not None:
            tmp1 = STEFAN_BOLTZMANN_DAY * ((tmax + 273.2) ** 4 +
                                           (tmin + 273.2) ** 4) / 2
        else:
            tmp1 = STEFAN_BOLTZMANN_DAY * (tmean + 273.2) ** 4

    tmp2 = 0.34 - 0.139 * sqrt(ea)  # OK
    tmp2 = clip(tmp2, 0.05, 1)
    tmp3 = a * solar_rat + b  # OK
    return tmp1 * tmp2 * tmp3


def calc_rad_short(rs=None, tindex=None, lat=None, alpha=0.23, n=None,
                   nn=None):
    """Net shortwave radiation [MJ m-2 d-1].

    Parameters
    ----------
    rs: pandas.Series, optional
        incoming solar radiation [MJ m-2 d-1]
    tindex: pandas..DatetimeIndex
    lat: float, optional
        the site latitude [rad]
    alpha: float, optional
        surface albedo [-]
    n: pandas.Series/float, optional
        actual duration of sunshine [hour]
    nn: pandas.Series/float, optional
        maximum possible duration of sunshine or daylight hours [hour]


    Returns
    -------
        pandas.Series containing the calculated net shortwave radiation

    Notes
    -----
        Based on equation 38 in [allen_1998]_.

    References
    -----
    .. [allen_1998] Allen, R. G., Pereira, L. S., Raes, D., & Smith, M. (1998).
       Crop evapotranspiration-Guidelines for computing crop water
       requirements-FAO Irrigation and drainage paper 56. Fao, Rome, 300.
       (http://www.fao.org/3/x0490e/x0490e06.htm#TopOfPage)

    """
    if rs is not None:
        return (1 - alpha) * rs
    else:
        return (1 - alpha) * calc_rad_sol_in(tindex, lat, n=n, nn=nn)


def calc_rad_sol_in(tindex, lat, as1=0.25, bs1=0.5, n=None, nn=None):
    """Incoming solar radiation [MJ m-2 d-1].

    Parameters
    ----------
    tindex: pandas.DatetimeIndex
    lat: float, optional
        the site latitude [rad]
    as1: float, optional
        regression constant,  expressing the fraction of extraterrestrial
        reaching the earth on overcast days (n = 0) [-]
    bs1: float, optional
        empirical coefficient for extraterrestrial radiation [-]
    n: pandas.Series/float, optional
        actual duration of sunshine [hour]
    nn: pandas.Series/float, optional
        maximum possible duration of sunshine or daylight hours [hour]

    Returns
    -------
        pandas.Series containing the calculated net shortwave radiation

    Notes
    -----
        Based on equation 35 in [allen_1998]_.

    References
    -----
    .. [allen_1998] Allen, R. G., Pereira, L. S., Raes, D., & Smith, M. (1998).
       Crop evapotranspiration-Guidelines for computing crop water
       requirements-FAO Irrigation and drainage paper 56. Fao, Rome, 300.
       (http://www.fao.org/3/x0490e/x0490e06.htm#TopOfPage)

    """
    ra = extraterrestrial_r(tindex, lat)
    if n is None:
        n = daylight_hours(tindex, lat)
    return (as1 + bs1 * n / nn) * ra


def calc_rso(ra, elevation):
    """Clear-sky solar radiation [MJ m-2 day-1].

    Parameters
    ----------
    ra: pandas.Series, optional
        Extraterrestrial daily radiation [MJ m-2 d-1]
    elevation: float, optional
        the site elevation [m]

    Returns
    -------
        pandas.Series containing the calculated Clear-sky solar radiation

    Notes
    -----
        Based on equation 37 in [allen_1998]_.

    """
    return (0.75 + (2 * 10 ** -5) * elevation) * ra


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
        Relative sensitivity of rl to Δ[CO2]
    co2: float
        CO2 concentration [ppm]

    Returns
    -------
        pandas.Series containing the calculated surface resistance

    References
    -----
    .. [allen_1998] Allen, R. G., Pereira, L. S., Raes, D., & Smith, M. (1998).
       Crop evapotranspiration-Guidelines for computing crop water
       requirements-FAO Irrigation and drainage paper 56. Fao, Rome, 300.
       (http://www.fao.org/3/x0490e/x0490e06.htm#TopOfPage)
    .. [zhang_2008] Zhang, B., Kang, S., Li, F., & Zhang, L. (2008). Comparison
       of three evapotranspiration models to Bowen ratio-energy balance method
       for a vineyard in an arid desert region of northwest China. Agricultural
        and Forest Meteorology, 148(10), 1629-1640.
    .. [schymanski_2016] Schymanski, S. J., & Or, D. (2017). Leaf-scale
       experiments reveal an important omission in the Penman–Monteith
       equation. Hydrology and Earth System Sciences, 21(2), 685-706.
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

    References
    -----
    .. [allen_1998] Allen, R. G., Pereira, L. S., Raes, D., & Smith, M. (1998).
       Crop evapotranspiration-Guidelines for computing crop water
       requirements-FAO Irrigation and drainage paper 56. Fao, Rome, 300.
       (http://www.fao.org/3/x0490e/x0490e06.htm#TopOfPage)
    .. [zhang_2008] Zhang, B., Kang, S., Li, F., & Zhang, L. (2008). Comparison
       of three evapotranspiration models to Bowen ratio-energy balance method
       for a vineyard in an arid desert region of northwest China. Agricultural
        and Forest Meteorology, 148(10), 1629-1640.
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

    References
    -----
    .. [allen_1998] Allen, R. G., Pereira, L. S., Raes, D., & Smith, M. (1998).
       Crop evapotranspiration-Guidelines for computing crop water
       requirements-FAO Irrigation and drainage paper 56. Fao, Rome, 300.
       (http://www.fao.org/3/x0490e/x0490e06.htm#TopOfPage)

    """
    if ra_method == 1:
        return 208 / wind
    else:
        d = 0.667 * croph
        zom = 0.123 * croph
        zoh = 0.0123 * croph
        return (log((zw - d) / zom)) * (log((zh - d) / zoh) /
                                        (0.41 ** 2) / wind)
