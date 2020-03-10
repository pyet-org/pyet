import numpy as np
import pandas as pd


def pm_hydrus(wind, elevation, latitude, rs=None, tmax=None, tmin=None,
              rhmax=None, rhmin=None, rh=None, croph=None):
    ta = (tmax + tmin) / 2
    lambd = lambda_calc(ta)
    pressure = press_calc(elevation)
    gamma = psy_calc(pressure, lambd)
    eamax = ea_calc(tmax)
    eamin = ea_calc(tmin)
    dlt = vpc_calc(ta, ea_calc(ta))
    ra = calc_ra(wind, croph, method=2)
    lai = lai_calc(croph=croph)
    rc = rc_calc(lai=lai)
    gamma1 = gamma1_calc(gamma)  # * (1 + rc / raa)
    ea = (eamax + eamin) / 2
    ed = ed_calc(tmax, tmin, rh=rh)

    gm_dl = gamma / (dlt + gamma1)
    AeroTCff = 0.622 * 3.486 * 86400. / ra / 1.01

    etaero = gm_dl * AeroTCff / (ta + 273.) * (ea - ed)

    dl_dl = dlt / (dlt + gamma1)
    s_flux = 0

    rns = calc_rns(rs=rs)  # in #  [MJ/m2/d]
    rnl = calc_rnl(rs, ed, latitude, tmax=tmax, tmin=tmin, rhmax=rhmax,
                   rhmin=rhmin, rh=rh)  # in #  [MJ/m2/d]
    rn = rns - rnl

    radterm = dl_dl * (rn - s_flux) / lambd
    pm = (etaero + radterm) * (
            1. - 7.37e-6 * (ta - 4.) ** 2 + 3.79e-8 * (ta - 4.) ** 3)
    return pm


def penman_monteith(wind, elevation, latitude, rs=None, tmax=None, tmin=None,
                    rhmax=None, rhmin=None, rh=None, croph=None):
    """Returns evapotranspiration calculated with the Penman-Monteith
    (FAO,1990) method.
    Parameters
    ----------
    wind: pandas.Series
        mean day wind speed [m/s]
    elevation: float/int
        the site elevation [m]
    latitude: float/int
        the site latitude [rad]
    rs: pandas.Series
        incoming measured solar radiation [MJ m-2 d-1]
    tmax: pandas.Series
        maximum day temperature [°C]
    tmin: pandas.Series
        minimum day temperature [°C]
    rhmax: pandas.Series
        maximum daily relative humidity [%]
    rhmin: pandas.Series
        mainimum daily relative humidity [%]
    rh: pandas.Series
        mean daily relative humidity [%]
    croph: float/int
        crop height [m]
    Returns
    -------
        pandas.Series containing the calculated evapotranspiration
    Examples
    --------
    >>> pm = et.penman_monteith(wind, elevation, latitude, rs=solar, tmax=tmax,
                                tmin=tmin, rh=rh, croph=0.6)
    """
    ta = (tmax + tmin) / 2
    lambd = lambda_calc(ta)
    pressure = press_calc(elevation)
    gamma = psy_calc(pressure, lambd)
    eamax = ea_calc(tmax)
    eamin = ea_calc(tmin)
    dlt = vpc_calc(ta, ea_calc(ta))
    ra = calc_ra(wind, croph, method=2)
    lai = lai_calc(croph=croph)
    rc = rc_calc(lai=lai)
    gamma1 = gamma1_calc(gamma)  # * (1 + rc / raa)
    ea = (eamax + eamin) / 2
    ed = ed_calc(tmax, tmin, rh=rh)
    cp = gamma * 0.622 * 2.45 / pressure
    qa = pressure / (1.01 * (ta + 273) * 0.287)
    rns = calc_rns(rs=rs)  # in #  [MJ/m2/d]
    rnl = calc_rnl(rs, ed, latitude, tmax=tmax, tmin=tmin, rhmax=rhmax,
                   rhmin=rhmin, rh=rh)  # in #  [MJ/m2/d]

    rn = rns - rnl
    den = (lambd * (dlt + gamma1))
    num1 = dlt * rn / den
    num2 = 86400 * cp * qa * (ea - ed) / ra / den
    pet = (num1 + num2)
    return pet


def calc_ra(wind=None, croph=None, method=1, ):
    if method == 1:
        return 208 / wind
    elif method == 2:
        return (np.log((2 - 0.667 * croph) / (0.123 * croph))) * \
               (np.log((2 - 0.667 * croph) / (0.0123 * croph))) / \
               (0.41 ** 2) / wind


def lai_calc(method=1, croph=None):
    if method == 1:
        return 0.24 * croph


def rc_calc(lai=None, method=1):
    if method == 1:
        return 200 / lai


def gamma1_calc(gamma, method=1):
    if method == 1:
        return gamma


def rn_calc(rs, ed, lat, tmax=None, tmin=None, rhmax=None,
            rhmin=None, rh=None):
    """
    Saturation Vapour Pressure  (ea)
    From FAO (1990), ANNEX V, eq. 50
    """
    rns = calc_rns(rs)
    rnl = calc_rnl(rs, ed, lat, tmax=tmax, tmin=tmin, rhmax=rhmax,
                   rhmin=rhmin, rh=rh)
    return rns - rnl


def calc_rnl1(rs, ed, lat, tmax=None, tmin=None, rhmax=None, rhmin=None,
              rh=None):
    """
    Net Shortwave Radiation Rns
    From FAO (1990), ANNEX V, eq. 56
    """
    steff = 4.9 * 10 ** (-9)  # MJm-2K-4d-1

    rso_0 = calc_rns(rs=rs_calc(rs.index, lat), alpha=0)
    ta = (tmax + tmin) / 2
    # sunshine_hours =
    # j = day_of_year(meteoindex)
    # sol_dec = solar_declination(j)
    # omega = sunset_angle(lat, sol_dec)
    # max_daylight = 24 * omega / np.pi
    # c_factor = cloudiness_factor(sunshine_hours, max_daylight)
    rso = rs_calc(rs.index, lat)
    c_factor = cloudiness_factor1(rs, rso)
    emiss = emissivity(ed)

    t_factor = ((tmax + 273.2) ** 4 + (tmin + 273) ** 4) / 2
    return steff * emiss * c_factor * (ta + 273) ** 4


def calc_rnl(rs, ed, lat, tmax=None, tmin=None, rhmax=None, rhmin=None,
             rh=None):
    """
    Net Shortwave Radiation Rns
    From FAO (1990), ANNEX V, eq. 56
    """
    steff = 4.9 * 10 ** (-9)  # MJm-2K-4d-1

    rso_0 = calc_rns(rs=rs_calc(rs.index, lat), alpha=0)

    # sunshine_hours =
    # j = day_of_year(meteoindex)
    # sol_dec = solar_declination(j)
    # omega = sunset_angle(lat, sol_dec)
    # max_daylight = 24 * omega / np.pi
    # c_factor = cloudiness_factor(sunshine_hours, max_daylight)
    rso = rs_calc(rs.index, lat)
    c_factor = cloudiness_factor1(rs, rso)
    emiss = emissivity(ed)

    t_factor = ((tmax + 273.2) ** 4 + (tmin + 273.2) ** 4) / 2
    return steff * emiss * c_factor * t_factor


def ed_calc(tmax, tmin, rhmax=None, rhmin=None, rh=None, ta=None):
    """
    Actual Vapour Pressure (ed)
    From FAO (1990), ANNEX V, eq. 11
    """
    eamax = ea_calc(tmax)
    eamin = ea_calc(tmin)
    if rhmax is not None:
        return (eamin * rhmax / 200) + (eamax * rhmin / 200)
    elif rh is not None:
        return rh / (50 / eamin + 50 / eamax)
    else:
        # Based on equation 19, page 39 in Allen et al (1998).
        return (rh / 100) * ea_calc(ta)


def cloudiness_factor(sunshine_hours, max_daylight, al=0.9, bl=0.1):
    """
    Cloudiness factor f
    From FAO (1990), ANNEX V, eq. 57
    """

    return al * sunshine_hours / max_daylight + bl


def cloudiness_factor1(rs, rso, ac=1.35, bc=-0.35):
    """
    Cloudiness factor f
    From FAO (1990), ANNEX V, eq. 57
    """

    return ac * rs / rso + bc


def relative_distance(j):
    """
    Relative distance Earth - Sun
    From FAO (1990), ANNEX V, eq. 21
    """
    return 1 + 0.033 * np.cos(2 * np.pi / 365 * j)


def day_of_year(meteoindex):
    return pd.to_numeric(meteoindex.strftime('%j'))


def sunset_angle(lat, sol_dec):
    """
    Sunset hour angle [rad] (omega)
    From FAO (1990), ANNEX V, eq. 20
    """
    return np.arccos(-np.tan(lat) * np.tan(sol_dec))


def solar_declination(j):
    """
    Solar declination [rad] (sol_dec)
    From FAO (1990), ANNEX V, eq. 22
    """
    return 0.4093 * np.sin(2 * np.pi / 365 * j - 1.39)


def lambda_calc(temperature):
    """
    From FAO (1990), ANNEX V, eq. 1
    """
    return 2.501 - 0.002361 * temperature


def psy_calc(pressure, lambd):
    """
    From FAO (1990), ANNEX V, eq. 4
    """
    return 0.0016286 * pressure / lambd


def press_calc(elevation):
    """
    From FAO (1990), ANNEX V, eq. 6
    """
    return 101.3 * ((293. - 0.0065 * elevation) / 293.) ** 5.253


def vpc_calc(temperature, ea):
    """
    From FAO (1990), ANNEX V, eq. 3
    """
    return 4098 * ea / (temperature + 237.3) ** 2


def vpc_calc1(tmax, tmin):
    """
    From FAO (1990), ANNEX V, eq. 3
    """
    eamax = ea_calc(tmax)
    eamin = ea_calc(tmin)
    return 2049 * eamax / (tmax + 237.3) ** 2 + \
           2049 * eamin / (tmin + 237.3) ** 2


def ea_calc(temperature):
    """
    Saturation Vapour Pressure  (ea)
    From FAO (1990), ANNEX V, eq. 10
    """
    return 0.6108 * np.exp((17.27 * temperature) / (temperature + 237.3))


def emissivity(ed, al=0.34, bl=-0.14):
    """
    Net Emissivity
    From FAO (1990), ANNEX V, eq. 60
    """
    return al + bl * np.sqrt(ed)


def calc_rns(rs=None, meteoindex=None, lat=None, alpha=0.23):
    """
    Net Shortwave Radiation Rns
    From FAO (1990), ANNEX V, eq. 51
    """
    if rs is not None:
        return (1 - alpha) * rs
    else:
        return (1 - alpha) * rs_calc(meteoindex, lat)


def rs_calc(meteoindex, lat, a_s=0.25, b_s=0.5):
    """
    Nncoming solar radiation rs
    From FAO (1990), ANNEX V, eq. 52
    """
    ra = ra_calc(meteoindex, lat)
    nn = 1
    return (a_s + b_s * nn) * ra


def ra_calc(meteoindex, lat):
    """
    Extraterrestrial Radiation (Ra)
    From FAO (1990), ANNEX V, eq. 18
    """

    j = day_of_year(meteoindex)
    dr = relative_distance(j)
    sol_dec = solar_declination(j)

    omega = sunset_angle(lat, sol_dec)
    gsc = 0.082 * 24 * 60
    # gsc = 1360
    return gsc / np.pi * dr * (omega * np.sin(sol_dec) * np.sin(lat) +
                               np.cos(sol_dec) * np.cos(lat) * np.sin(omega))
