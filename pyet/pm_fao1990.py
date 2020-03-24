import numpy as np
from pandas import to_numeric

pi=3.141592654

def pm_fao1990(wind, elevation, latitude, solar=None, tmax=None, tmin=None,
               rh=None, croph=None):
    """Returns evapotranspiration calculated with the FAO Penman-Monteith
    (Monteith, 1965; FAO, 1990) method.

    Based on equation 30 (FAO, 1990).

    Parameters
    ----------
    wind: pandas.Series
        mean day wind speed [m/s]
    elevation: float/int
        the site elevation [m]
    latitude: float/int
        the site latitude [rad]
    solar: pandas.Series
        incoming measured solar radiation [MJ m-2 d-1]
    tmax: pandas.Series
        maximum day temperature [°C]
    tmin: pandas.Series
        minimum day temperature [°C]
    rh: pandas.Series
        mean daily relative humidity [%]
    croph: float/int/pandas.series
        crop height [m]
    Returns
    -------
        pandas.Series containing the calculated evapotranspiration
    Examples
    --------
    >>> pm_fao1990_et = pm_fao1990(wind, elevation, latitude, solar=solar, \
                                   tmax=tmax, tmin=tmin, rh=rh, croph=0.6)
    """
    # aeroterm
    ta = (round(tmax, 8)+round(tmin, 8))/2
    lambd = round(lambda_calc(ta), 8)
    pressure = press_calc(elevation)
    gamma = psy_calc(pressure, lambd)
    eamax = e0_calc(tmax)
    eamin = e0_calc(tmin)
    dlt = vpc_calc(tmax, tmin)

    aerdyn = calc_raa(croph=croph, method=2)
    raa = aerdyn / wind

    gamma1 = gamma * (1 + 60. / raa)
    eamean = (eamax + eamin) / 2
    eadew = ed_calc(tmax, tmin, rh)

    gm_dl = gamma / (dlt + gamma1)
    aerotcff = 0.622 * 3.486 * 86400. / aerdyn / 1.01

    etaero = gm_dl * aerotcff / (ta + 273.) * wind * (eamean - eadew)

    dl_dl = dlt / (dlt + gamma1)
    # rad term
    rns = calc_rns(solar=solar)  # in #  [MJ/m2/d]

    rso = rs_calc(solar.index, latitude)  # radiation of clear sky
    cloudf = cloudiness_factor(solar, rso)
    rnl = calc_rnl(tmax, tmin, eadew, cloudf)  # in [MJ/m2/d]
    #rnl = calc_rnl(ta, ta, (e0_calc(ta)*rh/100), cloudf)  # in [MJ/m2/d]
    rn = rns - rnl

    radterm = dl_dl * (rn-0) / lambd
    pm = (etaero + radterm) * \
         (1. - 7.37e-6 * (ta - 4.) ** 2 + 3.79e-8 * (ta - 4.) ** 3)
    return rnl, rns, radterm, etaero, pm



def calc_rnl(tmax, tmin, ea, cloudf, longa=0.34, longb=-0.139):
    """
    Net Longwave Radiation Rnl
    From FAO (1990), ANNEX V, eq. 56
    Parameters
    ----------
    solar: Series
        incoming measured solar radiation [MJ m-2 d-1]
    ed: Series
        Actual Vapour Pressure (ed).
    lat: float/int
        the site latitude [rad]
    tmax: Series
        maximum day temperature [°C]
    tmin: Series
        minimum day temperature [°C]
    Returns
    -------
        Series containing the calculated net outgoing radiation
    """
    sigma = 0.00000000245*((tmax + 273.16) ** 4 + (tmin + 273.16) ** 4)
    emiss = longa + longb * round(np.sqrt(ea), 8)
    return sigma * cloudf * emiss


def cloudiness_factor(rs, rso, ac=1.35, bc=-0.35):
    """
    Cloudiness factor f
    From FAO (1990), ANNEX V, eq. 57
    """

    return ac * rs / rso + bc


def rs_calc(meteoindex, lat, a_s=0.25, b_s=0.5):
    """
    Nncoming solar radiation rs
    From FAO (1990), ANNEX V, eq. 52
    """
    ra = ra_calc(meteoindex, lat)
    nn = 1
    return (a_s + b_s * nn) * ra


def calc_raa(croph=None, method=1, windh=2, temph=2):
    if method == 1:
        return 208
    elif method == 2:
        return (np.log((windh - 0.667 * croph) / (0.123 * croph))) * \
               (np.log((temph - 0.667 * croph) / (0.0123 * croph))) / \
               (0.41 ** 2)


def lai_calc(method=1, croph=None):
    if method == 1:
        return 0.24 * croph


def rc_calc(lai=None, method=1):
    if method == 1:
        return 70
    elif method == 2:
        return 200 / lai


def ed_calc(tmax, tmin, rh):
    """
    Actual Vapour Pressure (ed).
    From FAO (1990), ANNEX V, eq. 11
    """
    eamax = e0_calc(tmax)
    eamin = e0_calc(tmin)
    return rh / (50. / eamin + 50. / eamax)


#def cloudiness_factor(sunshine_hours, max_daylight, al=0.9, bl=0.1):
#    """
#    Cloudiness factor f
#    From FAO (1990), ANNEX V, eq. 57
#    """
#    return al * sunshine_hours / max_daylight + bl


def relative_distance(j):
    """
    Relative distance Earth - Sun
    From FAO (1990), ANNEX V, eq. 21
    """
    return (1 + 0.033 * np.cos(2 * pi / 365 * j))


def day_of_year(meteoindex):
    return to_numeric(meteoindex.strftime('%j'))


def sunset_angle(lat, sol_dec):
    """
    Sunset hour angle [rad] (omega)
    From FAO (1990), ANNEX V, eq. 20
    """
    return np.arccos(np.tan(lat) * -np.tan(sol_dec))


def solar_declination(j):
    """
    Solar declination [rad] (sol_dec)
    From FAO (1990), ANNEX V, eq. 22
    """
    return 0.4093 * np.sin(2 * pi / 365 * j - 1.39)


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


def vpc_calc(tmax, tmin):
    """
    From FAO (1990), ANNEX V, eq. 3
    """
    eamax = e0_calc(tmax)
    eamin = e0_calc(tmin)
    return round((2049. * eamax / (tmax + 237.3) ** 2) +
                 (2049. * eamin / (tmin + 237.3) ** 2), 8)


def e0_calc(temperature):
    """
    Saturation Vapour Pressure  (ea)
    From FAO (1990), ANNEX V, eq. 10
    """
    return 0.6108 * np.exp((17.27 * temperature) / (temperature + 237.3))


def calc_rns(solar=None, meteoindex=None, lat=None, alpha=0.23):
    """
    Net Shortwave Radiation Rns
    From FAO (1990), ANNEX V, eq. 51
    """
    if solar is not None:
        return (1 - alpha) * solar
    else:
        return (1 - alpha) * rs_calc(meteoindex, lat)


def ra_calc(meteoindex, lat):
    """
    Extraterrestrial Radiation (Ra)
    From FAO (1990), ANNEX V, eq. 18
    """

    j = day_of_year(meteoindex)
    dr = relative_distance(j)
    sol_dec = solar_declination(j)

    omega = sunset_angle(lat, sol_dec)
    gsc = 0.082 * 24 * 60  # =118.08
    # gsc = 1360
    return gsc / pi * dr * (omega * np.sin(sol_dec) * np.sin(lat) +
                               np.cos(sol_dec) * np.cos(lat) * np.sin(omega))
