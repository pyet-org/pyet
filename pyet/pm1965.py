import numpy as np
import pandas as pd


def pm1965(wind, elevation, latitude, solar=None, net=None, sflux=0, tmax=None,
           tmin=None, rhmax=None, rhmin=None, rh=None, n=None, nn=None,
           rso=None):
    """Returns evapotranspiration calculated with the FAO Penman-Monteith
    (Monteith, 1965; FAO, 1990) method.
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
    net: pandas.Series
        net radiation [MJ m-2 d-1]
    sflux: Series/float/int
        soil heat flux [MJ m-2 d-1]
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
    n: Series/float
        actual duration of sunshine [hour]
    nn: Series/float
        maximum possible duration of sunshine or daylight hours [hour]
    rso: Series/float
        clear-sky solar radiation [MJ m-2 day-1]
    Returns
    -------
        pandas.Series containing the calculated evapotranspiration
    Examples
    --------
    >>> pm1965 = pm1965(wind, elevation, latitude, rs=solar, tmax=tmax, \
                                tmin=tmin, rh=rh)
    """
    ta = (tmax + tmin) / 2
    lambd = lambda_calc(ta)
    pressure = press_calc(elevation)
    gamma = psy_calc(pressure)
    dlt = vpc_calc(ta)
    cp = 1.01  # [Jkg-1°C-1]
    rho_a = calc_rhoa(pressure, ta)
    ra = calc_ra(wind, method=1)
    rc = rc_calc(method=1)
    gamma1 = gamma * (1 + rc / ra)

    ea = ea_calc(tmax=tmax, tmin=tmin, rhmax=rhmax, rhmin=rhmin, rh=rh)
    es = es_calc(tmax, tmin)
    if net is None:
        rns = shortwave_r(solar=solar, n=n, nn=nn)
        rnl = longwave_r(solar=solar, tmax=tmax, tmin=tmin, rhmax=rhmax,
                         rhmin=rhmin, rh=rh, rso=rso, elevation=elevation,
                         lat=latitude, ea=ea)
        net = rns - rnl

    den = (lambd * (dlt + gamma1))
    num1 = (dlt * (net - sflux) / den)
    num2 = (gamma * (es - ea) * rho_a * cp / den)
    pet = (num1 + num2)
    return pet


def calc_rhoa(pressure, ta):
    r = 287  # [Jkg-1K-1] universal gas constant for dry air
    return pressure / (1.01 * (ta + 273) * r)


def calc_ra(wind=None, croph=None, method=1):
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
        return 70
    elif method == 2:
        return 200 / lai


def lambda_calc(temperature):
    """
    From FAO (1990), ANNEX V, eq. 1
    """
    return 2.501 - 0.002361 * temperature


def vpc_calc(temperature):
    """
    Slope of saturation vapour pressure curve at air Temperature.

    Based on equation 13. in Allen et al 1998.
    The slope of the vapour pressure curve is in the FAO-56 method calculated
    using mean air temperature
    Parameters
    ----------
    temperature: Series
        mean day temperature [degC]
    Returns
    -------
        Series of Saturation vapour pressure [kPa degC-1]

    """
    ea = e0_calc(temperature)
    return 4098 * ea / (temperature + 237.3) ** 2


def e0_calc(temperature):
    """
    saturation vapour pressure at the air temperature T.

    Based on equations 11 in ALLen et al (1998).
    Parameters
    ----------Saturation Vapour Pressure  (es) from air temperature
    temperature: pandas.Series
         temperature [degC]
    Returns
    -------
        pandas.Series of saturation vapour pressure at the air temperature
        T [kPa]

    """
    return 0.6108 * np.exp((17.27 * temperature) / (temperature + 237.3))


def es_calc(tmax, tmin):
    """
    saturation vapour pressure at the air temperature T.

    Based on equations 11 in ALLen et al (1998).
    Parameters
    ----------Saturation Vapour Pressure  (es) from air temperature
    tmax: pandas.Series
        maximum day temperature [°C]
    tmin: pandas.Series
        minimum day temperature [°C]
    Returns
    -------
        pandas.Series of saturation vapour pressure at the air temperature
        T [kPa]

    """
    eamax = e0_calc(tmax)
    eamin = e0_calc(tmin)
    return (eamax + eamin) / 2


def ea_calc(tmax, tmin, rhmax=None, rhmin=None, rh=None):
    """Actual Vapour Pressure (ea) from air temperature.

    Based on equations 17, 18, 19, in ALLen et al (1998).
    Parameters
    ----------
    tmax: Series
        maximum day temperature [degC]
    tmin: Series
        minimum day temperature [degC]
    rhmax: Series
        maximum daily relative humidity [%]
    rhmin: Series
        mainimum daily relative humidity [%]
    rh: pandas.Series/int
        mean daily relative humidity [%]
    Returns
    -------
        Series of saturation vapour pressure at the air temperature
        T [kPa]
    """
    eamax = e0_calc(tmax)
    eamin = e0_calc(tmin)
    if rhmax is not None and rhmin is not None:  # eq. 17
        return (eamin * rhmax / 200) + (eamax * rhmin / 200)
    elif rhmax is not None and rhmin is None:  # eq.18
        return eamin * rhmax / 100
    elif rhmax is None and rhmin is not None:  # eq. 48
        return eamin
    elif rh is not None:  # eq. 19
        return rh / 200 * (eamax + eamin)
    else:
        print("error")


def day_of_year(meteoindex):
    """
    Return day of the year (1-365) based on pandas.series.index

    Parameters
    ----------
    meteoindex: pandas.Series.index

    Returns
    -------
        array of with ints specifyng day of year.

    """
    return pd.to_numeric(meteoindex.strftime('%j'))


def rso_calc(ra, elevation):
    """
    Actual Vapour Pressure (ea) from air temperature.

    Based on equations 37 in ALLen et al (1998).
    Parameters
    ----------
    ra: Series
        extraterrestrial radiation [MJ m-2 day-1]
    elevation: float/int
        the site elevation [m]
    Returns
    -------
        Series of clear-sky solar radiation [MJ m-2 day-1]

    """
    return (0.75 + (2 * 10 ** -5) * elevation) * ra


def relative_distance(j):
    """
    Inverse relative distance between earth and sun from day of the year.

    Based on equation 23 in Allen et al (1998).

    Parameters
    ----------
    j: array.py
        day of the year (1-365)
    Returns
    -------
        array.py specifyng day of year.
    """
    return 1 + 0.033 * np.cos(2 * np.pi / 365 * j)


def sunset_angle(sol_dec, lat):
    """
    Sunset hour angle from latitude and solar declination [rad].

    Based on equations 25 in ALLen et al (1998).
    Parameters
    ----------
    sol_dec: Series/float
        solar declination [rad]
    lat: float/int
        the site latitude [rad]
    Returns
    -------
        Series of sunset hour angle [rad].

    """
    return np.arccos(-np.tan(lat) * np.tan(sol_dec))


def solar_declination(j):
    """
    Solar declination [rad] from day of year [rad].

    Based on equations 24 in ALLen et al (1998).
    Parameters
    ----------
    j: array.py
        day of the year (1-365)
    Returns
    -------
        array.py of solar declination [rad].

    """
    return 0.4093 * np.sin(2 * np.pi / 365 * j - 1.39)


def extraterrestrial_r(meteoindex, lat):
    """
    Extraterrestrial Radiation (Ra)

    Based on equation 21 in Allen et al (1998).
    Parameters
    ----------
    meteoindex: pandas.Series.index
    lat: float/int
        the site latitude [rad]
    Returns
    -------
        array of ints of solar declination [rad].

    """
    j = day_of_year(meteoindex)
    dr = relative_distance(j)
    sol_dec = solar_declination(j)

    omega = sunset_angle(lat, sol_dec)
    gsc = 0.082  # solar constant [ MJ m-2 min-1]
    return gsc * 24 * 60 / np.pi * dr * (
            omega * np.sin(sol_dec) * np.sin(lat) +
            np.cos(sol_dec) * np.cos(lat) * np.sin(omega))


def psy_calc(pressure):
    """
    Psychrometric constant [kPa degC-1].

    Based on equation 8 in Allen et al (1998).
    Parameters
    ----------
    pressure: int/real
        atmospheric pressure [kPa].
    Returns
    -------
        pandas.series of Psychrometric constant [kPa degC-1].

    """
    return 0.000665 * pressure


def press_calc(elevation):
    """
    Atmospheric pressure.

    Based on equation 7 in Allen et al (1998).
    Parameters
    ----------
    elevation: int/real
        elevation above sea level [m].
    Returns
    -------
        int/real of atmospheric pressure [kPa].

    """
    return 101.3 * ((293. - 0.0065 * elevation) / 293.) ** 5.26


def longwave_r(solar, tmax=None, tmin=None, rhmax=None, rhmin=None,
               rh=None, rso=None, elevation=None, lat=None, ea=None):
    """Net outgoing longwave radiation.

    Based on equation 39 in Allen et al (1998).
        Parameters
    ----------
    rs: pandas.Series
        incoming measured solar radiation [MJ m-2 d-1]
    elevation: float/int
        the site elevation [m]
    lat: float/int
        the site latitude [rad]
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
    Returns
    -------
        pandas.Series containing the calculated net outgoing radiation
    """
    steff = 4.903 * 10 ** (-9)  # MJm-2K-4d-1
    if rso is None:
        ra = extraterrestrial_r(solar.index, lat)
        rso = rso_calc(ra, elevation)
    solar_rat = solar / rso
    if ea is None:
        ea = ea_calc(tmin=tmin, tmax=tmax, rhmin=rhmin, rhmax=rhmax, rh=rh)
    tmp1 = steff * ((tmax + 273.2) ** 4 + (tmin + 273.2) ** 4) / 2
    tmp2 = 0.34 - 0.14 * np.sqrt(ea)
    tmp3 = 1.35 * solar_rat - 0.35
    return tmp1 * tmp2 * tmp3


def shortwave_r(solar=None, meteoindex=None, lat=None, alpha=0.23, n=None,
                nn=None):
    """
    Net solar or shortwave radiation

    Based on equation 38 in Allen et al (1998).
    Parameters
    ----------
    meteoindex: pandas.Series.index
    rs: Series
        incoming measured solar radiation [MJ m-2 d-1]
    lat: float/int
        the site latitude [rad]
    alpha: float/int
        albedo or canopy reflection coefficient, which is 0.23 for the
        hypothethical grass reference crop [-]
    n: float/int
        actual duration of sunshine [hour]
    nn: float/int
        daylight hours [-]
    Returns
    -------
        Series containing the calculated net outgoing radiation
    """
    if solar is not None:
        return (1 - alpha) * solar
    else:
        return (1 - alpha) * in_solar_r(meteoindex, lat, n=n, nn=nn)


def in_solar_r(meteoindex, lat, a_s=0.25, b_s=0.5, n=None, nn=None):
    """
    Nncoming solar radiation rs
    Based on eq. 35 from FAO56.
    """
    ra = extraterrestrial_r(meteoindex, lat)
    if n is None:
        n = daylight_hours(meteoindex, lat)
    return (a_s + b_s * n / nn) * ra


def daylight_hours(meteoindex, lat):
    """
    Daylight hours
    Based on eq. 34 from FAO56.
    """
    j = day_of_year(meteoindex)
    sol_dec = solar_declination(j)
    sangle = sunset_angle(sol_dec, lat)
    return round(24 / np.pi * sangle, 1)