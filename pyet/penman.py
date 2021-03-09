from numpy import sqrt, log, exp, clip

from .utils import extraterrestrial_r, daylight_hours, extraterrestrial_r_hour


def penman(wind, elevation, lat, solar=None, net=None, sflux=0, tmax=None,
           tmin=None, rhmax=None, rhmin=None, rh=None, n=None, nn=None,
           rso=None, a=2.6, b=0.536):
    """Evapotranspiration calculated with the Penman (1948) method.

    Parameters
    ----------
    wind: pandas.Series
        mean day wind speed [m/s]
    elevation: float/int
        the site elevation [m]
    lat: float/int
        the site latitude [rad]
    solar: pandas.Series, optional
        incoming measured solar radiation [MJ m-2 d-1]
    net: pandas.Series, optional
        net radiation [MJ m-2 d-1]
    sflux: pandas.Series/int, optional
        soil heat flux [MJ m-2 d-1]
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
    n: pandas.Series/float, optional
        actual duration of sunshine [hour]
    nn: pandas.Series/float, optional
        maximum possible duration of sunshine or daylight hours [hour]
    rso: pandas.Series/float, optional
        clear-sky solar radiation [MJ m-2 day-1]
    a: float/int, optional
        wind coefficient [-]
    b: float/int, optional
        wind coefficient [-]

    Returns
    -------
    pandas.Series
        Series containing the calculated evapotranspiration

    Examples
    --------
    >>> penman_et = penman(wind, elevation, lat, solar=solar, tmax=tmax,
    >>>                    tmin=tmin, rh=rh)

    Notes
    -----
    Based on equation 6 in Allen et al (1998).

    """
    ta = (tmax + tmin) / 2
    pressure = press_calc(elevation=elevation, temperature=ta)
    gamma = psy_calc(pressure=pressure)
    dlt = vpc_calc(temperature=ta)
    lambd = lambda_calc(temperature=ta)

    ea = ea_calc(tmax=tmax, tmin=tmin, rhmax=rhmax, rhmin=rhmin, rh=rh)
    es = es_calc(tmax=tmax, tmin=tmin)

    if net is None:
        rns = shortwave_r(solar=solar, n=n, nn=nn)  # in #  [MJ/m2/d]
        rnl = longwave_r(solar=solar, tmax=tmax, tmin=tmin, rhmax=rhmax,
                         rhmin=rhmin, rh=rh, rso=rso, elevation=elevation,
                         lat=lat, ea=ea)  # in #  [MJ/m2/d]
        net = rns - rnl

    w = a * (1 + b * wind)

    den = lambd * (dlt + gamma)
    num1 = (dlt * (net - sflux) / den)
    num2 = (gamma * (es - ea) * w / den)
    pet = (num1 + num2)
    return pet


def pm_fao56(wind, elevation, lat, solar=None, net=None, sflux=0, tmax=None,
             tmin=None, rhmax=None, rhmin=None, rh=None, n=None, nn=None,
             rso=None):
    """Reference evapotranspiration using the FAO-56 Penman-Monteith
    equation (Monteith, 1965; Allen et al, 1998).

    Parameters
    ----------
    wind: Series
        mean day wind speed [m/s]
    elevation: float/int
        the site elevation [m]
    lat: float/int
        the site latitude [rad]
    solar: pandas.Series, optional
        incoming measured solar radiation [MJ m-2 d-1]
    net: pandas.Series, optional
        net radiation [MJ m-2 d-1]
    sflux: Series/float/int, optional
        soil heat flux [MJ m-2 d-1]
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
    n: Series/float, optional
        actual duration of sunshine [hour]
    nn: Series/float, optional
        maximum possible duration of sunshine or daylight hours [hour]
    rso: Series/float, optional
        clear-sky solar radiation [MJ m-2 day-1]

    Returns
    -------
        pandas.Series containing the calculated evapotranspiration

    Examples
    --------
    >>> pm_fao56_et = pm_fao56(wind, elevation, lat, solar=solar,
    >>>                        tmax=tmax, tmin=tmin, rh=rh)

    Notes
    -----
    Based on equation 6 in Allen et al (1998).

    """
    ta = (tmax + tmin) / 2
    pressure = press_calc(elevation=elevation, temperature=ta)
    gamma = psy_calc(pressure=pressure)
    dlt = vpc_calc(temperature=ta)

    gamma1 = (gamma * (1 + 0.34 * wind))

    ea = ea_calc(tmax=tmax, tmin=tmin, rhmax=rhmax, rhmin=rhmin, rh=rh)
    es = es_calc(tmax=tmax, tmin=tmin)
    if net is None:
        rns = shortwave_r(solar=solar, n=n, nn=nn)  # in [MJ/m2/d]
        rnl = longwave_r(solar=solar, tmax=tmax, tmin=tmin, rhmax=rhmax,
                         rhmin=rhmin, rh=rh, rso=rso, elevation=elevation,
                         lat=lat, ea=ea)  # in [MJ/m2/d]
        net = rns - rnl

    den = (dlt + gamma1)
    num1 = (0.408 * dlt * (net - sflux))
    num2 = (gamma * (es - ea) * 900 * wind / (ta + 273))
    return (num1 + num2) / den


def pm_asce(wind, elevation, lat, solar=None, net=None, sflux=0, tmax=None,
            tmin=None, rhmax=None, rhmin=None, rh=None, n=None, nn=None,
            rso=None, lai=None, croph=None, rs=1, ra=1, rl=100):
    """
    Returns evapotranspiration calculated with the ASCE Penman-Monteith
    (Monteith, 1965; ASCE, 2005) method.

    Parameters
    ----------
    wind: pandas.Series
        mean day wind speed [m/s]
    elevation: float/int
        the site elevation [m]
    lat: float/int
        the site latitude [rad]
    solar: pandas.Series, optional
        incoming measured solar radiation [MJ m-2 d-1]
    net: pandas.Series, optional
        net radiation [MJ m-2 d-1]
    sflux: Series/float/int, optional
        soil heat flux [MJ m-2 d-1]
    tmax: pandas.Series, optional
        maximum day temperature [°C]
    tmin: pandas.Series, optional
        minimum day temperature [°C]
    rhmax: pandas.Series, optional
        maximum daily relative humidity [%]
    rhmin: pandas.Series, optional
        minimum daily relative humidity [%]
    rh: pandas.Series, optional
        mean daily relative humidity [%]
    n: pandas.Series/float, optional
        actual duration of sunshine [hour]
    nn: pandas.Series/float, optional
        maximum possible duration of sunshine or daylight hours [hour]
    rso: pandas.Series/float, optional
        clear-sky solar radiation [MJ m-2 day-1]
    lai: pandas.Series/float, optional
        measured leaf area index [-]
    croph: float/int/pandas.series, optional
        crop height [m]
    rs: int, optional
        1 => rs = 70
        2 => rs = rl/LAI; rl = 200
    ra: int, optional
        1 => ra = 208/wind
        2 => ra is calculated based on equation 36 in FAO (1990), ANNEX V.

    Returns
    -------
        pandas.Series containing the calculated evapotranspiration

    Examples
    --------
    >>> pmasce = pm_asce(wind, elevation, lat, rs=solar, tmax=tmax,
    >>>                  tmin=tmin, rh=rh)

    """
    ta = (tmax + tmin) / 2
    lambd = lambda_calc(ta)
    pressure = press_calc(elevation, ta)
    gamma = psy_calc(pressure)
    dlt = vpc_calc(ta)
    cp = 1.013 * 10 ** (-3)
    r_a = aero_r(wind, method=ra, croph=croph)
    r_s = surface_r(method=rs, lai=lai, rl=rl)
    gamma1 = gamma * (1 + r_s / r_a)

    ea = ea_calc(tmax=tmax, tmin=tmin, rhmax=rhmax, rhmin=rhmin, rh=rh)
    es = es_calc(tmax, tmin)
    rho_a = calc_rhoa(pressure, ta, ea)
    if net is None:
        rns = shortwave_r(solar=solar, n=n, nn=nn)
        rnl = longwave_r(solar=solar, tmax=tmax, tmin=tmin, rhmax=rhmax,
                         rhmin=rhmin, rh=rh, rso=rso, elevation=elevation,
                         lat=lat, ea=ea)
        net = rns - rnl
    kmin = 86400
    den = (lambd * (dlt + gamma1))
    num1 = (dlt * (net - sflux) / den)
    num2 = (rho_a * cp * kmin * (es - ea) / r_a / den)
    return num1 + num2


def pm_corrected(wind, elevation, lat, solar=None, net=None, sflux=0,
                 tmean=None, tmax=None, tmin=None, rhmax=None, rhmin=None,
                 rh=None, n=None, nn=None, rso=None, lai=None, croph=None,
                 r_s=None, rs=1, ra=1, a_s=1, a_sh=1, rl=100, a=1.35, b=-0.35,
                 co2=300, srs=0.0009, laieff=0, flai=1, freq="D"):
    """Evapotranspiration calculated with the upscaled corrected
    Penman-Monteith equation from Schymanski (Schymanski, 2017).

    Parameters
    ----------
    wind: pandas.Series
        mean day wind speed [m/s]
    elevation: float
        the site elevation [m]
    lat: float
        the site latitude [rad]
    solar: pandas.Series, optional
        incoming measured solar radiation [MJ m-2 d-1]
    net: pandas.Series, optional
        net radiation [MJ m-2 d-1]
    sflux: Series/float, optional
        soil heat flux [MJ m-2 d-1]
    tmax: pandas.Series, optional
        maximum day temperature [°C]
    tmin: pandas.Series, optional
        minimum day temperature [°C]
    rhmax: pandas.Series, optional
        maximum daily relative humidity [%]
    rhmin: pandas.Series, optional
        minimum daily relative humidity [%]
    rh: pandas.Series, optional
        mean daily relative humidity [%]
    n: pandas.Series/float, optional
        actual duration of sunshine [hour]
    nn: pandas.Series/float, optional
        maximum possible duration of sunshine or daylight hours [hour]
    rso: pandas.Series/float, optional
        clear-sky solar radiation [MJ m-2 day-1]
    lai: pandas.Series/float, optional
        measured leaf area index [-]
    croph: float/int/pandas.series, optional
        crop height [m]
    rs: int, optional
        1 => rs = 70
        2 => rs = rl/LAI; rl = 200
    ra: int, optional
        1 => ra = 208/wind
        2 => ra is calculated based on equation 36 in FAO (1990), ANNEX V.
    a_s: int, optional
        Fraction of one-sided leaf area covered by stomata (1 if stomata are 1
        on one side only, 2 if they are on both sides)
    a_sh: int, optional
        Fraction of projected area exchanging sensible heat with the air (2)
    a: float/int, optional
        wind coefficient [-]
    b: float/int, optional
        wind coefficient [-]

    Returns
    -------
    pandas.Series
        Series containing the calculated evapotranspiration.

    """
    if "D" in freq:
        tmean = (tmax + tmin) / 2
        es = es_calc(tmax=tmax, tmin=tmin)
    else:
        es = e0_calc(temperature=tmean)
    lambd = lambda_calc(temperature=tmean)
    pressure = press_calc(elevation=elevation, temperature=tmean)
    gamma = psy_calc(pressure=pressure)
    dlt = vpc_calc(temperature=tmean)
    r_a = aero_r(wind, method=ra, croph=croph)
    ea = ea_calc(tmax=tmax, tmin=tmin, rhmax=rhmax, rhmin=rhmin, rh=rh)
    rho_a = calc_rhoa(pressure=pressure, ta=tmean, ea=ea)

    if net is None:
        rns = shortwave_r(solar=solar, n=n, nn=nn)
        rnl = longwave_r(solar=solar, tmax=tmax, tmin=tmin, rhmax=rhmax,
                         rhmin=rhmin, rh=rh, rso=rso, elevation=elevation,
                         lat=lat, ea=ea, a=a, b=b, freq=freq)
        net = rns - rnl * a_sh

    kmin = 86400
    if freq == "H":
        kmin /= 24

    if r_s is None:
        r_s = surface_r(method=rs, lai=lai, rl=rl, co2=co2, srs=srs,
                        laieff=laieff, flai=flai)

    cp = 1.013 * 10 ** (-3)
    gamma1 = gamma * a_sh / a_s * (1 + r_s / r_a)
    den = (lambd * (dlt + gamma1))
    num1 = (dlt * (net - sflux) / den)
    num2 = (rho_a * cp * kmin * (es - ea) * a_sh / r_a / den)
    return num1 + num2


def pm_fao1990(wind, elevation, lat, solar=None, tmax=None, tmin=None,
               rh=None, croph=None, ra=2, n=None, nn=None):
    """Evapotranspiration calculated with the FAO Penman-Monteith
    (Monteith, 1965; FAO, 1990) method.

    Parameters
    ----------
    wind: pandas.Series
        mean day wind speed [m/s]
    elevation: float/int
        the site elevation [m]
    lat: float/int
        the site latitude [rad]
    solar: pandas.Series, optional
        incoming measured solar radiation [MJ m-2 d-1]
    tmax: pandas.Series, optional
        maximum day temperature [°C]
    tmin: pandas.Series, optional
        minimum day temperature [°C]
    rh: pandas.Series, optional
        mean daily relative humidity [%]
    croph: float/int/pandas.series, optional
        crop height [m]

    Returns
    -------
    pandas.Series
        Series containing the calculated evapotranspiration

    Examples
    --------
    >>> pm_fao1990_et = pm_fao1990(wind, elevation, lat, solar=solar,
    >>>                            tmax=tmax, tmin=tmin, rh=rh, croph=0.6)

    Notes
    -----
    Based on equation 30 (FAO, 1990).

    """
    # aeroterm
    ta = (tmax + tmin) / 2.
    lambd = lambda_calc(temperature=ta)
    pressure = press_calc(elevation=elevation, temperature=ta)
    gamma = psy_calc(pressure=pressure, lambd=lambd)
    eamax = e0_calc(temperature=tmax)
    eamin = e0_calc(temperature=tmin)

    raa = aero_r(wind, method=ra, croph=croph)
    eadew = ed_calc(tmax=tmax, tmin=tmin, rh=rh)  # OK
    aerodyn = raa * wind
    aerotcff = 0.622 * 3.486 * 86400. / aerodyn / 1.01
    lai = croph * 24
    rs = 200 / lai
    gamma1 = gamma * (1 + rs / raa)
    dlt = vpc_calc(tmin=tmin, tmax=tmax, method=1)

    gm_dl = gamma / (dlt + gamma1)
    eamean = (eamax + eamin) / 2
    etaero = gm_dl * aerotcff / (ta + 273.) * wind * (eamean - eadew)

    dl_dl = dlt / (dlt + gamma)
    # rad term
    rso = rs_calc(tindex=solar.index, lat=lat)  # OK
    rns = shortwave_r(solar=solar, n=n, nn=nn)
    rnl = longwave_r(solar, tmax=tmax, tmin=tmin, rh=rh, rso=rso,
                     elevation=elevation, lat=lat, ea=eadew)
    net = rns - rnl

    radterm = dl_dl * (net) / lambd
    pm = (etaero + radterm)
    return pm, radterm, etaero, rnl, rns


def priestley_taylor(wind, elevation, lat, solar=None, net=None, tmax=None,
                     tmin=None, rhmax=None, rhmin=None, rh=None, rso=None,
                     n=None, nn=None, alpha=1.26):
    """Evapotranspiration calculated with Penman-Monteith (FAO, 1990) method.

    Parameters
    ----------
    wind: pandas.Series
        mean day wind speed [m/s]
    elevation: float/int
        the site elevation [m]
    lat: float/int
        the site latitude [rad]
    solar: pandas.Series
        incoming measured solar radiation [MJ m-2 d-1]
    net: pandas.Series
        net radiation [MJ m-2 d-1]
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
    alpha: Series/float
        calibration coefficient

    Returns
    -------
        pandas.Series containing the calculated evapotranspiration

    Examples
    --------
    >>> pm = priestley_taylor(wind, elevation, lat, solar=solar,
    >>>                       tmax=tmax, tmin=tmin, rh=rh, croph=0.6)

    Notes
    -----
    Based on equation 6 in Allen et al (1998).

    """
    ta = (tmax + tmin) / 2
    lambd = lambda_calc(temperature=ta)
    pressure = press_calc(elevation=elevation, temperature=ta)
    gamma = psy_calc(pressure=pressure, lambd=None)
    dlt = vpc_calc(temperature=ta, method=0)

    ea = ea_calc(tmax, tmin, rhmax=rhmax, rhmin=rhmin, rh=rh)

    if net is None:
        rns = shortwave_r(solar=solar, n=n, nn=nn)  # in [MJ/m2/d]
        rnl = longwave_r(solar=solar, tmax=tmax, tmin=tmin, rhmax=rhmax,
                         rhmin=rhmin, rh=rh, rso=rso, elevation=elevation,
                         lat=lat, ea=ea)  # in [MJ/m2/d]
        net = rns - rnl

    return (alpha * dlt * net) / (lambd * (dlt + gamma))


def makkink(tmax, tmin, rs, elevation, f=1):
    """Evapotranspiration calculated with the Makkink (1957) method.

    Parameters
    ----------
    tmax: pandas.Series
        maximum day temperature [°C]
    tmin: pandas.Series
        minimum day temperature [°C]
    rs: pandas.Series
        incoming measured solar radiation [MJ m-2 d-1]
    elevation: float
        the site elevation [m]
    f: float, optional
        crop coefficient [-]

    Returns
    -------
        Series containing the calculated evapotranspiration

    Examples
    --------
    >>> mak = makkink(tmax, tmin, rs, elevation)

    """
    ta = (tmax + tmin) / 2
    pressure = press_calc(elevation, ta)
    gamma = psy_calc(pressure=pressure, lambd=None)
    dlt = vpc_calc(temperature=ta, method=0)

    return f / 2.45 * 0.61 * rs * dlt / (dlt + gamma) - 0.12


##% Utility functions (TODO: Make private?)


def longwave_r(solar, tmax=None, tmin=None, rhmax=None, rhmin=None,
               rh=None, rso=None, elevation=None, lat=None, ea=None,
               a=1.35, b=-0.35, freq="D", ta=None):
    """Net outgoing longwave radiation.

    Parameters
    ----------
    solar: Series
        incoming measured solar radiation [MJ m-2 d-1]
    elevation: float
        the site elevation above sea level [m]
    lat: float/int
        the site latitude [rad]
    tmax: Series
        maximum day temperature [°C]
    tmin: Series
        minimum day temperature [°C]
    rhmax: Series
        maximum daily relative humidity [%]
    rhmin: Series
        mainimum daily relative humidity [%]
    rh: Series
        mean daily relative humidity [%]
    rso: Series/float
        clear-sky solar radiation [MJ m-2 day-1]
    ea: Series
        actual vapour pressure.

    Returns
    -------
        pandas.Series containing the calculated net outgoing radiation

    Notes
    -----
    Based on equation 39 in Allen et al (1998).

    """
    if ea is None:
        ea = ea_calc(tmin=tmin, tmax=tmax, rhmin=rhmin, rhmax=rhmax, rh=rh)

    if freq == "H":
        steff = 2.042 * 10 ** (-10)  # MJm-2K-4h-1
        if rso is None:
            ra = extraterrestrial_r_hour(tindex=solar.index, lat=lat)
            rso = rso_calc(ra=ra, elevation=elevation)
        solar_rat = clip(solar / rso, 0.3, 1)
        tmp1 = steff * (ta + 273.2) ** 4
    else:
        steff = 4.903 * 10 ** (-9)  # MJm-2K-4d-1
        if rso is None:
            ra = extraterrestrial_r(tindex=solar.index, lat=lat)
            rso = rso_calc(ra=ra, elevation=elevation)
        solar_rat = clip(solar / rso, 0.3, 1)
        tmp1 = steff * ((tmax + 273.2) ** 4 + (tmin + 273.2) ** 4) / 2

    tmp2 = 0.34 - 0.139 * sqrt(ea)  # OK
    tmp2 = clip(tmp2, 0.05, 1)
    tmp3 = a * solar_rat + b  # OK
    return tmp1 * tmp2 * tmp3


def vpc_calc(temperature=None, tmin=None, tmax=None, method=0):
    """
    Slope of saturation vapour pressure curve at air Temperature.

    Parameters
    ----------
    temperature: Series
        mean day temperature [degC].

    Returns
    -------
        Series of Saturation vapour pressure [kPa degC-1].

    Notes
    -----
    if method is 0:
        Based on equation 13. in Allen et al 1998. The slope of the vapour
        pressure curve is in the FAO-56 method calculated using mean air
        temperature
    if method is 1:
        From FAO (1990), ANNEX V, eq. 3

    """
    if method == 0:
        ea = e0_calc(temperature)
        return 4098 * ea / (temperature + 237.3) ** 2
    elif method == 1:
        eamax = e0_calc(tmax)
        eamin = e0_calc(tmin)
        return round((2049. * eamax / (tmax + 237.3) ** 2) +
                     (2049. * eamin / (tmin + 237.3) ** 2), 8)
    elif method == 2:
        return 2503 * exp((17.27 * temperature) / (temperature + 237.3)) / (
                temperature + 237.3) ** 2


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
    return 0.6108 * exp((17.27 * temperature) / (temperature + 237.3))


def es_calc(tmax, tmin):
    """
    saturation vapour pressure at the air temperature T.

    Based on equations 11 in Allen et al (1998).

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


def rso_calc(ra, elevation):
    """
    Actual Vapour Pressure (ea) from air temperature.

    Based on equations 37 in ALLen et al (1998).
    Parameters
    ----------
    ra: Series
        extraterrestrial radiation [MJ m-2 day-1]
    elevation: float
        the site elevation above sea level [m]

    Returns
    -------
        Series of clear-sky solar radiation [MJ m-2 day-1]

    """
    return (0.75 + (2 * 10 ** -5) * elevation) * ra


def psy_calc(pressure, lambd=None):
    """Psychrometric constant [kPa degC-1].

    Parameters
    ----------
    pressure: int/real
        atmospheric pressure [kPa].
    lambd: float,m optional
        Divide the pressure by this value.

    Returns
    -------
        pandas.series of Psychrometric constant [kPa degC-1].

    Notes
    -----
    if lambd is none:
        From FAO (1990), ANNEX V, eq. 4
    else:
        Based on equation 8 in Allen et al (1998).

    """
    if lambd is None:
        return 0.000665 * pressure
    else:
        return 0.0016286 * pressure / lambd


def press_calc(elevation, temperature):
    """Atmospheric pressure. Based on equation 7 in Allen et al (1998).

    Parameters
    ----------
    elevation: int/real
        elevation above sea level [m].
    temperature

    Returns
    -------
        int/real of atmospheric pressure [kPa].

    """
    return 101.3 * (((273.16 + temperature) - 0.0065 * elevation) /
                    (273.16 + temperature)) ** (9.807 / (0.0065 * 287))


def shortwave_r(solar=None, tindex=None, lat=None, alpha=0.23, n=None,
                nn=None):
    """Net solar or shortwave radiation.

    Based on equation 38 in Allen et al (1998).

    Parameters
    ----------
    tindex: pandas.Series.index
    solar: Series
        incoming measured solar radiation [MJ m-2 d-1]
    lat: float/int
        the site latitude [rad]
    alpha: float/int
        albedo or canopy reflection coefficient, which is 0.23 for the
        hypothetical grass reference crop [-]
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
        return (1 - alpha) * in_solar_r(tindex, lat, n=n, nn=nn)


def in_solar_r(tindex, lat, a_s=0.25, b_s=0.5, n=None, nn=None):
    """Incoming solar radiation. Based on eq. 35 from FAO56.
    """
    ra = extraterrestrial_r(tindex, lat)
    if n is None:
        n = daylight_hours(tindex, lat)
    return (a_s + b_s * n / nn) * ra


def lai_calc(method=1, croph=None):
    if method == 1:
        return 0.24 * croph


def surface_r(lai=None, method=1, laieff=0, rl=100, co2=300, srs=0.0009,
              flai=1):
    if method == 1:
        return 70
    elif method == 2:
        lai_eff = calc_laieff(lai=lai, method=laieff)
        return rl / lai_eff
    elif method == 3:
        lai_eff = calc_laieff(lai=lai, method=laieff)
        fco2 = (1 + srs * (co2 - 300))
        return rl / lai_eff * fco2
    elif method == 4:
        lai_eff = calc_laieff(lai=lai, method=laieff)
        flai1 = lai_eff / lai_eff.max() * flai
        fco2 = (1 + srs * (co2 - 300))
        return rl / lai_eff * fco2 * flai1


def calc_laieff(lai=None, method=0):
    if method == 0:
        return 0.5 * lai
    if method == 1:
        return lai / (0.3 * lai + 1.2)
    if method == 2:
        laie = lai.copy()
        laie[(lai > 2) & (lai < 4)] = 2
        laie[lai > 4] = 0.5 * lai
        return laie
    if method == 3:
        laie = lai.copy()
        laie[lai > 4] = 4
        return laie * 0.5


def lambda_calc(temperature):
    """
    From ASCE (2001), eq. B.7
    """
    return 2.501 - 0.002361 * temperature


def calc_rhoa(pressure, ta, ea):
    tkv = (273.16 + ta) * (1 - 0.378 * ea / pressure) ** (-1)
    return 3.486 * pressure / tkv


def aero_r(wind, croph=None, zw=2, zh=2, method=1):
    """The aerodynamic resistance, applied for neutral stability conditions
     from ASCE (2001), eq. B.2

    Parameters
    ----------
    wind: Series
         wind speed at height z [m/s]
    zw: float
        height of wind measurements [m]
    zh: float
         height of humidity and or air temperature measurements [m]

    Returns
    -------
        pandas.Series containing the calculated aerodynamic resistance.

    """
    if method == 1:
        return 208 / wind
    elif method == 2:
        d = 0.667 * croph
        zom = 0.123 * croph
        zoh = 0.0123 * croph
        return (log((zw - d) / zom)) * \
               (log((zh - d) / zoh) /
                (0.41 ** 2) / wind)


def cloudiness_factor(rs, rso, ac=1.35, bc=-0.35):
    """
    Cloudiness factor f
    From FAO (1990), ANNEX V, eq. 57
    """
    return ac * rs / rso + bc


def rs_calc(tindex, lat, a_s=0.25, b_s=0.5):
    """
    Nncoming solar radiation rs
    From FAO (1990), ANNEX V, eq. 52
    """
    ra = extraterrestrial_r(tindex, lat)
    nn = 1
    return (a_s + b_s * nn) * ra


def ed_calc(tmax, tmin, rh):
    """
    Actual Vapour Pressure (ed).
    From FAO (1990), ANNEX V, eq. 11
    """
    eamax = e0_calc(tmax)
    eamin = e0_calc(tmin)
    return rh / (50. / eamin + 50. / eamax)


def calc_rns(solar=None, tindex=None, lat=None, alpha=0.23):
    """
    Net Shortwave Radiation Rns
    From FAO (1990), ANNEX V, eq. 51
    """
    if solar is not None:
        return (1 - alpha) * solar
    else:
        return (1 - alpha) * rs_calc(tindex, lat)


def calc_rnl(tmax, tmin, ea, cloudf, longa=0.34, longb=-0.139):
    """
    Net Longwave Radiation Rnl from FAO (1990), ANNEX V, eq. 56

    Parameters
    ----------
    tmax: Series
        maximum day temperature [°C]
    tmin: Series
        minimum day temperature [°C]

    Returns
    -------
        pandas.Series containing the calculated net outgoing radiation.

    """
    sigma = 0.00000000245 * ((tmax + 273.16) ** 4 + (tmin + 273.16) ** 4)
    emiss = longa + longb * round(sqrt(ea), 8)
    return sigma * cloudf * emiss
