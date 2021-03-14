from numpy import sqrt, log, exp, clip

from .utils import extraterrestrial_r, daylight_hours, \
    extraterrestrial_r_hour, day_of_year

# Specific heat of air [MJ kg-1 °C-1]
CP = 1.013 * 10 ** -3
# Stefan Boltzmann constant - hourly [MJm-2K-4h-1]
STEFAN_BOLTZMANN_HOUR = 2.042 * 10 ** -10
# Stefan Boltzmann constant - daily [MJm-2K-4d-1]
STEFAN_BOLTZMANN_DAY = 4.903 * 10 ** -9


def penman(wind, Rs=None, Rn=None, G=0, tmean=None, tmax=None, tmin=None,
           rhmax=None, rhmin=None, rh=None, pressure=None, elevation=None,
           lat=None, n=None, nn=None, Rso=None, aw=2.6, bw=0.536, a=1.35,
           b=-0.35):
    """Evaporation calculated according to [penman_1948]_.

    Parameters
    ----------
    wind: pandas.Series
        mean day wind speed [m/s]
    Rs: pandas.Series, optional
        incoming measured solar radiation [MJ m-2 d-1]
    Rn: pandas.Series, optional
        net radiation [MJ m-2 d-1]
    G: pandas.Series/int, optional
        soil heat flux [MJ m-2 d-1]
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
    pressure: float, optional
        atmospheric pressure [kPa]
    elevation: float, optional
        the site elevation [m]
    lat: float, optional
        the site latitude [rad]
    n: pandas.Series/float, optional
        actual duration of sunshine [hour]
    nn: pandas.Series/float, optional
        maximum possible duration of sunshine or daylight hours [hour]
    Rso: pandas.Series/float, optional
        clear-sky solar radiation [MJ m-2 day-1]
    aw: float, optional
        wind coefficient [-]
    bw: float, optional
        wind coefficient [-]
    a: float, optional
        empirical coefficient for Net Long-Wave radiation [-]
    b: float, optional
        empirical coefficient for Net Long-Wave radiation [-]

    Returns
    -------
    pandas.Series containing the calculated evaporation

    Examples
    --------
    >>> et_penman = penman(wind, Rn=Rn, tmean=tmean, rh=rh)

    .. math::
    -----
        E = \\frac{\\Delta (R_{n}-G)+ 6.43 \\gamma (a_w+b_w u_2)
        (e_{s}-e_{a})}{\\lambda (\\Delta +\\gamma)}

    References
    -----
    .. [penman_1948] Penman, H. L. (1948). Natural evaporation from open water,
       bare soil and grass. Proceedings of the Royal Society of London. Series
       A. Mathematical and Physical Sciences, 193, 120-145.
    .. [valiantzas_2006] Valiantzas, J. D. (2006). Simplified versions for the
       Penman evaporation equation using routine weather data. Journal of
       Hydrology, 331(3-4), 690-702.

    """
    if tmean is None:
        tmean = (tmax + tmin) / 2
    if pressure is None:
        pressure = calc_press(elevation)
    gamma = calc_psy(pressure)
    dlt = calc_vpc(tmean)
    lambd = calc_lambda(tmean)

    ea = calc_ea(tmean=tmean, tmax=tmax, tmin=tmin, rhmax=rhmax, rhmin=rhmin,
                 rh=rh)
    es = calc_es(tmean=tmean, tmax=tmax, tmin=tmin)

    if Rn is None:
        Rns = calc_rad_short(Rs=Rs, n=n, nn=nn)  # [MJ/m2/d]
        Rnl = calc_rad_long(Rs=Rs, tmean=tmean, tmax=tmax, tmin=tmin,
                            rhmax=rhmax, rhmin=rhmin, rh=rh,
                            elevation=elevation, lat=lat, Rso=Rso, a=a, b=b,
                            ea=ea)  # [MJ/m2/d]
        Rn = Rns - Rnl

    fu = aw * (1 + bw * wind)

    den = lambd * (dlt + gamma)
    num1 = dlt * (Rn - G) / den
    num2 = gamma * (es - ea) * fu / den
    return num1 + num2


def pm_fao56(wind, Rs=None, Rn=None, G=0, tmean=None, tmax=None, tmin=None,
             rhmax=None, rhmin=None, rh=None, pressure=None, elevation=None,
             lat=None, n=None, nn=None, Rso=None, a=1.35, b=-0.35):
    """Evaporation calculated according to [allen_1998]_..

    Parameters
    ----------
    wind: pandas.Series
        mean day wind speed [m/s]
    Rs: pandas.Series, optional
        incoming measured solar radiation [MJ m-2 d-1]
    Rn: pandas.Series, optional
        net radiation [MJ m-2 d-1]
    G: pandas.Series/int, optional
        soil heat flux [MJ m-2 d-1]
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
    pressure: float, optional
        atmospheric pressure [kPa]
    elevation: float, optional
        the site elevation [m]
    lat: float, optional
        the site latitude [rad]
    n: pandas.Series/float, optional
        actual duration of sunshine [hour]
    nn: pandas.Series/float, optional
        maximum possible duration of sunshine or daylight hours [hour]
    Rso: pandas.Series/float, optional
        clear-sky solar radiation [MJ m-2 day-1]
    a: float, optional
        empirical coefficient for Net Long-Wave radiation [-]
    b: float, optional
        empirical coefficient for Net Long-Wave radiation [-]

    Returns
    -------
        pandas.Series containing the calculated evaporation

    Examples
    --------
    >>> et_fao56 = pm_fao56(wind, Rn=Rn, tmean=tmean, rh=rh)

    .. math::
    -----
        $E = \\frac{0.408 \\Delta (R_{n}-G)+ \\gamma \\frac{900}{T+273}
        (e_{s}-e_{a})}{\\Delta +\\gamma(1+0.34u_2)}$

    References
    -----
    .. [allen_1998] Allen, R. G., Pereira, L. S., Raes, D., & Smith, M. (1998).
       Crop evapotranspiration-Guidelines for computing crop water
       requirements-FAO Irrigation and drainage paper 56. Fao, Rome, 300.
       (http://www.fao.org/3/x0490e/x0490e06.htm#TopOfPage)

    """
    if tmean is None:
        tmean = (tmax + tmin) / 2
    if pressure is None:
        pressure = calc_press(elevation)
    gamma = calc_psy(pressure)
    dlt = calc_vpc(tmean)

    gamma1 = (gamma * (1 + 0.34 * wind))

    ea = calc_ea(tmean=tmean, tmax=tmax, tmin=tmin, rhmax=rhmax, rhmin=rhmin,
                 rh=rh)
    es = calc_es(tmean=tmean, tmax=tmax, tmin=tmin)

    if Rn is None:
        Rns = calc_rad_short(Rs=Rs, n=n, nn=nn)  # [MJ/m2/d]
        Rnl = calc_rad_long(Rs=Rs, tmean=tmean, tmax=tmax, tmin=tmin,
                            rhmax=rhmax, rhmin=rhmin, rh=rh,
                            elevation=elevation, lat=lat, Rso=Rso, a=a, b=b,
                            ea=ea)  # [MJ/m2/d]
        Rn = Rns - Rnl

    den = dlt + gamma1
    num1 = (0.408 * dlt * (Rn - G)) / den
    num2 = (gamma * (es - ea) * 900 * wind / (tmean + 273)) / den
    return num1 + num2


def pm(wind, Rs=None, Rn=None, G=0, tmean=None, tmax=None, tmin=None,
       rhmax=None, rhmin=None, rh=None, pressure=None, elevation=None,
       lat=None, n=None, nn=None, Rso=None, a=1.35, b=-0.35, lai=None,
       croph=None, rl=100, rs=70, ra_method=1, a_sh=1, a_s=1, lai_eff=1,
       srs=0.0009, co2=300):
    """Evaporation calculated according to [monteith_1965]_.

    Parameters
    ----------
    wind: pandas.Series
        mean day wind speed [m/s]
    Rs: pandas.Series, optional
        incoming measured solar radiation [MJ m-2 d-1]
    Rn: pandas.Series, optional
        net radiation [MJ m-2 d-1]
    G: pandas.Series/int, optional
        soil heat flux [MJ m-2 d-1]
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
    pressure: float, optional
        atmospheric pressure [kPa]
    elevation: float, optional
        the site elevation [m]
    lat: float, optional
        the site latitude [rad]
    n: pandas.Series/float, optional
        actual duration of sunshine [hour]
    nn: pandas.Series/float, optional
        maximum possible duration of sunshine or daylight hours [hour]
    Rso: pandas.Series/float, optional
        clear-sky solar radiation [MJ m-2 day-1]
    a: float, optional
        empirical coefficient for Net Long-Wave radiation [-]
    b: float, optional
        empirical coefficient for Net Long-Wave radiation [-]
    lai: pandas.Series/float, optional
        measured leaf area index [-]
    croph: pandas.series/float, optional
        crop height [m]
    rl: pandas.series/float, optional
        bulk stomatal resistance [s m-1]
    rs: pandas.series/float, optional
        bulk surface resistance [s m-1]
    ra_method: float, optional
        1 => ra = 208/wind
        2 => ra is calculated based on equation 36 in FAO (1990), ANNEX V.
    a_s: float, optional
        Fraction of one-sided leaf area covered by stomata (1 if stomata are 1
        on one side only, 2 if they are on both sides)
    a_sh: float, optional
        Fraction of projected area exchanging sensible heat with the air (2)
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
        pandas.Series containing the calculated evaporation

    .. math::
    -----
        $E = \\frac{\\Delta (R_{n}-G)+ \\rho_a c_p K_{min} \\frac{e_{s}-e_{a}}
        {r_a}}{\\lambda(\\Delta +\\gamma(1+\\frac{r_s}{r_a}))}$

    References
    -----
    .. [monteith_1965] Monteith, J. L. (1965). Evaporation and environment.
       In Symposia of the society for experimental biology (Vol. 19, pp.
       205-234). Cambridge University Press (CUP) Cambridge.
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
    if tmean is None:
        tmean = (tmax + tmin) / 2
    if pressure is None:
        pressure = calc_press(elevation)
    gamma = calc_psy(pressure)
    dlt = calc_vpc(tmean)
    lambd = calc_lambda(tmean)

    ea = calc_ea(tmean=tmean, tmax=tmax, tmin=tmin, rhmax=rhmax, rhmin=rhmin,
                 rh=rh)
    es = calc_es(tmean=tmean, tmax=tmax, tmin=tmin)

    res_a = calc_res_aero(wind, ra_method=ra_method, croph=croph)
    res_s = calc_res_surf(lai=lai, rs=rs, rl=rl, lai_eff=lai_eff, srs=srs,
                          co2=co2)
    gamma1 = gamma * a_sh / a_s * (1 + res_s / res_a)

    if Rn is None:
        Rns = calc_rad_short(Rs=Rs, n=n, nn=nn)  # [MJ/m2/d]
        Rnl = calc_rad_long(Rs=Rs, tmean=tmean, tmax=tmax, tmin=tmin,
                            rhmax=rhmax, rhmin=rhmin, rh=rh,
                            elevation=elevation, lat=lat, Rso=Rso, a=a, b=b,
                            ea=ea)  # [MJ/m2/d]
        Rn = Rns - Rnl

    kmin = 86400  # unit conversion s d-1
    rho_a = calc_rho(pressure, tmean, ea)

    den = lambd * (dlt + gamma1)
    num1 = dlt * (Rn - G) / den
    num2 = rho_a * CP * kmin * (es - ea) * a_sh / res_a / den
    return num1 + num2


def kimberly_penman(wind, Rs=None, Rn=None, G=0, tmean=None, tmax=None,
                    tmin=None, rhmax=None, rhmin=None, rh=None, pressure=None,
                    elevation=None, lat=None, n=None, nn=None, Rso=None,
                    a=1.35, b=-0.35):
    """Evaporation calculated according to [wright_1982]_.

    Parameters
    ----------
    wind: pandas.Series
        mean day wind speed [m/s]
    Rs: pandas.Series, optional
        incoming measured solar radiation [MJ m-2 d-1]
    Rn: pandas.Series, optional
        net radiation [MJ m-2 d-1]
    G: pandas.Series/int, optional
        soil heat flux [MJ m-2 d-1]
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
    pressure: float, optional
        atmospheric pressure [kPa]
    elevation: float, optional
        the site elevation [m]
    lat: float, optional
        the site latitude [rad]
    n: pandas.Series/float, optional
        actual duration of sunshine [hour]
    nn: pandas.Series/float, optional
        maximum possible duration of sunshine or daylight hours [hour]
    Rso: pandas.Series/float, optional
        clear-sky solar radiation [MJ m-2 day-1]
    a: float, optional
        empirical coefficient for Net Long-Wave radiation [-]
    b: float, optional
        empirical coefficient for Net Long-Wave radiation [-]

    Returns
    -------
        pandas.Series containing the calculated evaporation

    .. math::
    -----
        $E = \\frac{\\Delta (R_{n}-G)+ \\rho_a c_p K_{min} \\frac{e_{s}-e_{a}}
        {r_a}}{\\lambda(\\Delta +\\gamma(1+\\frac{r_s}{r_a}))}$

    References
    -----
    .. [wright_1982] Wright, J. L. (1982). New evapotranspiration crop
       coefficients. Proceedings of the American Society of Civil Engineers,
       Journal of the Irrigation and Drainage Division, 108(IR2), 57-74.
    .. [oudin_2005] Oudin, L., Hervieu, F., Michel, C., Perrin, C.,
       Andréassian, V., Anctil, F., & Loumagne, C. (2005). Which potential
       evapotranspiration input for a lumped rainfall–runoff model?:
       Part 2—Towards a simple and efficient potential evapotranspiration model
       for rainfall–runoff modelling. Journal of hydrology, 303(1-4), 290-306.

    """
    if tmean is None:
        tmean = (tmax + tmin) / 2
    if pressure is None:
        pressure = calc_press(elevation)
    gamma = calc_psy(pressure)
    dlt = calc_vpc(tmean)
    lambd = calc_lambda(tmean)

    ea = calc_ea(tmean=tmean, tmax=tmax, tmin=tmin, rhmax=rhmax, rhmin=rhmin,
                 rh=rh)
    es = calc_es(tmean=tmean, tmax=tmax, tmin=tmin)

    if Rn is None:
        Rns = calc_rad_short(Rs=Rs, n=n, nn=nn)  # [MJ/m2/d]
        Rnl = calc_rad_long(Rs=Rs, tmean=tmean, tmax=tmax, tmin=tmin,
                            rhmax=rhmax, rhmin=rhmin, rh=rh,
                            elevation=elevation, lat=lat, Rso=Rso, a=a, b=b,
                            ea=ea)  # [MJ/m2/d]
        Rn = Rns - Rnl

    j = day_of_year(tmean.index)
    w = wind * (0.4 + 0.14 * exp(-((j - 173)/58)**2) + (0.605 + 0.345 *
                                                        exp((j - 243)/80)**2))

    den = lambd * (dlt + gamma)
    num1 = dlt * (Rn - G) / den
    num2 = gamma * (es - ea) * w / den
    return num1 + num2


def fao_24(wind, Rs=None, Rn=None, tmean=None, tmax=None, tmin=None, rh=None,
           pressure=None, elevation=None, albedo=0.23):
    """Evaporation calculated according to [priestley_and_taylor_1965]_.

    Parameters
    ----------
    wind: pandas.Series
        mean day wind speed [m/s]
    Rs: pandas.Series, optional
        incoming measured solar radiation [MJ m-2 d-1]
    tmean: pandas.Series, optional
        average day temperature [°C]
    tmax: pandas.Series, optional
        maximum day temperature [°C]
    tmin: pandas.Series, optional
        minimum day temperature [°C]
    rh: pandas.Series, optional
        mean daily relative humidity [%]
    pressure: float, optional
        atmospheric pressure [kPa]
    elevation: float, optional
        the site elevation [m]
    albedo: float, optional
        surface albedo [-]

    Returns
    -------
        pandas.Series containing the calculated evapotranspiration

    Examples
    --------
    >>> pt = priestley_taylor(wind, Rn=Rn, tmean=tmean, rh=rh)

    .. math::
    -----
        $E = \\frac{\\Delta (R_{n}-G)+ \\rho_a c_p K_{min} \\frac{e_{s}-e_{a}}
        {r_a}}{\\lambda(\\Delta +\\gamma(1+\\frac{r_s}{r_a}))}$

    References
    -----
    .. [priestley_and_taylor_1965] Priestley, C. H. B., & TAYLOR, R. J. (1972).
       On the assessment of surface heat flux and evaporation using large-scale
       parameters. Monthly weather review, 100(2), 81-92.

    """
    if tmean is None:
        tmean = (tmax + tmin) / 2
    if pressure is None:
        pressure = calc_press(elevation)
    gamma = calc_psy(pressure)
    dlt = calc_vpc(tmean)
    lambd = calc_lambda(tmean)

    w = 1.066 - 0.13 * rh/100 + 0.045 * wind - 0.023 * rh/100 * wind - 3.15 * \
        (rh/100)**2 - 0.0011 * wind

    return -0.3 + dlt * (dlt + gamma) * Rs * (1 - albedo) * w / lambd


def thom_oliver(wind, Rs=None, Rn=None, G=0, tmean=None, tmax=None, tmin=None,
                rhmax=None, rhmin=None, rh=None, pressure=None, elevation=None,
                lat=None, n=None, nn=None, Rso=None, a=1.35, b=-0.35, lai=None,
                croph=None, rl=100, rs=70, ra_method=1, lai_eff=1, srs=0.0009,
                co2=300):
    """Evaporation calculated according to [thom_1977]_.

    Parameters
    ----------
    wind: pandas.Series
        mean day wind speed [m/s]
    Rs: pandas.Series, optional
        incoming measured solar radiation [MJ m-2 d-1]
    Rn: pandas.Series, optional
        net radiation [MJ m-2 d-1]
    G: pandas.Series/int, optional
        soil heat flux [MJ m-2 d-1]
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
    pressure: float, optional
        atmospheric pressure [kPa]
    elevation: float, optional
        the site elevation [m]
    lat: float, optional
        the site latitude [rad]
    n: pandas.Series/float, optional
        actual duration of sunshine [hour]
    nn: pandas.Series/float, optional
        maximum possible duration of sunshine or daylight hours [hour]
    Rso: pandas.Series/float, optional
        clear-sky solar radiation [MJ m-2 day-1]
    a: float, optional
        empirical coefficient for Net Long-Wave radiation [-]
    b: float, optional
        empirical coefficient for Net Long-Wave radiation [-]
    lai: pandas.Series/float, optional
        measured leaf area index [-]
    croph: pandas.series/float, optional
        crop height [m]
    rl: pandas.series/float, optional
        bulk stomatal resistance [s m-1]
    rs: pandas.series/float, optional
        bulk surface resistance [s m-1]
    ra_method: float, optional
        1 => ra = 208/wind
        2 => ra is calculated based on equation 36 in FAO (1990), ANNEX V.
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
        pandas.Series containing the calculated evaporation

    .. math::
    -----
        $E = \\frac{\\Delta (R_{n}-G)+ \\rho_a c_p K_{min} \\frac{e_{s}-e_{a}}
        {r_a}}{\\lambda(\\Delta +\\gamma(1+\\frac{r_s}{r_a}))}$

    References
    -----
    .. [thom_1977] Thom, A. S., & Oliver, H. R. (1977). On Penman's equation
       for estimating regional evaporation. Quarterly Journal of the Royal
       Meteorological Society, 103(436), 345-357.
    .. [oudin_2005] Oudin, L., Hervieu, F., Michel, C., Perrin, C.,
       Andréassian, V., Anctil, F., & Loumagne, C. (2005). Which potential
       evapotranspiration input for a lumped rainfall–runoff model?:
       Part 2—Towards a simple and efficient potential evapotranspiration model
       for rainfall–runoff modelling. Journal of hydrology, 303(1-4), 290-306.

    """
    if tmean is None:
        tmean = (tmax + tmin) / 2
    if pressure is None:
        pressure = calc_press(elevation)
    gamma = calc_psy(pressure)
    dlt = calc_vpc(tmean)
    lambd = calc_lambda(tmean)

    ea = calc_ea(tmean=tmean, tmax=tmax, tmin=tmin, rhmax=rhmax, rhmin=rhmin,
                 rh=rh)
    es = calc_es(tmean=tmean, tmax=tmax, tmin=tmin)

    res_a = calc_res_aero(wind, ra_method=ra_method, croph=croph)
    res_s = calc_res_surf(lai=lai, rs=rs, rl=rl, lai_eff=lai_eff, srs=srs,
                          co2=co2)
    gamma1 = gamma * (1 + res_s / res_a)

    if Rn is None:
        Rns = calc_rad_short(Rs=Rs, n=n, nn=nn)  # [MJ/m2/d]
        Rnl = calc_rad_long(Rs=Rs, tmean=tmean, tmax=tmax, tmin=tmin,
                            rhmax=rhmax, rhmin=rhmin, rh=rh,
                            elevation=elevation, lat=lat, Rso=Rso, a=a, b=b,
                            ea=ea)  # [MJ/m2/d]
        Rn = Rns - Rnl

    w = 2.6 * (1 + 0.536 * wind)

    den = lambd * (dlt + gamma1)
    num1 = dlt * (Rn - G) / den
    num2 = 2.5 * gamma * (es - ea) * w / den
    return num1 + num2


def priestley_taylor(wind, Rs=None, Rn=None, G=0, tmean=None, tmax=None,
                     tmin=None, rhmax=None, rhmin=None, rh=None, pressure=None,
                     elevation=None, lat=None, n=None, nn=None, Rso=None,
                     a=1.35, b=-0.35, alpha=1.26):
    """Evaporation calculated according to [priestley_and_taylor_1965]_.

    Parameters
    ----------
    wind: pandas.Series
        mean day wind speed [m/s]
    Rs: pandas.Series, optional
        incoming measured solar radiation [MJ m-2 d-1]
    Rn: pandas.Series, optional
        net radiation [MJ m-2 d-1]
    G: pandas.Series/int, optional
        soil heat flux [MJ m-2 d-1]
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
    pressure: float, optional
        atmospheric pressure [kPa]
    elevation: float, optional
        the site elevation [m]
    lat: float, optional
        the site latitude [rad]
    n: pandas.Series/float, optional
        actual duration of sunshine [hour]
    nn: pandas.Series/float, optional
        maximum possible duration of sunshine or daylight hours [hour]
    Rso: pandas.Series/float, optional
        clear-sky solar radiation [MJ m-2 day-1]
    a: float, optional
        empirical coefficient for Net Long-Wave radiation [-]
    b: float, optional
        empirical coefficient for Net Long-Wave radiation [-]
    alpha: float, optional
        calibration coeffiecient [-]

    Returns
    -------
        pandas.Series containing the calculated evapotranspiration

    Examples
    --------
    >>> pt = priestley_taylor(wind, Rn=Rn, tmean=tmean, rh=rh)

    .. math::
    -----
        $E = \\frac{\\Delta (R_{n}-G)+ \\rho_a c_p K_{min} \\frac{e_{s}-e_{a}}
        {r_a}}{\\lambda(\\Delta +\\gamma(1+\\frac{r_s}{r_a}))}$

    References
    -----
    .. [priestley_and_taylor_1965] Priestley, C. H. B., & TAYLOR, R. J. (1972).
       On the assessment of surface heat flux and evaporation using large-scale
       parameters. Monthly weather review, 100(2), 81-92.

    """
    if tmean is None:
        tmean = (tmax + tmin) / 2
    if pressure is None:
        pressure = calc_press(elevation)
    gamma = calc_psy(pressure)
    dlt = calc_vpc(tmean)
    lambd = calc_lambda(tmean)

    if Rn is None:
        Rns = calc_rad_short(Rs=Rs, n=n, nn=nn)  # [MJ/m2/d]
        Rnl = calc_rad_long(Rs=Rs, tmean=tmean, tmax=tmax, tmin=tmin,
                            rhmax=rhmax, rhmin=rhmin, rh=rh,
                            elevation=elevation, lat=lat, Rso=Rso, a=a, b=b)
        Rn = Rns - Rnl

    return (alpha * dlt * (Rn-G)) / (lambd * (dlt + gamma))


def makkink(Rs, tmean=None, tmax=None, tmin=None, pressure=None,
            elevation=None):
    """"Evaporation calculated according to [makkink_1965]_.

    Parameters
    ----------
    Rs: pandas.Series, optional
        incoming measured solar radiation [MJ m-2 d-1]
    tmean: pandas.Series, optional
        average day temperature [°C]
    tmax: pandas.Series, optional
        maximum day temperature [°C]
    tmin: pandas.Series, optional
        minimum day temperature [°C]
    pressure: float, optional
        atmospheric pressure [kPa]
    elevation: float, optional
        the site elevation [m]

    Returns
    -------
        pandas.Series containing the calculated evaporation

    Examples
    --------
    >>> mak = makkink(Rs=Rs, tmean=tmean)

    .. math::
    -----
        $E = \\frac{\\Delta (R_{n}-G)+ \\rho_a c_p K_{min} \\frac{e_{s}-e_{a}}
        {r_a}}{\\lambda(\\Delta +\\gamma(1+\\frac{r_s}{r_a}))}$

    References
    -----
    .. [priestley_and_taylor_1965] Priestley, C. H. B., & TAYLOR, R. J. (1972).
       On the assessment of surface heat flux and evaporation using large-scale
       parameters. Monthly weather review, 100(2), 81-92.

    """
    if tmean is None:
        tmean = (tmax + tmin) / 2
    if pressure is None:
        pressure = calc_press(elevation)
    gamma = calc_psy(pressure)
    dlt = calc_vpc(tmean)
    lambd = calc_lambda(tmean)

    return 0.65 * dlt / (dlt + gamma) * Rs / lambd


##% Utility functions (TODO: Make private?)


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
    -----
    .. [allen_1998] Allen, R. G., Pereira, L. S., Raes, D., & Smith, M. (1998).
       Crop evapotranspiration-Guidelines for computing crop water
       requirements-FAO Irrigation and drainage paper 56. Fao, Rome, 300.
       (http://www.fao.org/3/x0490e/x0490e06.htm#TopOfPage)

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

    References
    -----
    .. [allen_1998] Allen, R. G., Pereira, L. S., Raes, D., & Smith, M. (1998).
       Crop evapotranspiration-Guidelines for computing crop water
       requirements-FAO Irrigation and drainage paper 56. Fao, Rome, 300.
       (http://www.fao.org/3/x0490e/x0490e06.htm#TopOfPage)

    """
    es = calc_e0(tmean)
    return 4098 * es / (tmean + 237.3) ** 2


def calc_lambda(tmean):
    """ Latent Heat of Vaporization [MJ kg-1].

    Parameters
    ----------
    tmean: pandas.Series, optional
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

    References
    -----
    .. [allen_1998] Allen, R. G., Pereira, L. S., Raes, D., & Smith, M. (1998).
       Crop evapotranspiration-Guidelines for computing crop water
       requirements-FAO Irrigation and drainage paper 56. Fao, Rome, 300.
       (http://www.fao.org/3/x0490e/x0490e06.htm#TopOfPage)

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

    References
    -----
    .. [allen_1998] Allen, R. G., Pereira, L. S., Raes, D., & Smith, M. (1998).
       Crop evapotranspiration-Guidelines for computing crop water
       requirements-FAO Irrigation and drainage paper 56. Fao, Rome, 300.
       (http://www.fao.org/3/x0490e/x0490e06.htm#TopOfPage)

    """
    return 101.3 * ((293 - 0.0065 * elevation) / 293) ** 5.26


def calc_rho(pressure, tmean, ea):
    """atmospheric air density calculated according to [allen_1998]_..

    Parameters
    ----------
    pressure: pandas.Series
        atmospheric pressure [kPa]
    tmean: pandas.Series, optional
        average day temperature [°C]
    ea: pandas.Series, optional
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

    .. math::
    -----
        rho = 3.486 \\frac{P}{T_{KV}}

    References
    -----
    .. [allen_1998] Allen, R. G., Pereira, L. S., Raes, D., & Smith, M. (1998).
       Crop evapotranspiration-Guidelines for computing crop water
       requirements-FAO Irrigation and drainage paper 56. Fao, Rome, 300.
       (http://www.fao.org/3/x0490e/x0490e06.htm#TopOfPage)

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
        pandas.Series containing the calculated saturation vapor pressure at
        the air temperature tmean [kPa].

    Examples
    --------
    >>> e0 = calc_e0(tmean)

    Notes
    -----
        Based on equation 11 in [allen_1998]_.

    References
    -----
    .. [allen_1998] Allen, R. G., Pereira, L. S., Raes, D., & Smith, M. (1998).
       Crop evapotranspiration-Guidelines for computing crop water
       requirements-FAO Irrigation and drainage paper 56. Fao, Rome, 300.
       (http://www.fao.org/3/x0490e/x0490e06.htm#TopOfPage)

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
        pandas.Series containing the calculated saturation vapor pressure
        [kPa].

    Examples
    --------
    >>> es = calc_es(tmean)

    Notes
    -----
        Based on equation 11, 12 in [allen_1998]_.

    References
    -----
    .. [allen_1998] Allen, R. G., Pereira, L. S., Raes, D., & Smith, M. (1998).
       Crop evapotranspiration-Guidelines for computing crop water
       requirements-FAO Irrigation and drainage paper 56. Fao, Rome, 300.
       (http://www.fao.org/3/x0490e/x0490e06.htm#TopOfPage)

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
        pandas.Series containing the calculated actual vapor pressure
        [kPa].

    Examples
    --------
    >>> ea = calc_ea(tmean, rh)

    Notes
    -----
        Based on equation 17, 19 in [allen_1998]_.

    References
    -----
    .. [allen_1998] Allen, R. G., Pereira, L. S., Raes, D., & Smith, M. (1998).
       Crop evapotranspiration-Guidelines for computing crop water
       requirements-FAO Irrigation and drainage paper 56. Fao, Rome, 300.
       (http://www.fao.org/3/x0490e/x0490e06.htm#TopOfPage)

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


def calc_rad_long(Rs, tmean=None, tmax=None, tmin=None, rhmax=None,
                  rhmin=None, rh=None, elevation=None, lat=None, Rso=None,
                  a=1.35, b=-0.35, ea=None, freq="D"):
    """Net longwave radiation [MJ m-2 d-1].

    Parameters
    ----------
    Rs: pandas.Series, optional
        incoming measured solar radiation [MJ m-2 d-1]
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
    Rso: pandas.Series/float, optional
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
        if Rso is None:
            Ra = extraterrestrial_r_hour(tindex=Rs.index, lat=lat)
            Rso = calc_rso(Ra=Ra, elevation=elevation)
        solar_rat = clip(Rs / Rso, 0.3, 1)
        tmp1 = STEFAN_BOLTZMANN_HOUR * (tmean + 273.2) ** 4
    else:
        if Rso is None:
            Ra = extraterrestrial_r(tindex=Rs.index, lat=lat)
            Rso = calc_rso(Ra=Ra, elevation=elevation)
        solar_rat = clip(Rs / Rso, 0.3, 1)
        if tmax is not None:
            tmp1 = STEFAN_BOLTZMANN_DAY * ((tmax + 273.2) ** 4 +
                                           (tmin + 273.2) ** 4) / 2
        else:
            tmp1 = STEFAN_BOLTZMANN_DAY * (tmean + 273.2) ** 4

    tmp2 = 0.34 - 0.139 * sqrt(ea)  # OK
    tmp2 = clip(tmp2, 0.05, 1)
    tmp3 = a * solar_rat + b  # OK
    return tmp1 * tmp2 * tmp3


def calc_rad_short(Rs=None, tindex=None, lat=None, alpha=0.23, n=None,
                   nn=None):
    """Net shortwave radiation [MJ m-2 d-1].

    Parameters
    ----------
    Rs: pandas.Series, optional
        incoming measured solar radiation [MJ m-2 d-1]
    tindex: pandas.Series.index
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
    if Rs is not None:
        return (1 - alpha) * Rs
    else:
        return (1 - alpha) * calc_rad_sol_in(tindex, lat, n=n, nn=nn)


def calc_rad_sol_in(tindex, lat, as1=0.25, bs1=0.5, n=None, nn=None):
    """Incoming solar radiation [MJ m-2 d-1].

    Parameters
    ----------
    tindex: pandas.Series.index
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
    Ra = extraterrestrial_r(tindex, lat)
    if n is None:
        n = daylight_hours(tindex, lat)
    return (as1 + bs1 * n / nn) * Ra


def calc_rso(Ra, elevation):
    """Clear-sky solar radiation [MJ m-2 day-1].

    Parameters
    ----------
    Ra: pandas.Series, optional
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
    return (0.75 + (2 * 10 ** -5) * elevation) * Ra


def calc_res_surf(lai=None, rs=70, rl=100, lai_eff=0, srs=None, co2=None):
    """Surface resistance [s m-1].

    Parameters
    ----------
    lai: pandas.Series/float, optional
        measured leaf area index [-]
    rs: pandas.series/float, optional
        surface resistance [s m-1]
    rl: float, optional
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
    if lai:
        fco2 = (1 + srs * (co2 - 300))
        return fco2 * rl / calc_laieff(lai=lai, lai_eff=lai_eff)
    else:
        return rs


def calc_laieff(lai=None, lai_eff=0):
    """Effective leaf area index [-].

    Parameters
    ----------
    lai: pandas.Series/float, optional
        measured leaf area index [-]
    lai_eff: float, optional
        1 => LAI_eff = 0.5 * LAI
        2 => LAI_eff = lai / (0.3 * lai + 1.2)
        3 => LAI_eff = 0.5 * LAI; (LAI>4=4)
        4 => see [zhang_2008]_.

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
    .. [schymanski_2016] Schymanski, S. J., & Or, D. (2017). Leaf-scale
       experiments reveal an important omission in the Penman–Monteith
       equation. Hydrology and Earth System Sciences, 21(2), 685-706.
    .. [yang_2019] Yang, Y., Roderick, M. L., Zhang, S., McVicar, T. R., &
       Donohue, R. J. (2019). Hydrologic implications of vegetation response to
       elevated CO 2 in climate projections. Nature Climate Change, 9, 44-48.
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
