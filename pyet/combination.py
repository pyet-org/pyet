"""The combination module contains functions of combination PE methods

"""
import pandas

from .meteo_utils import *
from .rad_utils import *
from .temperature import *
from .radiation import *
from .utils import *

# Specific heat of air [MJ kg-1 °C-1]
CP = 1.013 * 10 ** -3
# Stefan Boltzmann constant - hourly [MJm-2K-4h-1]
STEFAN_BOLTZMANN_HOUR = 2.042 * 10 ** -10
# Stefan Boltzmann constant - daily [MJm-2K-4d-1]
STEFAN_BOLTZMANN_DAY = 4.903 * 10 ** -9


def penman(tmean, wind, rs=None, rn=None, g=0, tmax=None, tmin=None,
           rhmax=None, rhmin=None, rh=None, pressure=None, elevation=None,
           lat=None, n=None, nn=None, rso=None, aw=2.6, bw=0.536, a=1.35,
           b=-0.35, ea=None, albedo=0.23, clip_zero=True):
    """Potential evaporation calculated according to [penman_1948]_.

    Parameters
    ----------
    tmean: pandas.Series/xarray.DataArray
        average day temperature [°C]
    wind: float/pandas.Series/xarray.DataArray
        mean day wind speed [m/s]
    rs: float/pandas.Series/xarray.DataArray, optional
        incoming solar radiation [MJ m-2 d-1]
    rn: float/pandas.Series/xarray.DataArray, optional
        net radiation [MJ m-2 d-1]
    g: float/pandas.Series/xarray.DataArray, optional
        soil heat flux [MJ m-2 d-1]
    tmax: float/pandas.Series/xarray.DataArray, optional
        maximum day temperature [°C]
    tmin: float/pandas.Series/xarray.DataArray, optional
        minimum day temperature [°C]
    rhmax: float/pandas.Series/xarray.DataArray, optional
        maximum daily relative humidity [%]
    rhmin: float/pandas.Series/xarray.DataArray, optional
        mainimum daily relative humidity [%]
    rh: float/pandas.Series/xarray.DataArray, optional
        mean daily relative humidity [%]
    pressure: float/pandas.Series/xarray.DataArray, optional
        atmospheric pressure [kPa]
    elevation: float/xarray.DataArray, optional
        the site elevation [m]
    lat: float/xarray.DataArray, optional
        the site latitude [rad]
    n: float/pandas.Series/xarray.DataArray, optional
        actual duration of sunshine [hour]
    nn: float/pandas.Series/xarray.DataArray, optional
        maximum possible duration of sunshine or daylight hours [hour]
    rso: float/pandas.Series/xarray.DataArray, optional
        clear-sky solar radiation [MJ m-2 day-1]
    aw: float, optional
        wind coefficient [-]
    bw: float, optional
        wind coefficient [-]
    a: float, optional
        empirical coefficient for Net Long-Wave radiation [-]
    b: float, optional
        empirical coefficient for Net Long-Wave radiation [-]
    ea: float/pandas.Series/xarray.DataArray, optional
        actual vapor pressure [kPa]
    albedo: float, optional
        surface albedo [-]
    clip_zero: bool, optional
        if True, replace all negative values with 0.

    Returns
    -------
    pandas.Series/xarray.DataArray containing the calculated
            potential evaporation [mm d-1].

    Examples
    --------
    >>> et_penman = penman(tmean, wind, rn=rn, rh=rh)

    Notes
    -----
    Following [penman_1948]_ and [valiantzas_2006]_.

    .. math:: PE = \\frac{\\Delta (R_n-G) + \\gamma 2.6 (1 + 0.536 u_2)
        (e_s-e_a)}{\\lambda (\\Delta +\\gamma)}

    References
    ----------
    .. [penman_1948] Penman, H. L. (1948). Natural evaporation from open water,
       bare soil and grass. Proceedings of the Royal Society of London. Series
       A. Mathematical and Physical Sciences, 193, 120-145.
    .. [valiantzas_2006] Valiantzas, J. D. (2006). Simplified versions for the
       Penman evaporation equation using routine weather data. Journal of
       Hydrology, 331(3-4), 690-702.

    """
    pressure = calc_press(elevation, pressure)
    gamma = calc_psy(pressure)
    dlt = calc_vpc(tmean)
    lambd = calc_lambda(tmean)

    ea = calc_ea(tmean=tmean, tmax=tmax, tmin=tmin, rhmax=rhmax,
                 rhmin=rhmin, rh=rh, ea=ea)
    es = calc_es(tmean=tmean, tmax=tmax, tmin=tmin)

    if rn is None:
        rn = get_rn(tmean, rs, check_lat(lat), n, nn, tmax, tmin, rhmax, rhmin,
                    rh, elevation, rso, a, b, ea, albedo)

    fu = aw * (1 + bw * wind)

    den = lambd * (dlt + gamma)
    num1 = dlt * (rn - g) / den
    num2 = gamma * (es - ea) * fu / den
    pe = num1 + num2
    pe = clip_zeros(pe, clip_zero)
    return pe.rename("Penman")


def pm_asce(tmean, wind, rs=None, rn=None, g=0, tmax=None, tmin=None,
            rhmax=None, rhmin=None, rh=None, pressure=None, elevation=None,
            lat=None, n=None, nn=None, rso=None, a=1.35, b=-0.35, cn=900,
            cd=0.34, ea=None, albedo=0.23, kab=None, as1=0.25, bs1=0.5,
            clip_zero=True, etype="os"):
    """Potential evaporation calculated according to [monteith_1965]_.

    Parameters
    ----------
    tmean: pandas.Series/xarray.DataArray
        average day temperature [°C]
    wind: float/pandas.Series/xarray.DataArray
        mean day wind speed [m/s]
    rs: float/pandas.Series/xarray.DataArray, optional
        incoming solar radiation [MJ m-2 d-1]
    rn: float/pandas.Series/xarray.DataArray, optional
        net radiation [MJ m-2 d-1]
    g: float/pandas.Series/xarray.DataArray, optional
        soil heat flux [MJ m-2 d-1]
    tmax: float/pandas.Series/xarray.DataArray, optional
        maximum day temperature [°C]
    tmin: float/pandas.Series/xarray.DataArray, optional
        minimum day temperature [°C]
    rhmax: float/pandas.Series/xarray.DataArray, optional
        maximum daily relative humidity [%]
    rhmin: float/pandas.Series/xarray.DataArray, optional
        mainimum daily relative humidity [%]
    rh: float/pandas.Series/xarray.DataArray, optional
        mean daily relative humidity [%]
    pressure: float/xarray.DataArray, optional
        atmospheric pressure [kPa]
    elevation: float/xarray.DataArray, optional
        the site elevation [m]
    lat: float/xarray.DataArray, optional
        the site latitude [rad]
    n: float/pandas.Series/xarray.DataArray, optional
        actual duration of sunshine [hour]
    nn: float/pandas.Series/xarray.DataArray, optional
        maximum possible duration of sunshine or daylight hours [hour]
    rso: float/pandas.Series/xarray.DataArray, optional
        clear-sky solar radiation [MJ m-2 day-1]
    a: float, optional
        empirical coefficient for Net Long-Wave radiation [-]
    b: float, optional
        empirical coefficient for Net Long-Wave radiation [-]
    cn: float, optional
        numerator constant [K mm s3 Mg-1 d-1]
    cd: float, optional
        denominator constant [s m-1]
    ea: pandas.Series/float, optional
        actual vapor pressure [kPa]
    albedo: float, optional
        surface albedo [-]
    kab: float, optional
        coefficient derived from as1, bs1 for estimating clear-sky radiation
        [degrees].
    as1: float, optional
        regression constant,  expressing the fraction of extraterrestrial
        reaching the earth on overcast days (n = 0) [-]
    bs1: float, optional
        empirical coefficient for extraterrestrial radiation [-]
    clip_zero: bool, optional
        if True, replace all negative values with 0.
    etype: str, optional
        "os" => ASCE-PM method is applied for a reference surfaces
        representing clipped grass (a short, smooth crop)
        "rs" => ASCE-PM method is applied for a reference surfaces
        representing alfalfa (a taller, rougher agricultural crop),)

    Returns
    -------
    pandas.Series/xarray.DataArray containing the calculated
            potential evaporation [mm d-1].

    Examples
    --------
    >>> et_pm = pm_asce(tmean, wind, rn=rn, rh=rh)

    Notes
    -----
    Following [monteith_1965]_ and [walter_2000]_

    .. math:: PE = \\frac{\\Delta (R_{n}-G)+ \\rho_a c_p K_{min}
        \\frac{e_s-e_a}{r_a}}{\\lambda(\\Delta +\\gamma(1+\\frac{r_s}{r_a}))}

    References
    ----------
    .. [monteith_1965] Monteith, J. L. (1965). Evaporation and environment.
       In Symposia of the society for experimental biology (Vol. 19, pp.
       205-234). Cambridge University Press (CUP) Cambridge.
    .. [walter_2000] Walter, I. A., Allen, R. G., Elliott, R., Jensen, M. E.,
       Itenfisu, D., Mecham, B., ... & Martin, D. (2000). ASCE's standardized
       reference evapotranspiration equation. In Watershed management and
       operations management 2000 (pp. 1-11).
    """
    pressure = calc_press(elevation, pressure)
    gamma = calc_psy(pressure)
    dlt = calc_vpc(tmean)
    ea = calc_ea(tmean=tmean, tmax=tmax, tmin=tmin, rhmax=rhmax,
                 rhmin=rhmin, rh=rh, ea=ea)
    es = calc_es(tmean=tmean, tmax=tmax, tmin=tmin)
    if rn is None:
        rn = get_rn(tmean, rs, check_lat(lat), n, nn, tmax, tmin, rhmax, rhmin,
                    rh, elevation, rso, a, b, ea, albedo, as1, bs1, kab)

    if etype == "rs":
        cn = 1600
        cd = 0.38

    den = dlt + gamma * (1 + cd * wind)
    num1 = (0.408 * dlt * (rn - g)) / den
    num2 = gamma * cn / (tmean + 273) * wind * (es - ea) / den
    pe = num1 + num2
    pe = clip_zeros(pe, clip_zero)
    return pe.rename("PM_ASCE")


def pm(tmean, wind, rs=None, rn=None, g=0, tmax=None, tmin=None, rhmax=None,
       rhmin=None, rh=None, pressure=None, elevation=None, lat=None, n=None,
       nn=None, rso=None, ea=None, a=1.35, b=-0.35, lai=None, croph=0.12,
       r_l=100, r_s=None, ra_method=0, a_sh=1, a_s=1, lai_eff=0, srs=0.0009,
       co2=300, albedo=0.23, kab=None, as1=0.25, bs1=0.5, clip_zero=True):
    """Potential evaporation calculated according to [monteith_1965]_.

    Parameters
    ----------
    tmean: float/xarray.DataArray
        average day temperature [°C]
    wind: float/pandas.Series/xarray.DataArray
        mean day wind speed [m/s]
    rs: float/pandas.Series/xarray.DataArray, optional
        incoming solar radiation [MJ m-2 d-1]
    rn: float/pandas.Series/xarray.DataArray, optional
        net radiation [MJ m-2 d-1]
    g: float/pandas.Series/xarray.DataArray, optional
        soil heat flux [MJ m-2 d-1]
    tmax: float/pandas.Series/xarray.DataArray, optional
        maximum day temperature [°C]
    tmin: float/pandas.Series/xarray.DataArray, optional
        minimum day temperature [°C]
    rhmax: float/pandas.Series/xarray.DataArray, optional
        maximum daily relative humidity [%]
    rhmin: float/pandas.Series/xarray.DataArray, optional
        mainimum daily relative humidity [%]
    rh: float/pandas.Series/xarray.DataArray, optional
        mean daily relative humidity [%]
    pressure: float/xarray.DataArray, optional
        atmospheric pressure [kPa]
    elevation: float/xarray.DataArray, optional
        the site elevation [m]
    lat: float/xarray.DataArray, optional
        the site latitude [rad]
    n: float/pandas.Series/xarray.DataArray, optional
        actual duration of sunshine [hour]
    nn: float/pandas.Series/xarray.DataArray, optional
        maximum possible duration of sunshine or daylight hours [hour]
    rso: float/pandas.Series/xarray.DataArray, optional
        clear-sky solar radiation [MJ m-2 day-1]
    ea: float/pandas.Series/xarray.DataArray, optional
        actual vapor pressure [kPa]
    a: float, optional
        empirical coefficient for Net Long-Wave radiation [-]
    b: float, optional
        empirical coefficient for Net Long-Wave radiation [-]
    lai: float/pandas.Series/xarray.DataArray, optional
        leaf area index [-]
    croph: float/pandas.Series/xarray.DataArray, optional
        crop height [m]
    r_l: pandas.series/float, optional
        bulk stomatal resistance [s m-1]
    r_s: pandas.series/float, optional
        bulk surface resistance [s m-1]
    ra_method: float, optional
        0 => ra = 208/wind
        1 => ra is calculated based on equation 36 in FAO (1990), ANNEX V.
    a_s: float, optional
        Fraction of one-sided leaf area covered by stomata (1 if stomata are 1
        on one side only, 2 if they are on both sides)
    a_sh: float, optional
        Fraction of projected area exchanging sensible heat with the air (2)
    lai_eff: float, optional
        0 => LAI_eff = 0.5 * LAI
        1 => LAI_eff = lai / (0.3 * lai + 1.2)
        2 => LAI_eff = 0.5 * LAI; (LAI>4=4)
        3 => see [zhang_2008]_.
    srs: float/pandas.Series/xarray.DataArray, optional
        Relative sensitivity of rl to Δ[CO2]
    co2: float/pandas.Series/xarray.DataArray, optional
        CO2 concentration [ppm]
    albedo: float, optional
        surface albedo [-]
    kab: float, optional
        coefficient derived from as1, bs1 for estimating clear-sky radiation
        [degrees].
    as1: float, optional
        regression constant,  expressing the fraction of extraterrestrial
        reaching the earth on overcast days (n = 0) [-]
    bs1: float, optional
        empirical coefficient for extraterrestrial radiation [-]
    clip_zero: bool, optional
        if True, replace all negative values with 0.

    Returns
    -------
    pandas.Series/xarray.DataArray containing the calculated
            potential evaporation [mm d-1].

    Examples
    --------
    >>> et_pm = pm(tmean, wind, rn=rn, rh=rh)

    Notes
    -----
    Following [monteith_1965]_, [allen_1998]_, [zhang_2008]_,
        [schymanski_2016]_ and [yang_2019]_.

    .. math:: PE = \\frac{\\Delta (R_{n}-G)+ \\rho_a c_p K_{min}
        \\frac{e_s-e_a}{r_a}}{\\lambda(\\Delta +\\gamma(1+\\frac{r_s}{r_a}))}

    References
    ----------
    .. [monteith_1965] Monteith, J. L. (1965). Evaporation and environment.
       In Symposia of the society for experimental biology (Vol. 19, pp.
       205-234). Cambridge University Press (CUP) Cambridge.
    .. [schymanski_2016] Schymanski, S. J., & Or, D. (2017). Leaf-scale
       experiments reveal an important omission in the Penman–Monteith
       equation. Hydrology and Earth System Sciences, 21(2), 685-706.
    """
    pressure = calc_press(elevation, pressure)
    gamma = calc_psy(pressure)
    dlt = calc_vpc(tmean)
    lambd = calc_lambda(tmean)

    ea = calc_ea(tmean=tmean, tmax=tmax, tmin=tmin, rhmax=rhmax,
                 rhmin=rhmin, rh=rh, ea=ea)
    es = calc_es(tmean=tmean, tmax=tmax, tmin=tmin)

    res_a = calc_res_aero(wind, ra_method=ra_method, croph=croph)
    res_s = calc_res_surf(lai=lai, r_s=r_s, r_l=r_l, lai_eff=lai_eff, srs=srs,
                          co2=co2, croph=croph)
    gamma1 = gamma * a_sh / a_s * (1 + res_s / res_a)

    if rn is None:
        rn = get_rn(tmean, rs, check_lat(lat), n, nn, tmax, tmin, rhmax, rhmin,
                    rh, elevation, rso, a, b, ea, albedo, as1, bs1, kab)

    kmin = 86400  # unit conversion s d-1
    rho_a = calc_rho(pressure, tmean, ea)

    den = lambd * (dlt + gamma1)
    num1 = dlt * (rn - g) / den
    num2 = rho_a * CP * kmin * (es - ea) * a_sh / res_a / den
    pe = num1 + num2
    pe = clip_zeros(pe, clip_zero)
    return pe.rename("Penman_Monteith")


def pm_fao56(tmean, wind, rs=None, rn=None, g=0, tmax=None, tmin=None,
             rhmax=None, rhmin=None, rh=None, pressure=None, elevation=None,
             lat=None, n=None, nn=None, rso=None, a=1.35, b=-0.35, ea=None,
             albedo=0.23, kab=None, as1=0.25, bs1=0.5, clip_zero=True):
    """Potential evaporation calculated according to [allen_1998]_.

    Parameters
    ----------
    tmean: pandas.Series/xarray.DataArray
        average day temperature [°C]
    wind: float/pandas.Series/xarray.DataArray
        mean day wind speed [m/s]
    rs: float/pandas.Series/xarray.DataArray, optional
        incoming solar radiation [MJ m-2 d-1]
    rn: float/pandas.Series/xarray.DataArray, optional
        net radiation [MJ m-2 d-1]
    g: float/pandas.Series/xarray.DataArray, optional
        soil heat flux [MJ m-2 d-1]
    tmax: float/pandas.Series/xarray.DataArray, optional
        maximum day temperature [°C]
    tmin: float/pandas.Series/xarray.DataArray, optional
        minimum day temperature [°C]
    rhmax: float/pandas.Series/xarray.DataArray, optional
        maximum daily relative humidity [%]
    rhmin: float/pandas.Series/xarray.DataArray, optional
        mainimum daily relative humidity [%]
    rh: float/pandas.Series/xarray.DataArray optional
        mean daily relative humidity [%]
    pressure: float/pandas.Series/xarray.DataArray, optional
        atmospheric pressure [kPa]
    elevation: float/xarray.DataArray, optional
        the site elevation [m]
    lat: float/xarray.DataArray, optional
        the site latitude [rad]
    n: float/pandas.Series/xarray.DataArray, optional
        actual duration of sunshine [hour]
    nn: float/pandas.Series/xarray.DataArray, optional
        maximum possible duration of sunshine or daylight hours [hour]
    rso: float/pandas.Series/xarray.DataArray, optional
        clear-sky solar radiation [MJ m-2 day-1]
    a: float, optional
        empirical coefficient for Net Long-Wave radiation [-]
    b: float, optional
        empirical coefficient for Net Long-Wave radiation [-]
    ea: float/pandas.Series/xarray.DataArray, optional
        actual vapor pressure [kPa]
    albedo: float, optional
        surface albedo [-]
    kab: float, optional
        coefficient derived from as1, bs1 for estimating clear-sky radiation
        [degrees].
    as1: float, optional
        regression constant,  expressing the fraction of extraterrestrial
        reaching the earth on overcast days (n = 0) [-]
    bs1: float, optional
        empirical coefficient for extraterrestrial radiation [-]
    clip_zero: bool, optional
        if True, replace all negative values with 0.

    Returns
    -------
    pandas.Series/xarray.DataArray containing the calculated
            potential evaporation [mm d-1].

    Examples
    --------
    >>> et_fao56 = pm_fao56(tmean, wind, rn=rn, rh=rh)

    Notes
    -----
    .. math:: PE = \\frac{0.408 \\Delta (R_{n}-G)+\\gamma \\frac{900}{T+273}
        (e_s-e_a) u_2}{\\Delta+\\gamma(1+0.34 u_2)}

    """
    if tmean is None:
        tmean = (tmax + tmin) / 2
    pressure = calc_press(elevation, pressure)
    gamma = calc_psy(pressure)
    dlt = calc_vpc(tmean)

    gamma1 = (gamma * (1 + 0.34 * wind))

    ea = calc_ea(tmean=tmean, tmax=tmax, tmin=tmin, rhmax=rhmax,
                 rhmin=rhmin, rh=rh, ea=ea)
    es = calc_es(tmean=tmean, tmax=tmax, tmin=tmin)

    if rn is None:
        rn = get_rn(tmean, rs, check_lat(lat), n, nn, tmax, tmin, rhmax, rhmin,
                    rh, elevation, rso, a, b, ea, albedo, as1, bs1, kab)

    den = dlt + gamma1
    num1 = (0.408 * dlt * (rn - g)) / den
    num2 = (gamma * (es - ea) * 900 * wind / (tmean + 273)) / den
    pe = num1 + num2
    pe = clip_zeros(pe, clip_zero)
    return pe.rename("PM_FAO_56")


def priestley_taylor(tmean, rs=None, rn=None, g=0, tmax=None, tmin=None,
                     rhmax=None, rhmin=None, rh=None, pressure=None,
                     elevation=None, lat=None, n=None, nn=None, rso=None,
                     a=1.35, b=-0.35, alpha=1.26, albedo=0.23, clip_zero=True):
    """Potential evaporation calculated according to
        [priestley_and_taylor_1965]_.

    Parameters
    ----------
    tmean: pandas.Series/xarray.DataArray
        average day temperature [°C]
    rs: float/pandas.Series/xarray.DataArray, optional
        incoming solar radiation [MJ m-2 d-1]
    rn: float/pandas.Series/xarray.DataArray, optional
        net radiation [MJ m-2 d-1]
    g: float/pandas.Series/xarray.DataArray, optional
        soil heat flux [MJ m-2 d-1]
    tmax: float/pandas.Series/xarray.DataArray, optional
        maximum day temperature [°C]
    tmin: float/pandas.Series/xarray.DataArray, optional
        minimum day temperature [°C]
    rhmax: float/pandas.Series/xarray.DataArray, optional
        maximum daily relative humidity [%]
    rhmin: float/pandas.Series/xarray.DataArray, optional
        mainimum daily relative humidity [%]
    rh: float/pandas.Series/xarray.DataArray, optional
        mean daily relative humidity [%]
    pressure: float/pandas.Series/xarray.DataArray, optional
        atmospheric pressure [kPa]
    elevation: float/xarray.DataArray, optional
        the site elevation [m]
    lat: float/xarray.DataArray, optional
        the site latitude [rad]
    n: float/pandas.Series/xarray.DataArray, optional
        actual duration of sunshine [hour]
    nn: float/pandas.Series/xarray.DataArray, optional
        maximum possible duration of sunshine or daylight hours [hour]
    rso: float/pandas.Series/xarray.DataArray, optional
        clear-sky solar radiation [MJ m-2 day-1]
    a: float, optional
        empirical coefficient for Net Long-Wave radiation [-]
    b: float, optional
        empirical coefficient for Net Long-Wave radiation [-]
    alpha: float, optional
        calibration coeffiecient [-]
    albedo: float, optional
        surface albedo [-]
    clip_zero: bool, optional
        if True, replace all negative values with 0.

    Returns
    -------
    pandas.Series/xarray.DataArray containing the calculated
            potential evaporation [mm d-1].
    Examples
    --------
    >>> pt = priestley_taylor(tmean, rn=rn, rh=rh)

    Notes
    -----

    .. math:: PE = \\frac{\\alpha_{PT} \\Delta (R_n-G)}
        {\\lambda(\\Delta +\\gamma)}

    References
    ----------
    .. [priestley_and_taylor_1965] Priestley, C. H. B., & TAYLOR, R. J. (1972).
       On the assessment of surface heat flux and evaporation using large-scale
       parameters. Monthly weather review, 100(2), 81-92.

    """
    pressure = calc_press(elevation, pressure)
    gamma = calc_psy(pressure)
    dlt = calc_vpc(tmean)
    lambd = calc_lambda(tmean)

    if rn is None:
        rn = get_rn(tmean, rs, check_lat(lat), n, nn, tmax, tmin, rhmax, rhmin,
                    rh, elevation, rso, a, b, albedo=albedo)

    pe = (alpha * dlt * (rn - g)) / (lambd * (dlt + gamma))
    pe = clip_zeros(pe, clip_zero)
    return pe.rename("Priestley_Taylor")


def kimberly_penman(tmean, wind, rs=None, rn=None, g=0, tmax=None, tmin=None,
                    rhmax=None, rhmin=None, rh=None, pressure=None,
                    elevation=None, lat=None, n=None, nn=None, rso=None,
                    a=1.35, b=-0.35, ea=None, albedo=0.23, clip_zero=True):
    """Potential evaporation calculated according to [wright_1982]_.

    Parameters
    ----------
    tmean: pandas.Series/xarray.DataArray
        average day temperature [°C]
    wind: pandas.Series/xarray.DataArray
        mean day wind speed [m/s]
    rs: pandas.Series/xarray.DataArray, optional
        incoming solar radiation [MJ m-2 d-1]
    rn: pandas.Series/xarray.DataArray, optional
        net radiation [MJ m-2 d-1]
    g: pandas.Series/int/xarray.DataArray, optional
        soil heat flux [MJ m-2 d-1]
    tmax: pandas.Series/xarray.DataArray, optional
        maximum day temperature [°C]
    tmin: pandas.Series/xarray.DataArray, optional
        minimum day temperature [°C]
    rhmax: pandas.Series/xarray.DataArray, optional
        maximum daily relative humidity [%]
    rhmin: pandas.Series/xarray.DataArray, optional
        mainimum daily relative humidity [%]
    rh: pandas.Series/xarray.DataArray, optional
        mean daily relative humidity [%]
    pressure: float/xarray.DataArray, optional
        atmospheric pressure [kPa]
    elevation: float/xarray.DataArray, optional
        the site elevation [m]
    lat: float/xarray.DataArray, optional
        the site latitude [rad]
    n: pandas.Series/float, optional
        actual duration of sunshine [hour]
    nn: pandas.Series/float, optional
        maximum possible duration of sunshine or daylight hours [hour]
    rso: pandas.Series/float, optional
        clear-sky solar radiation [MJ m-2 day-1]
    a: float, optional
        empirical coefficient for Net Long-Wave radiation [-]
    b: float, optional
        empirical coefficient for Net Long-Wave radiation [-]
    ea: float/pandas.Series/xarray.DataArray, optional
        actual vapor pressure [kPa]
    albedo: float, optional
        surface albedo [-]
    clip_zero: bool, optional
        if True, replace all negative values with 0.

    Returns
    -------
    pandas.Series/xarray.DataArray containing the calculated
            potential evaporation [mm d-1].

    Notes
    -----
    Following [oudin_2005]_.

    .. math:: PE = \\frac{\\Delta (R_n-G)+ \\gamma (e_s-e_a) w}
        {\\lambda(\\Delta +\\gamma)}
    .. math:: w =  u_2 * (0.4 + 0.14 * exp(-(\\frac{J_D-173}{58})^2)) +
            (0.605 + 0.345 * exp(-(\\frac{J_D-243}{80})^2))

    References
    ----------
    .. [wright_1982] Wright, J. L. (1982). New evapotranspiration crop
       coefficients. Proceedings of the American Society of Civil Engineers,
       Journal of the Irrigation and Drainage Division, 108(IR2), 57-74.
    """
    pressure = calc_press(elevation, pressure)
    gamma = calc_psy(pressure)
    dlt = calc_vpc(tmean)
    lambd = calc_lambda(tmean)

    ea = calc_ea(tmean=tmean, tmax=tmax, tmin=tmin, rhmax=rhmax,
                 rhmin=rhmin, rh=rh, ea=ea)
    es = calc_es(tmean=tmean, tmax=tmax, tmin=tmin)

    if rn is None:
        rn = get_rn(tmean, rs, check_lat(lat), n, nn, tmax, tmin, rhmax, rhmin,
                    rh, elevation, rso, a, b, ea, albedo)

    j = day_of_year(tmean.index)
    w = wind * (0.4 + 0.14 * exp(-((j - 173) / 58) ** 2) + (
            0.605 + 0.345 * exp((j - 243) / 80) ** 2))

    den = lambd * (dlt + gamma)
    num1 = dlt * (rn - g) / den
    num2 = gamma * (es - ea) * w / den
    pe = num1 + num2
    pe = clip_zeros(pe, clip_zero)
    return pe.rename("Kimberly_Penman")


def thom_oliver(tmean, wind, rs=None, rn=None, g=0, tmax=None, tmin=None,
                rhmax=None, rhmin=None, rh=None, pressure=None, elevation=None,
                lat=None, n=None, nn=None, rso=None, aw=2.6, bw=0.536, a=1.35,
                b=-0.35, lai=None, croph=0.12, r_l=100, r_s=None, ra_method=0,
                lai_eff=0, srs=0.0009, co2=300, ea=None, albedo=0.23,
                clip_zero=True):
    """Potential evaporation calculated according to [thom_1977]_.

    Parameters
    ----------
    tmean: pandas.Series
        average day temperature [°C]
    wind: pandas.Series
        mean day wind speed [m/s]
    rs: pandas.Series, optional
        incoming solar radiation [MJ m-2 d-1]
    rn: pandas.Series, optional
        net radiation [MJ m-2 d-1]
    g: pandas.Series/int, optional
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
    rso: pandas.Series/float, optional
        clear-sky solar radiation [MJ m-2 day-1]
    aw: float, optional
        wind coefficient [-]
    bw: float, optional
        wind coefficient [-]
    a: float, optional
        empirical coefficient for Net Long-Wave radiation [-]
    b: float, optional
        empirical coefficient for Net Long-Wave radiation [-]
    lai: pandas.Series/float, optional
        leaf area index [-]
    croph: pandas.series/float, optional
        crop height [m]
    r_l: pandas.series/float, optional
        bulk stomatal resistance [s m-1]
    r_s: pandas.series/float, optional
        bulk surface resistance [s m-1]
    ra_method: float, optional
        1 => ra = 208/wind
        2 => ra is calculated based on equation 36 in FAO (1990), ANNEX V.
    lai_eff: float, optional
        0 => LAI_eff = 0.5 * LAI
        1 => LAI_eff = lai / (0.3 * lai + 1.2)
        2 => LAI_eff = 0.5 * LAI; (LAI>4=4)
        3 => see [zhang_2008]_.
    srs: float, optional
        Relative sensitivity of rl to Δ[CO2]
    co2: float
        CO2 concentration [ppm]
    ea: float/pandas.Series/xarray.DataArray, optional
        actual vapor pressure [kPa]
    albedo: float, optional
        surface albedo [-]
    clip_zero: bool, optional
        if True, replace all negative values with 0.

    Returns
    -------
    pandas.Series/xarray.DataArray containing the calculated
            potential evaporation [mm d-1].

    Notes
    -----
    Following [oudin_2005]_.

    .. math:: PE = \\frac{\\Delta (R_{n}-G)+ 2.5 \\gamma (e_s-e_a) w}
        {\\lambda(\\Delta +\\gamma(1+\\frac{r_s}{r_a}))}$
        $w=2.6(1+0.53u_2)

    References
    ----------
    .. [thom_1977] Thom, A. S., & Oliver, H. R. (1977). On Penman's equation
       for estimating regional evaporation. Quarterly Journal of the Royal
       Meteorological Society, 103(436), 345-357.
    """
    pressure = calc_press(elevation, pressure)
    gamma = calc_psy(pressure)
    dlt = calc_vpc(tmean)
    lambd = calc_lambda(tmean)

    ea = calc_ea(tmean=tmean, tmax=tmax, tmin=tmin, rhmax=rhmax,
                 rhmin=rhmin, rh=rh, ea=ea)
    es = calc_es(tmean=tmean, tmax=tmax, tmin=tmin)

    res_a = calc_res_aero(wind, ra_method=ra_method, croph=croph)
    res_s = calc_res_surf(lai=lai, r_s=r_s, r_l=r_l, lai_eff=lai_eff, srs=srs,
                          co2=co2, croph=croph)
    gamma1 = gamma * (1 + res_s / res_a)

    if rn is None:
        rn = get_rn(tmean, rs, check_lat(lat), n, nn, tmax, tmin, rhmax, rhmin,
                    rh, elevation, rso, a, b, ea, albedo)

    w = aw * (1 + bw * wind)

    den = lambd * (dlt + gamma1)
    num1 = dlt * (rn - g) / den
    num2 = 2.5 * gamma * (es - ea) * w / den
    pe = num1 + num2
    pe = clip_zeros(pe, clip_zero)
    return pe.rename("Thom_Oliver")


def get_rn(tmean, rs=None, lat=None, n=None, nn=None, tmax=None, tmin=None,
           rhmax=None, rhmin=None, rh=None, elevation=None, rso=None,
           a=1.35, b=-0.35, ea=None, albedo=0.23, as1=0.25, bs1=0.5, kab=None):
    if rs is None:
        rs = calc_rad_sol_in(n, lat, as1=as1, bs1=bs1, nn=nn)
    rns = calc_rad_short(rs=rs, lat=lat, n=n, nn=nn, albedo=albedo, as1=as1,
                         bs1=bs1)  # [MJ/m2/d]
    rnl = calc_rad_long(rs=rs, tmean=tmean, tmax=tmax, tmin=tmin, rhmax=rhmax,
                        rhmin=rhmin, rh=rh, elevation=elevation, lat=lat,
                        rso=rso, a=a, b=b, ea=ea, kab=kab)  # [MJ/m2/d]
    rn = rns - rnl
    return rn


def calculate_all(tmean, wind, rs, elevation, lat, tmax, tmin, rh):
    """Potential evaporation estimated based on all available methods in pyet
     with time series data.

     Parameters
     ----------
     tmean: pandas.Series/xarray.DataArray
         average day temperature [°C]
     wind: float/pandas.Series/xarray.DataArray
         mean day wind speed [m/s]
     rs: float/pandas.Series/xarray.DataArray
         incoming solar radiation [MJ m-2 d-1]
     elevation: float/xarray.DataArray
         the site elevation [m]
     lat: float/xarray.DataArray
         the site latitude [rad]
     tmax: float/pandas.Series/xarray.DataArray
         maximum day temperature [°C]
     tmin: float/pandas.Series/xarray.DataArray
         minimum day temperature [°C]
     rh: float/pandas.Series/xarray.DataArray
         mean daily relative humidity [%]

     Returns
     -------
     pandas.DataFrame containing the calculated potential evaporation

     Examples
     --------
     >>> pe_all = calculate_all(tmean, wind, rs, elevation, lat, tmax=tmax,
                                tmin=tmin, rh=rh)
     """
    pe_df = pandas.DataFrame()
    pe_df["Penman"] = penman(tmean, wind, rs=rs, elevation=elevation,
                             lat=lat, tmax=tmax, tmin=tmin, rh=rh)
    pe_df["FAO-56"] = pm_fao56(tmean, wind, rs=rs, elevation=elevation,
                               lat=lat, tmax=tmax, tmin=tmin, rh=rh)
    pe_df["Priestley-Taylor"] = priestley_taylor(tmean, rs=rs,
                                                 elevation=elevation, lat=lat,
                                                 tmax=tmax, tmin=tmin, rh=rh)
    pe_df["Kimberly-Penman"] = kimberly_penman(tmean, wind, rs=rs,
                                               elevation=elevation,
                                               lat=lat, tmax=tmax, tmin=tmin,
                                               rh=rh)
    pe_df["Thom-Oliver"] = thom_oliver(tmean, wind, rs=rs,
                                       elevation=elevation,
                                       lat=lat, tmax=tmax, tmin=tmin, rh=rh)

    pe_df["Blaney-Criddle"] = blaney_criddle(tmean, lat)
    pe_df["Hamon"] = hamon(tmean, lat=lat, method=1)
    pe_df["Romanenko"] = romanenko(tmean, rh=rh)
    pe_df["Linacre"] = linacre(tmean, elevation, lat, tmax=tmax, tmin=tmin)
    pe_df["Haude"] = haude(tmax, rh)

    pe_df["Turc"] = turc(tmean, rs, rh)
    pe_df["Jensen-Haise"] = jensen_haise(tmean, rs=rs, method=0)
    pe_df["Mcguinness-Bordne"] = mcguinness_bordne(tmean, lat=lat)
    pe_df["Hargreaves"] = hargreaves(tmean, tmax, tmin, lat)
    pe_df["FAO-24"] = fao_24(tmean, wind, rs=rs, rh=rh, elevation=elevation)
    pe_df["Abtew"] = abtew(tmean, rs)
    pe_df["Makkink"] = makkink(tmean, rs, elevation=elevation)
    pe_df["Oudin"] = oudin(tmean, lat=lat)

    return pe_df
