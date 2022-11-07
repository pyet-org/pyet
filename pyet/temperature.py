"""The temperature module contains functions of temperature PET methods

"""

from numpy import exp

from .meteo_utils import daylight_hours, calc_ea, calc_es, calc_e0

from pyet.utils import get_index, check_lat, clip_zeros


def blaney_criddle(tmean, lat, a=-1.55, b=0.96, k=0.65, wind=None, rhmin=None,
                   n=None, nn=None, py=None, method=0, clip_zero=True):
    """Potential evaporation calculated according to
    :cite:t:`blaney_determining_1952`.

    Parameters
    ----------
    tmean: pandas.Series/xarray.DataArray
        average day temperature [°C]
    lat: float/xarray.DataArray, optional
        the site latitude [rad]
    a: float, optional
        calibration coefficient for method 0 [-]
    b: float, optional
        calibration coefficient for method 0 [-]
    k: float, optional
        calibration coefficient for method 1 [-]
    wind: float/pandas.Series/xarray.DataArray, optional
        mean day wind speed [m/s]
    rhmin: float/pandas.Series/xarray.DataArray, optional
        mainimum daily relative humidity [%]
    n: float/pandas.Series/xarray.DataArray, optional
        actual duration of sunshine [hour]
    nn: float/pandas.Series/xarray.DataArray, optional
        maximum possible duration of sunshine or daylight hours [hour]
    py: float/pandas.Series/xarray.DataArray, optional
        percentage of actual day-light hours for the day compared to the
        number of day-light hour during the entire year [-]
    method: float, optional
        0 => Blaney Criddle after :cite:t:`schrodter_hinweise_1985`
        1 => Blaney Criddle after :cite:t:`xu_evaluation_2001`
        2 => FAO-24 Blaney Criddle after :cite:t:`mcmahon_estimating_2013`
    clip_zero: bool, optional
        if True, replace all negative values with 0.

    Returns
    -------
    float/pandas.Series/xarray.DataArray containing the calculated
            potential evaporation [mm d-1].

    Examples
    --------
    >>> et_blaney_criddle = blaney_criddle(tmean, lat)

    Notes
    -----
    Based on equation 6 in :cite:p:`xu_evaluation_2001`.

    .. math:: PE=kp(0.46 * T_a + 8.13)

    """
    index = get_index(tmean)
    if nn is None:
        nn = daylight_hours(index, check_lat(lat))
    if py is None:
        py = nn / (365 * 12) * 100
    if method == 0:
        pe = a + b * (py * (0.457 * tmean + 8.128))
    if method == 1:
        pe = k * py * (0.46 * tmean + 8.13)
    elif method == 2:

        k1 = (0.0043 * rhmin - n / nn - 1.41)
        e0, e1, e2, e3, e4, e5 = (0.81917, -0.0040922, 1.0705, 0.065649,
                                  -0.0059684, -0.0005967)
        bvar = e0 + e1 * rhmin + e2 * n / nn + e3 * wind + e4 * rhmin * n / \
               nn + e5 * rhmin * wind
        pe = k1 + bvar * py * (0.46 * tmean + 8.13)
    pe = clip_zeros(pe, clip_zero)
    return pe.rename("Blaney_Criddle")


def haude(tmean, rh, ea=None, k=1, clip_zero=True):
    """Potential evaporation calculated according to
    :cite:t:`haude_determination_1955`.

    Parameters
    ----------
    tmean: pandas.Series/xarray.DataArray
        temperature at 2pm or maximum dailty temperature [°C]
    rh: float/pandas.Series/xarray.DataArray
        average relative humidity at 2pm [%]
    ea: float/pandas.Series/xarray.DataArray, optional
        actual vapor pressure [kPa]
    k: float, optional
        calibration coefficient [-]
    clip_zero: bool, optional
        if True, replace all negative values with 0.

    Returns
    -------
    float/pandas.Series/xarray.DataArray containing the calculated
            potential evaporation [mm d-1].

    Examples
    --------
    >>> et_haude = haude(tmean, rh)

    Notes
    -----
    Following :cite:t:`haude_determination_1955` and
    :cite:t:`schiff_berechnung_1975`.

    .. math:: PE = f * (e_s-e_a)

    """
    e0 = calc_e0(tmean)
    if ea is None:
        ea = rh * e0 / 100
    # Haude coefficients from :cite:t:`schiff_berechnung_1975`
    fk = [0.27, 0.27, 0.28, 0.39, 0.39, 0.37, 0.35, 0.33, 0.31, 0.29, 0.27,
          0.27]
    index = get_index(tmean)
    f = ([fk[x - 1] for x in index.month] * (tmean / tmean).T).T
    pe = k * f * (e0 - ea) * 10  # kPa to hPa
    pe = clip_zeros(pe, clip_zero)
    return pe.rename("Haude")


def hamon(tmean, lat, k=1, c=13.97, cc=218.527, method=0, clip_zero=True):
    """Potential evaporation calculated according to
    :cite:t:`hamon_estimating_1963`.

    Parameters
    ----------
    tmean: pandas.Series/xarray.DataArray
        average day temperature [°C]
    lat: float/xarray.DataArray
        the site latitude [rad]
    k: float, optional
        calibration coefficient if method = 0 [-]
    c: float, optional
        c is a constant for calculation in mm per day if method = 1.
    cc: float, optional
        calibration coefficient if method = 2 [-].
    method: float, optional
        0 => Hamon after :cite:t:`oudin_which_2005`
        1 => Hamon after equation 7 in :cite:t:`ansorge_performance_2019`
        2 => Hamon after equation 12 in :cite:t:`ansorge_performance_2019`.
    clip_zero: bool, optional
        if True, replace all negative values with 0.

    Returns
    -------
    float/pandas.Series/xarray.DataArray containing the calculated
            potential evaporation [mm d-1].

    Examples
    --------
    >>> et_hamon = hamon(tmean, lat)

    Notes
    -----
    Following :cite:t:`hamon_estimating_1963` and :cite:t:`oudin_which_2005`.

    .. math:: PE = (\\frac{DL}{12})^2 exp(\\frac{T_a}{16})

    """
    index = get_index(tmean)
    # Use transpose to work with lat either as int or xarray.DataArray
    dl = daylight_hours(index, check_lat(lat))
    if method == 0:
        pe = k * (dl / 12) ** 2 * exp(tmean / 16)
    if method == 1:
        # saturated water content after Xu and Singh (2001)
        pt = 4.95 * exp(0.062 * tmean) / 100
        pe = c * (dl / 12) ** 2 * pt
    elif method == 2:
        pe = cc * (dl / 12) * 1 / (tmean + 273.3) * exp(
            (17.26939 * tmean) / (tmean + 273.3))
    pe = clip_zeros(pe, clip_zero)
    return pe.rename("Hamon")


def romanenko(tmean, rh, k=4.5, clip_zero=True):
    """Potential evaporation calculated according to
    :cite:t:`romanenko_computation_1961`.

    Parameters
    ----------
    tmean: float/pandas.Series/xarray.DataArray
        average day temperature [°C]
    rh: float/pandas.Series/xarray.DataArray
        mean daily relative humidity [%]
    k: float, optional
        calibration coefficient [-]
    clip_zero: bool, optional
        if True, replace all negative values with 0.

    Returns
    -------
    float/pandas.Series/xarray.DataArray containing the calculated
            potential evaporation [mm d-1].

    Examples
    --------
    >>> et_romanenko = romanenko(tmean, rh)

    Notes
    -----
    Based on equation 11 in :cite:p:`xu_evaluation_2001`.

    .. math:: PE=4.5(1 + (\\frac{T_a}{25})^2 (1  \\frac{e_a}{e_s})

    """
    ea = calc_ea(tmean=tmean, rh=rh)
    es = calc_es(tmean=tmean)
    pe = k * (1 + tmean / 25) ** 2 * (1 - ea / es)
    pe = clip_zeros(pe, clip_zero)
    return pe.rename("Romanenko")


def linacre(tmean, elevation, lat, tdew=None, tmax=None, tmin=None,
            clip_zero=True):
    """Potential evaporation calculated according to
    :cite:t:`linacre_simple_1977`.

    Parameters
    ----------
    tmean: float/pandas.Series/xarray.DataArray
        average day temperature [°C]
    elevation: float/xarray.DataArray
        the site elevation [m]
    lat: float/xarray.DataArray, optional
        the site latitude [°]
    tdew: float/pandas.Series/xarray.DataArray, optional
        mean dew-point temperature [°C]
    tmax: float/pandas.Series/xarray.DataArray, optional
        maximum day temperature [°C]
    tmin: float/pandas.Series/xarray.DataArray, optional
        minimum day temperature [°C]
    clip_zero: bool, optional
        if True, replace all negative values with 0.

    Returns
    -------
    float/pandas.Series/xarray.DataArray containing the calculated
            potential evaporation [mm d-1].

    Examples
    --------
    >>> et_linacre = linacre(tmean, elevation, lat)

    Notes
    -----
    Based on equation 5 in :cite:p:`xu_evaluation_2001`.

    .. math:: PE = \\frac{\\frac{500 T_m}{(100-A)}+15 (T_a-T_d)}{80-T_a}

    """
    if tdew is None:
        tdew = 0.52 * tmin + 0.6 * tmax - 0.009 * tmax ** 2 - 2
    tm = tmean + 0.006 * elevation
    pe = (500 * tm / (100 - check_lat(lat)) + 15 * (tmean - tdew)) / (
            80 - tmean)
    pe = clip_zeros(pe, clip_zero)
    return pe.rename("Linacre")
