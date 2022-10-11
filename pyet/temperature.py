"""The temeprature module contains functions of temeprature PET methods

"""

from numpy import exp

from .meteo_utils import daylight_hours, calc_ea, calc_es, calc_e0

from .utils import get_index


def blaney_criddle(tmean, lat, a=-1.55, b=0.96, k=0.65, wind=None, rhmin=None,
                   n=None, nn=None, method=0):
    """Potential evaporation calculated according to [blaney_1952]_.

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
    method: float, optional
        0 => Blaney Criddle after [schrodter_2013]_
        1 => Blaney Criddle after [Xu_2001]_
        2 => FAO-24 Blaney Criddle after [mcmahon_2013]_

    Returns
    -------
    float/pandas.Series/xarray.DataArray containing the calculated
            potential evaporation [mm d-1].

    Examples
    --------
    >>> et_blaney_criddle = blaney_criddle(tmean, lat)

    Notes
    -----
    Based on equation 6 in [xu_2001]_.

    .. math:: PE=kp(0.46 * T_a + 8.13)

    References
    ----------
    .. [schrodter_2013] Schrödter, H. (2013). Verdunstung:
        Anwendungsorientierte Meßverfahren und Bestimmungsmethoden.
        Springer-Verlag.
    .. [blaney_1952] Blaney, H. F. (1952). Determining water requirements in
       irrigated areas from climatological and irrigation data.
    .. [xu_2001] Xu, C. Y., & Singh, V. P. (2001). Evaluation and
       generalization of temperature‐based methods for calculating evaporation.
       Hydrological processes, 15(2), 305-319.
    .. [mcmahon_2013] McMahon, T. A., Peel, M. C., Lowe, L., Srikanthan, R.,
        and McVicar, T. R. (2013): Estimating actual, potential, reference crop
        and pan evaporation using standard meteorological data: a pragmatic
        synthesis, Hydrol. Earth Syst. Sci., 17, 1331–1363.
    """
    index = get_index(tmean)
    dl = daylight_hours(index, lat)
    py = dl / (365 * 12) * 100
    if method == 0:
        pet = a + b * (py * (0.457 * tmean + 8.128))
    if method == 1:
        pet = k * py * (0.46 * tmean + 8.13)
    if method == 2:
        if nn is None:
            nn = daylight_hours(index, lat)
        k1 = (0.0043 * rhmin - n / nn - 1.41)
        e0, e1, e2, e3, e4, e5 = (0.81917, -0.0040922, 1.0705, 0.065649,
                                  -0.0059684, -0.0005967)
        bvar = e0 + e1 * rhmin + e2 * n / nn + e3 * wind + e4 * rhmin * n / \
               nn + e5 * rhmin * wind
        pet = k1 + bvar * py * (0.46 * tmean + 8.13)
    return pet.rename("Blaney_Criddle")


def haude(tmean, rh, k=1):
    """Potential evaporation calculated according to [haude]_.

    Parameters
    ----------
    tmean: pandas.Series/xarray.DataArray
        temperature at 2pm or maximum dailty temperature [°C]
    rh: float/pandas.Series/xarray.DataArray
        average relative humidity at 2pm [%]
    k: float, optional
        calibration coefficient [-]

    Returns
    -------
    float/pandas.Series/xarray.DataArray containing the calculated
            potential evaporation [mm d-1].

    Examples
    --------
    >>> et_haude = haude(tmean, rh)

    Notes
    -----
    Following [haude_1955]_ and [schiff_1975]_.

    .. math:: PE = f * (e_s-e_a)

    References
    ----------
    .. [haude_1955] Haude, W. (1955). Determination of evapotranspiration by
        an approach as simple as possible. Mitt Dt Wetterdienst, 2(11).
    .. [schiff_1975] Schiff, H. (1975). Berechnung der potentiellen Verdunstung
        und deren Vergleich mit aktuellen Verdunstungswerten von Lysimetern.
        Archiv für Meteorologie, Geophysik und Bioklimatologie, Serie B, 23(4),
        331-342.
    """
    e0 = calc_e0(tmean)
    ea = rh * e0 / 100
    # Haude coefficients from [schiff_1975]_
    fk = [0.27, 0.27, 0.28, 0.39, 0.39, 0.37, 0.35, 0.33, 0.31, 0.29, 0.27,
          0.27]
    index = get_index(tmean)
    f = ([fk[x - 1] for x in index.month] * (tmean / tmean).T).T
    pe = k * f * (e0 - ea) * 10  # kPa to hPa
    return pe.rename("Haude")


def hamon(tmean, lat, k=1, c=13.97, cc=218.527, method=0):
    """Potential evaporation calculated according to [hamon_1961]_.

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
        0 => Hamon after [oudin_2005]_
        1 => Hamon after equation 7 in [ansorge_2019]_
        2 => Hamon after equation 12 in [ansorge_2019]_.

    Returns
    -------
    float/pandas.Series/xarray.DataArray containing the calculated
            potential evaporation [mm d-1].

    Examples
    --------
    >>> et_hamon = hamon(tmean, lat)

    Notes
    -----
    Following [hamon_1961]_ and [oudin_2005]_.

    .. math:: PE = (\\frac{DL}{12})^2 exp(\\frac{T_a}{16})

    References
    ----------
    .. [hamon_1961] Hamon, W. R. (1963). Estimating potential
       evapotranspiration. Transactions of the American Society of Civil
       Engineers, 128(1), 324-338.
    .. [oudin_2005] Oudin, L., Hervieu, F., Michel, C., Perrin, C.,
       Andréassian, V., Anctil, F., & Loumagne, C. (2005). Which potential
       evapotranspiration input for a lumped rainfall–runoff model?:
       Part 2—Towards a simple and efficient potential evapotranspiration model
       for rainfall–runoff modelling. Journal of hydrology, 303(1-4), 290-306.
    .. [ansorge_2019] Ansorge, L., & Beran, A. (2019). Performance of simple
       temperature-based evaporation methods compared with a time series of pan
       evaporation measures from a standard 20 m 2 tank. Journal of Water and
       Land Development.
    """
    index = get_index(tmean)
    # Use transpose to work with lat either as int or xarray.DataArray
    dl = daylight_hours(index, lat)
    if method == 0:
        pe = k * (dl / 12) ** 2 * exp(tmean / 16)
    if method == 1:
        # saturated water content after Xu and Singh (2001)
        pt = 4.95 * exp(0.062 * tmean) / 100
        pe = c * (dl / 12) ** 2 * pt
        pe = pe.where(tmean > 0, 0)
    if method == 2:
        pe = cc * (dl / 12) * 1 / (tmean + 273.3) * exp(
            (17.26939 * tmean) / (tmean + 273.3))
        pe = pe.where(tmean > 0, 0)
    return pe.rename("Hamon")


def romanenko(tmean, rh, k=4.5):
    """Potential evaporation calculated according to [romanenko_1961]_.

    Parameters
    ----------
    tmean: float/pandas.Series/xarray.DataArray
        average day temperature [°C]
    rh: float/pandas.Series/xarray.DataArray
        mean daily relative humidity [%]
    k: float, optional
        calibration coefficient [-]

    Returns
    -------
    float/pandas.Series/xarray.DataArray containing the calculated
            potential evaporation [mm d-1].

    Examples
    --------
    >>> et_romanenko = romanenko(tmean, rh)

    Notes
    -----
    Based on equation 11 in [xu_2001]_.

    .. math:: PE=4.5(1 + (\\frac{T_a}{25})^2 (1  \\frac{e_a}{e_s})

    References
    ----------
    .. [romanenko_1961] Romanenko, V. A. (1961). Computation of the autumn soil
       moisture using a universal relationship for a large area. Proc. of
       Ukrainian Hydrometeorological Research Institute, 3, 12-25.
    """
    ea = calc_ea(tmean=tmean, rh=rh)
    es = calc_es(tmean=tmean)
    pe = k * (1 + tmean / 25) ** 2 * (1 - ea / es)
    return pe.rename("Romanenko")


def linacre(tmean, elevation, lat, tdew=None, tmax=None, tmin=None):
    """Potential evaporation calculated according to [linacre_1977]_.

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

    Returns
    -------
    float/pandas.Series/xarray.DataArray containing the calculated
            potential evaporation [mm d-1].

    Examples
    --------
    >>> et_linacre = linacre(tmean, elevation, lat)

    Notes
    -----
    Based on equation 5 in [xu_2001]_.

    .. math:: PE = \\frac{\\frac{500 T_m}{(100-A)}+15 (T_a-T_d)}{80-T_a}

    References
    -----
    .. [linacre_1977] Linacre, E. T. (1977). A simple formula for estimating
       evaporation rates in various climates, using temperature data alone.
       Agricultural meteorology, 18(6), 409-424.
    """
    if tdew is None:
        tdew = 0.52 * tmin + 0.6 * tmax - 0.009 * tmax ** 2 - 2
    tm = tmean + 0.006 * elevation
    pe = (500 * tm / (100 - lat) + 15 * (tmean - tdew)) / (80 - tmean)
    return pe.rename("Linacre")
