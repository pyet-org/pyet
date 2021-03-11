from numpy import tan, cos, pi, sin, arccos, clip, maximum
from pandas import to_numeric


def day_of_year(tindex):
    """Day of the year (1-365) based on pandas.series.index

    Parameters
    ----------
    tindex: pandas.Series.index

    Returns
    -------
        array of with ints specifying day of year.

    """
    return to_numeric(tindex.strftime('%j'))


def daylight_hours(tindex, lat):
    """Daylight hours [hour].

    Parameters
    ----------
    tindex: pandas.Series.index
    lat: float
        the site latitude [rad]

    Returns
    -------
        pandas.Series containing the calculated daylight hours [hour]

    Notes
    -----
        Based on equation 34 in [allen_1998]_.

    References
    -----
    .. [allen_1998] Allen, R. G., Pereira, L. S., Raes, D., & Smith, M. (1998).
       Crop evapotranspiration-Guidelines for computing crop water
       requirements-FAO Irrigation and drainage paper 56. Fao, Rome, 300.
       (http://www.fao.org/3/x0490e/x0490e06.htm#TopOfPage)

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

    References
    -----
    .. [allen_1998] Allen, R. G., Pereira, L. S., Raes, D., & Smith, M. (1998).
       Crop evapotranspiration-Guidelines for computing crop water
       requirements-FAO Irrigation and drainage paper 56. Fao, Rome, 300.
       (http://www.fao.org/3/x0490e/x0490e06.htm#TopOfPage)

    """
    return arccos(-tan(sol_dec) * tan(lat))


def sunset_angle_hour(tindex, sol_dec, lat, lz, lm):
    """Sunset hour angle from latitude and solar declination - hourly [rad].

    Parameters
    ----------
    tindex: pandas.Series.index
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
        pandas.Series containing the calculated sunset hour angle - hourly
        [rad]

    Notes
    -----
        Based on equations 29, 30, 31, 32, 33 in [allen_1998]_.

    References
    -----
    .. [allen_1998] Allen, R. G., Pereira, L. S., Raes, D., & Smith, M. (1998).
       Crop evapotranspiration-Guidelines for computing crop water
       requirements-FAO Irrigation and drainage paper 56. Fao, Rome, 300.
       (http://www.fao.org/3/x0490e/x0490e06.htm#TopOfPage)

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

    References
    -----
    .. [allen_1998] Allen, R. G., Pereira, L. S., Raes, D., & Smith, M. (1998).
       Crop evapotranspiration-Guidelines for computing crop water
       requirements-FAO Irrigation and drainage paper 56. Fao, Rome, 300.
       (http://www.fao.org/3/x0490e/x0490e06.htm#TopOfPage)

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

    References
    -----
    .. [allen_1998] Allen, R. G., Pereira, L. S., Raes, D., & Smith, M. (1998).
       Crop evapotranspiration-Guidelines for computing crop water
       requirements-FAO Irrigation and drainage paper 56. Fao, Rome, 300.
       (http://www.fao.org/3/x0490e/x0490e06.htm#TopOfPage)

    """
    return 1 + 0.033 * cos(2. * pi / 365. * j)


def extraterrestrial_r(tindex, lat):
    """Extraterrestrial daily radiation [MJ m-2 d-1].

    Parameters
    ----------
    tindex: pandas.Series.index
    lat: float
        the site latitude [rad]

    Returns
    -------
        pandas.Series containing the calculated extraterrestrial radiation

    Notes
    -----
        Based on equation 21 in [allen_1998]_.

    References
    -----
    .. [allen_1998] Allen, R. G., Pereira, L. S., Raes, D., & Smith, M. (1998).
       Crop evapotranspiration-Guidelines for computing crop water
       requirements-FAO Irrigation and drainage paper 56. Fao, Rome, 300.
       (http://www.fao.org/3/x0490e/x0490e06.htm#TopOfPage)

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
    tindex: pandas.Series.index
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

    References
    -----
    .. [allen_1998] Allen, R. G., Pereira, L. S., Raes, D., & Smith, M. (1998).
       Crop evapotranspiration-Guidelines for computing crop water
       requirements-FAO Irrigation and drainage paper 56. Fao, Rome, 300.
       (http://www.fao.org/3/x0490e/x0490e06.htm#TopOfPage)

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
