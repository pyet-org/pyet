import numpy as np


def mk(tmax, tmin, rs, elevation, f=1):
    """Returns evapotranspiration calculated with the Makkink (1957) method.

    Based on Jensen and Haise (2013).

    Parameters
    ----------
    tmax: Series
        maximum day temperature [°C]
    tmin: Series
        minimum day temperature [°C]
    rs: Series
        incoming measured solar radiation [MJ m-2 d-1]
    elevation: float/int
        the site elevation [m]
    f: float/int
        crop coefficient [-]
    Returns
    -------
        Series containing the calculated evapotranspiration
    Examples
    --------
    >>> mak = mk(tmax, tmin, rs, elevation)
    """
    ta = (tmax + tmin) / 2
    pressure = press_calc(elevation)
    gamma = psy_calc(pressure)
    dlt = vpc_calc(ta)

    return f / 2.45 * 0.61 * rs * dlt / (dlt + gamma) - 0.12


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
    temperature: Series
         temperature [degC]
    Returns
    -------
        pandas.Series of saturation vapour pressure at the air temperature
        T [kPa]

    """
    return 0.6108 * np.exp((17.27 * temperature) / (temperature + 237.3))
