import numpy as np


def mk(tmax, tmin, rs, elevation, f=1):
    """Returns evapotranspiration calculated with the Makkink (1957) method.

    Based on equation 6 in Allen et al (1998).

    Parameters
    ----------
    tmax: pandas.Series
        maximum day temperature [°C]
    tmin: pandas.Series
        minimum day temperature [°C]
    rs: pandas.Series
        incoming measured solar radiation [MJ m-2 d-1]
    elevation: float/int
        the site elevation [m]
    f: float/int
        crop coefficient [-]
    Returns
    -------
        pandas.Series containing the calculated evapotranspiration
    Examples
    --------
    >>> mak = et.mk(tmax, tmin, rs, elevation)
    """
    ta = (tmax + tmin) / 2
    lambd = lambda_calc(ta)
    pressure = press_calc(elevation)
    gamma = psy_calc(pressure, lambd)
    dlt = vpc_calc(ta, ea_calc(ta))

    return f * 1 / lambd * 0.61 * rs * dlt / (dlt + gamma) - 0.12


def ea_calc(temperature):
    """
    Saturation Vapour Pressure  (ea)
    From FAO (1990), ANNEX V, eq. 10
    """
    return 0.6108 * np.exp((17.27 * temperature) / (temperature + 237.3))


def vpc_calc(temperature, ea):
    """
    From FAO (1990), ANNEX V, eq. 3
    """
    return 4098 * ea / (temperature + 237.3) ** 2


def press_calc(elevation):
    """
    From FAO (1990), ANNEX V, eq. 6
    """
    return 101.3 * ((293. - 0.0065 * elevation) / 293.) ** 5.253


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
