def jh(tmax, tmin, rs, cr=0.025, tx=-3):
    """Returns evapotranspiration calculated with the Jensen and Haise (1963)
    method.

    Based on equation 6 in Allen et al (1998).

    Parameters
    ----------
    tmax: pandas.Series
        maximum day temperature [°C]
    tmin: pandas.Series
        minimum day temperature [°C]
    rs: pandas.Series
        incoming measured solar radiation [MJ m-2 d-1]
    cr: float/int
        temperature coefficient [-]
    tx: float/int
        intercept of the temperature axis [°C]
    Returns
    -------
        pandas.Series containing the calculated evapotranspiration
    Examples
    --------
    >>> jh = et.jh(tmax, tmin, rs)
    """
    ta = (tmax + tmin) / 2
    lambd = lambda_calc(ta)
    return 1 / lambd * cr * (ta - tx) * rs


def lambda_calc(temperature):
    """
    From FAO (1990), ANNEX V, eq. 1
    """
    return 2.501 - 0.002361 * temperature
