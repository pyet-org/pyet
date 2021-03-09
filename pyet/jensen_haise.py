def jensen_haise(tmax, tmin, solar, cr=0.025, tx=-3.0):
    """Evapotranspiration calculated with the Jensen and Haise (1963) method.

    Parameters
    ----------
    tmax: pandas.Series
        maximum day temperature [°C]
    tmin: pandas.Series
        minimum day temperature [°C]
    solar: Series
        incoming measured solar radiation [MJ m-2 d-1]
    cr: float
        temperature coefficient [-]
    tx: float
        intercept of the temperature axis [°C]

    Returns
    -------
    pandas.Series
        Series containing the calculated evapotranspiration.

    Examples
    --------
    >>> jh_et = jensen_haise(tmax, tmin, solar)

    Notes
    -----
    Based on equation 6 in Allen et al (1998).

    """
    ta = (tmax + tmin) / 2
    lambd = lambda_calc(temperature=ta)
    return 1 / lambd * cr * (ta - tx) * solar


def lambda_calc(temperature):
    """From FAO (1990), ANNEX V, eq. 1
    """
    return 2.501 - 0.002361 * temperature
