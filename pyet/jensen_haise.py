from .penman import calc_lambda


def jensen_haise(tmean, rs, cr=0.025, tx=-3):
    """Evaporation calculated accordinf to [jensen_haise_1963]_.

    Parameters
    ----------
    tmean: pandas.Series, optional
        average day temperature [°C]
    rs: pandas.Series, optional
        incoming solar radiation [MJ m-2 d-1]
    cr: float, optional
        temperature coefficient [-]
    tx: float, optional
        intercept of the temperature axis [°C]

    Returns
    -------
    pandas.Series containing the calculated evaporation.

    Examples
    --------
    >>> jh_et = jensen_haise(tmean, rs)

    Notes
    -----
        Based on equation (K-15) in [jensen_allen_2016]_.

    .. math::
    -----
        E = \\frac{C_r(T-T_x)R_s}{\\lambda}

    References
    -----
    .. [jensen_allen_2016] Task Committee on Revision of Manual 70. (2016).
       Evaporation, evapotranspiration, and irrigation water requirements.
       American Society of Civil Engineers.

    """
    lambd = calc_lambda(tmean)
    return rs / lambd * cr * (tmean - tx)
