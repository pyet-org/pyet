from .penman import calc_lambda

from .utils import extraterrestrial_r


def jensen_haise(tmean, rs=None, cr=0.025, tx=-3, lat=None, method=1):
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
    lat: float
        the site latitude [rad]
    method: float, optional
        1 => after [jensen_allen_2016]
        2 => after [oudin_2005]

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
        $ET = \\frac{C_r(T-T_x)R_s}{\\lambda}$

    References
    -----
    .. [jensen_allen_2016] Task Committee on Revision of Manual 70. (2016).
       Evaporation, evapotranspiration, and irrigation water requirements.
       American Society of Civil Engineers.
    .. [oudin_2005] Oudin, L., Hervieu, F., Michel, C., Perrin, C.,
       Andréassian, V., Anctil, F., & Loumagne, C. (2005). Which potential
       evapotranspiration input for a lumped rainfall–runoff model?:
       Part 2—Towards a simple and efficient potential evapotranspiration model
       for rainfall–runoff modelling. Journal of hydrology, 303(1-4), 290-306.

    """
    lambd = calc_lambda(tmean)
    if method == 1:
        return rs / lambd * cr * (tmean - tx)
    else:
        ra = extraterrestrial_r(tmean.index, lat)
        return ra * (tmean + 5) / 68 / lambd
