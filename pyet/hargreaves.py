from numpy import sqrt

from .utils import extraterrestrial_r

from .penman import calc_lambda


def hargreaves(tindex, tmax, tmin, lat):
    """Evaporation calculated according to [hargreaves_samani_1982]_.

        Parameters
        ----------
        tindex: pandas.DatetimeIndex
        tmax: pandas.Series, optional
            maximum day temperature [°C]
        tmin: pandas.Series, optional
            minimum day temperature [°C]
        lat: float, optional
            the site latitude [rad]

        Returns
        -------
        pandas.Series containing the calculated evaporation.

        Examples
        --------
        >>> et_har = hargreaves(tindex, tmax, tmin, lat)

        Notes
        -----
            Based on equation (8-16) in [jensen_allen_2016]_.

        .. math::
        -----
            $ET = 0.0023 \\frac{R_a (T_a+17.8)\\sqrt{(T_{max}-T_{min})}}
            {\\lambda}$

        References
        -----
        .. [hargreaves_samani_1982] Hargreaves, G. H., & Samani, Z. A. (1982).
           Estimating potential evapotranspiration. Journal of the irrigation
           and Drainage Division, 108(3), 225-230.
        .. [jensen_allen_2016] Task Committee on Revision of Manual 70. (2016).
           Evaporation, evapotranspiration, and irrigation water requirements.
           American Society of Civil Engineers.

        """
    tmean = (tmax + tmin) / 2
    lambd = calc_lambda(tmean)
    ra = extraterrestrial_r(tindex=tindex, lat=lat)
    return 0.0023 * (tmean + 17.8) * sqrt(tmax - tmin) * ra / lambd
