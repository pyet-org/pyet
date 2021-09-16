"""The rad_utils module contains utility functions for radiation data

"""

from numpy import sqrt, clip

from .meteo_utils import calc_ea, extraterrestrial_r_hour, \
    extraterrestrial_r, daylight_hours

# Stefan Boltzmann constant - hourly [MJm-2K-4h-1]
STEFAN_BOLTZMANN_HOUR = 2.042 * 10 ** -10
# Stefan Boltzmann constant - daily [MJm-2K-4d-1]
STEFAN_BOLTZMANN_DAY = 4.903 * 10 ** -9


def calc_rad_long(rs, tmean=None, tmax=None, tmin=None, rhmax=None,
                  rhmin=None, rh=None, elevation=None, lat=None, rso=None,
                  a=1.35, b=-0.35, ea=None, freq="D"):
    """Net longwave radiation [MJ m-2 d-1].

    Parameters
    ----------
    rs: pandas.Series, optional
        incoming solar radiation [MJ m-2 d-1]
    tmean: pandas.Series, optional
        average day temperature [°C]
    tmax: pandas.Series, optional
        maximum day temperature [°C]
    tmin: pandas.Series, optional
        minimum day temperature [°C]
    rhmax: pandas.Series, optional
        maximum daily relative humidity [%]
    rhmin: pandas.Series, optional
        mainimum daily relative humidity [%]
    rh: pandas.Series, optional
        mean daily relative humidity [%]
    elevation: float, optional
        the site elevation [m]
    lat: float, optional
        the site latitude [rad]
    rso: pandas.Series/float, optional
        clear-sky solar radiation [MJ m-2 day-1]
    a: float, optional
        empirical coefficient for Net Long-Wave radiation [-]
    b: float, optional
        empirical coefficient for Net Long-Wave radiation [-]
    ea: pandas.Series, optional
        actual vapor pressure [kPa]
    freq: string, optional
        "D" => daily estimation
        "H" => hourly estimation

    Returns
    -------
        pandas.Series containing the calculated net longwave radiation

    Notes
    -----
    Based on equation 39 in [allen_1998]_.

    References
    ----------
    """
    if ea is None:
        ea = calc_ea(tmean=tmean, tmax=tmax, tmin=tmin, rhmax=rhmax,
                     rhmin=rhmin, rh=rh)

    if freq == "H":
        if rso is None:
            ra = extraterrestrial_r_hour(tindex=rs.index, lat=lat)
            rso = calc_rso(ra=ra, elevation=elevation)
        solar_rat = clip(rs / rso, 0.3, 1)
        tmp1 = STEFAN_BOLTZMANN_HOUR * (tmean + 273.2) ** 4
    else:
        if rso is None:
            ra = extraterrestrial_r(tindex=rs.index, lat=lat)
            rso = calc_rso(ra=ra, elevation=elevation)
        solar_rat = clip(rs / rso, 0.3, 1)
        if tmax is not None:
            tmp1 = STEFAN_BOLTZMANN_DAY * ((tmax + 273.2) ** 4 +
                                           (tmin + 273.2) ** 4) / 2
        else:
            tmp1 = STEFAN_BOLTZMANN_DAY * (tmean + 273.2) ** 4

    tmp2 = 0.34 - 0.139 * sqrt(ea)  # OK
    tmp2 = clip(tmp2, 0.05, 1)
    tmp3 = a * solar_rat + b  # OK
    return tmp1 * tmp2 * tmp3


def calc_rad_short(rs=None, tindex=None, lat=None, alpha=0.23, n=None,
                   nn=None):
    """Net shortwave radiation [MJ m-2 d-1].

    Parameters
    ----------
    rs: pandas.Series, optional
        incoming solar radiation [MJ m-2 d-1]
    tindex: pandas..DatetimeIndex
    lat: float, optional
        the site latitude [rad]
    alpha: float, optional
        surface albedo [-]
    n: pandas.Series/float, optional
        actual duration of sunshine [hour]
    nn: pandas.Series/float, optional
        maximum possible duration of sunshine or daylight hours [hour]


    Returns
    -------
    pandas.Series containing the calculated net shortwave radiation

    Notes
    -----
    Based on equation 38 in [allen_1998]_.
    """
    if rs is not None:
        return (1 - alpha) * rs
    else:
        return (1 - alpha) * calc_rad_sol_in(tindex, lat, n, nn=nn)


def calc_rad_sol_in(tindex, lat, n, as1=0.25, bs1=0.5, nn=None):
    """Incoming solar radiation [MJ m-2 d-1].

    Parameters
    ----------
    tindex: pandas.DatetimeIndex
    lat: float, optional
        the site latitude [rad]
    n: pandas.Series/float
        actual duration of sunshine [hour]
    as1: float, optional
        regression constant,  expressing the fraction of extraterrestrial
        reaching the earth on overcast days (n = 0) [-]
    bs1: float, optional
        empirical coefficient for extraterrestrial radiation [-]
    nn: pandas.Series/float, optional
        maximum possible duration of sunshine or daylight hours [hour]

    Returns
    -------
    pandas.Series containing the calculated net shortwave radiation

    Notes
    -----
    Based on equation 35 in [allen_1998]_.
    """
    ra = extraterrestrial_r(tindex, lat)
    if nn is None:
        nn = daylight_hours(tindex, lat)
    return (as1 + bs1 * n / nn) * ra


def calc_rso(ra, elevation):
    """Clear-sky solar radiation [MJ m-2 day-1].

    Parameters
    ----------
    ra: pandas.Series, optional
        Extraterrestrial daily radiation [MJ m-2 d-1]
    elevation: float, optional
        the site elevation [m]

    Returns
    -------
    pandas.Series containing the calculated Clear-sky solar radiation

    Notes
    -----
    Based on equation 37 in [allen_1998]_.

    """
    return (0.75 + (2 * 10 ** -5) * elevation) * ra
