"""The rad_utils module contains utility functions for radiation data

"""

from numpy import sqrt, clip
from pandas import Series

from .meteo_utils import calc_ea, extraterrestrial_r_hour, \
    extraterrestrial_r, daylight_hours

# Stefan Boltzmann constant - hourly [MJm-2K-4h-1]
STEFAN_BOLTZMANN_HOUR = 2.042 * 10 ** -10
# Stefan Boltzmann constant - daily [MJm-2K-4d-1]
STEFAN_BOLTZMANN_DAY = 4.903 * 10 ** -9


def calc_rad_long(rs, tmean=None, tmax=None, tmin=None, rhmax=None,
                  rhmin=None, rh=None, elevation=None, lat=None, rso=None,
                  a=1.35, b=-0.35, ea=None, lz=0, lon=0, kab=None, freq="D"):
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
    lz: float, optional
        longitude of the center of the local time zone [expressed as positive
        dergrees west of Greenwich, England]. lz = 0 for Greenwich and 345 deg.
        for Paris (France).
    lon: float, optioonal
        longitude [expressed as positive degrees west of Greenwich, England].
    kab: float, optional
        coefficient that can be derived from the as and bs coefficients of the
        Angstrom formula, where Kab = as + bs, and where Kab represents the
        fraction of extraterrestrial radiation reaching the earth on clear-sky
        days [-]
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
            ra = extraterrestrial_r_hour(tindex=rs.index, lat=lat, lz=lz,
                                         lon=lon)
            rso = calc_rso(ra=ra, elevation=elevation, kab=kab)
        rso[rso == 0] = 0.00001  # prevent division with 0
        rso = Series(rso, rs.index)
        solar_rat = Series(clip(rs / rso, 0.25, 1), rs.index)
        idx_night = rso[rso == 0.00001].index
        for idx, value in solar_rat.iteritems():
            rso17 = solar_rat.iloc[0]
            if "17:00" in str(idx):
                rso17 = value
            if idx in idx_night:
                solar_rat.loc[idx] = rso17
        tmp1 = STEFAN_BOLTZMANN_HOUR * (tmean + 273.16) ** 4
    else:
        if rso is None:
            ra = extraterrestrial_r(tindex=rs.index, lat=lat)
            rso = calc_rso(ra=ra, elevation=elevation, kab=kab)
        solar_rat = clip(rs / rso, 0.3, 1)
        if tmax is not None:
            tmp1 = STEFAN_BOLTZMANN_DAY * ((tmax + 273.16) ** 4 +
                                           (tmin + 273.16) ** 4) / 2
        else:
            tmp1 = STEFAN_BOLTZMANN_DAY * (tmean + 273.16) ** 4

    tmp2 = 0.34 - 0.14 * sqrt(ea)  # OK
    tmp3 = a * solar_rat + b  # OK
    tmp3 = clip(tmp3, 0.05, 1)
    return tmp1 * tmp2 * tmp3


def calc_rad_short(rs=None, tindex=None, lat=None, albedo=0.23, n=None, lz=0,
                   lon=0, nn=None, as1=0.25, bs1=0.5, freq="D"):
    """Net shortwave radiation [MJ m-2 d-1].

    Parameters
    ----------
    rs: pandas.Series, optional
        incoming solar radiation [MJ m-2 d-1]
    tindex: pandas..DatetimeIndex
    lat: float, optional
        the site latitude [rad]
    albedo: float, optional
        surface albedo [-]
    n: pandas.Series/float, optional
        actual duration of sunshine [hour]
    lz: float
        longitude of the center of the local time zone [expressed as positive
        degrees west of Greenwich, England]. In the United States, Lz = 75, 90,
        105 and 120° for the Eastern, Central, Rocky Mountain and Pacific time
        zones, respectively, and Lz = 0° for Greenwich, 345° for Paris
        (France), and 255° for Bangkok (Thailand)
    lon: float
        longitude of the solar radiation measurement site [expressed as
        positive degrees west of Greenwich, England]
    as1: float, optional
        regression constant,  expressing the fraction of extraterrestrial
        reaching the earth on overcast days (n = 0) [-]
    bs1: float, optional
        empirical coefficient for extraterrestrial radiation [-]
    nn: pandas.Series/float, optional
        maximum possible duration of sunshine or daylight hours [hour]
    freq: string, optional
        "D" => daily estimation
        "H" => hourly estimation

    Returns
    -------
    pandas.Series containing the calculated net shortwave radiation

    Notes
    -----
    Based on equation 38 in [allen_1998]_.
    """
    if rs is not None:
        return (1 - albedo) * rs
    else:
        return (1 - albedo) * calc_rad_sol_in(tindex, lat, n, lz=lz, lon=lon,
                                              as1=as1, bs1=bs1, nn=nn,
                                              freq=freq)


def calc_rad_sol_in(tindex, lat, n, lz=0, lon=0, as1=0.25, bs1=0.5, nn=None,
                    freq="D"):
    """Incoming solar radiation [MJ m-2 d-1].

    Parameters
    ----------
    tindex: pandas.DatetimeIndex
    lat: float, optional
        the site latitude [rad]
    n: pandas.Series/float
        actual duration of sunshine [hour]
    lz: float
        longitude of the center of the local time zone [expressed as positive
        degrees west of Greenwich, England]. In the United States, Lz = 75, 90,
        105 and 120° for the Eastern, Central, Rocky Mountain and Pacific time
        zones, respectively, and Lz = 0° for Greenwich, 345° for Paris
        (France), and 255° for Bangkok (Thailand)
    lon: float
        longitude of the solar radiation measurement site [expressed as
        positive degrees west of Greenwich, England]
    as1: float, optional
        regression constant,  expressing the fraction of extraterrestrial
        reaching the earth on overcast days (n = 0) [-]
    bs1: float, optional
        empirical coefficient for extraterrestrial radiation [-]
    nn: pandas.Series/float, optional
        maximum possible duration of sunshine or daylight hours [hour]
    freq: string, optional
        "D" => daily estimation
        "H" => hourly estimation

    Returns
    -------
    pandas.Series containing the calculated net shortwave radiation

    Notes
    -----
    Based on equation 35 in [allen_1998]_.
    """
    if freq == "D":
        ra = extraterrestrial_r(tindex, lat)
    else:
        ra = extraterrestrial_r_hour(tindex, lat, lz, lon)
    if nn is None:
        nn = daylight_hours(tindex, lat)
    return (as1 + bs1 * n / nn) * ra


def calc_rso(ra, elevation, kab=None):
    """Clear-sky solar radiation [MJ m-2 day-1].

    Parameters
    ----------
    ra: pandas.Series, optional
        Extraterrestrial daily radiation [MJ m-2 d-1]
    elevation: float, optional
        the site elevation [m]
    kab: float, optional
        coefficient that can be derived from the as and bs coefficients of the
        Angstrom formula, where Kab = as + bs, and where Kab represents the
        fraction of extraterrestrial radiation reaching the earth on clear-sky
        days [-]

    Returns
    -------
    pandas.Series containing the calculated Clear-sky solar radiation

    Notes
    -----
    Based on equation 37 in [allen_1998]_.

    """
    if kab is None:
        return (0.75 + (2 * 10 ** -5) * elevation) * ra
    else:
        return kab * ra
