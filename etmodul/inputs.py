import numpy as np
import pandas as pd
def ea_calc(t=None, rh=None, tmin=None, tmax=None, rhmin=None, rhmax=None):
    if t is None and rh is None:
        #saturation vapour pressure at daily minimum temperature [kPa]
        emin = 0.6108 * np.exp((17.27 * tmin) / (tmin + 237.3))
        #saturation vapour pressure at daily maximum temperature [kPa]
        emax = 0.6108 * np.exp((17.27 * tmax) / (tmax + 237.3))
        #actual vapour pressure [kPa]
        #Based on equation 17, page 38 in Allen et al (1998).
        ea = ((emin * rhmax/100) + (emax * rhmin/100)) / 2
    else:
        #actual vapour pressure [kPa]
        #Based on equation 19, page 39 in Allen et al (1998).
        ea = (rh/100) * es_calc(t=t)
    return ea

def es_calc(t=None, tmin=None, tmax=None):
    if t is None:
        #saturation vapour pressure at daily minimum temperature [kPa]
        emin = 0.6108 * np.exp((17.27 * tmin) / (tmin + 237.3))
        #saturation vapour pressure at daily maximum temperature [kPa]
        emax = 0.6108 * np.exp((17.27 * tmax) / (tmax + 237.3))
        #mean saturation vapour pressure for a day, week...
        #Based on equation 12, page 36 in Allen et al (1998).
        return (emax + emin) / 2
    else:
        # mean saturation vapour pressure hour...
        return e0_calc(t)

def Ra_calc(meteoindex, latitude):
    # inverse relative distance Earth-Sun
    # Based on equation 23, page 46 in Allen et al (1998).
    dr = 1 + 0.033 * np.cos(2 * np.pi / 365 * day_of_year(meteoindex))
    # solar declination [rad]
    # Based on equation 24, page 46 in Allen et al (1998).
    sol_dec = 0.409 * np.sin(2 * np.pi / 365 * day_of_year(meteoindex) - 1.39)
    # Calculate sunset hour angle (*Ws*) from latitude and solar
    # declination.
    lat = np.pi / 180 * latitude
    # solar constant = 0.0820 MJ m-2 min-1
    gsc = 0.082
    # Based on equation 21, page 46 in Allen et al (1998).
    omega = sunset_hangle_calc(meteoindex, latitude)
    Ra = 24 * 60 / np.pi * gsc * dr * \
        (omega * np.sin(sol_dec) * np.sin(lat) +
         np.cos(sol_dec) * np.cos(lat) * np.sin(omega))
    return Ra

def sunset_hangle_calc(meteoindex, latitude):
    # inverse relative distance Earth-Sun
    # Based on equation 23, page 46 in Allen et al (1998).
    dr = 1 + 0.033 * np.cos(2 * np.pi / 365 * day_of_year(meteoindex))
    # solar declination [rad]
    # Based on equation 24, page 46 in Allen et al (1998).
    sol_dec = 0.409 * np.sin(2 * np.pi / 365 * day_of_year(meteoindex) - 1.39)
    # Calculate sunset hour angle (*Ws*) from latitude and solar
    # declination.
    lat = np.pi / 180 * latitude
    # solar constant = 0.0820 MJ m-2 min-1
    Gsc = 0.082
    # extraterrestrial radiation [MJ m-2 tunit-1]
    # Based on equation 21, page 46 in Allen et al (1998).
    omega = np.arccos(-np.tan(lat) * np.tan(sol_dec))
    # Based on equation 21, page 46 in Allen et al (1998).
    omega = np.arccos(-np.tan(lat) * np.tan(sol_dec))
    return omega
def e0_calc(t):
    # saturation vapour pressure at the air temperature T [kPa].
    # Based on equation 11, page 36 in Allen et al (1998).
    e0 = 0.6108 * np.exp((17.27 * t) / (t + 237.3))
    return e0


def vpc_calc(temperature):
    #slope of saturation vapour pressure curve at air temperature T [kPa °C-1].
    #Based on equation 13, page 37 in Allen et al (1998).
    vpc = 4098 * e0_calc(temperature) / (temperature + 237.3)**2
    return vpc

def psi_calc(elevation):
    # atmospheric pressure [kPa]. Based on equation 7, page 31 in Allen et al (1998).
    pressure = 101.3 * ((293 - 0.0065 * elevation) / 293) ** (5.26)
    # psychrometric constant [kPa °C-1]. Based on equation 8, page 32 in Allen et al (1998).
    psi = 0.665 * 10 ** (-3) * pressure
    return psi

def day_of_year(meteoindex):
    j = pd.to_numeric(meteoindex.strftime('%j'))
    return j

def Rn_calc(solar, tmax, tmin, rhmin, rhmax, elevation, latitude, meteoindex,
            G = 0.23):
    # clear-sky solar radiation [MJ m-2 tunit-1]
    # Based on equation 37, page 51 in Allen et al (1998).
    r_so = (0.75 + (2 * 10 ** -5) * elevation) * Ra_calc(meteoindex, latitude)
    # net solar or shortwave radiation [MJ m-2 tunit-1]
    # Based on equation 38, page 51 in Allen et al (1998).
    Rns = (1 - G) * solar
    # Based on equation 39, page 52 in Allen et al (1998).
    stef = 4.903e-09
    solar_rat = solar / r_so

    tmp1 = stef * ((tmax + 273.2) ** 4 + (tmin + 273.2) ** 4) / 2
    tmp2 = 0.34 - 0.14 * (ea_calc(tmin=tmin, tmax=tmax, rhmin=rhmin, rhmax=rhmax) ** 0.5)
    tmp3 = 1.35 * solar_rat - 0.35

    # net outgoing longwave radiation [MJ m-2 tunit-1]
    Rnl = tmp1 * tmp2 * tmp3
    # net radiation
    return (Rns - Rnl)

def lambda_calc(temperature):
    lambd = 2.501-(2.361 * 10**-3)*temperature
    return lambd