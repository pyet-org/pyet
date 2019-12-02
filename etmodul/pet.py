import numpy as np
import pandas as pd
from etmodul.input import Rn_calc, vpc_calc, ea_calc, es_calc, psi_calc, \
    lambda_calc, day_of_year, sunset_hangle_calc, Ra_calc


def penman(tmax, tmin, rhmin, rhmax, elevation, latitude, meteoindex, u2,
           solar=None, net=None):
    """
    u2 = wind speed at 2m
    Rn = net solar radiation (MJ m-2 day-1)
    vpc = Slope of vapor pressure curve (kPa degC-1)
    psi  = psychrometric constant (kPa degC-1)
    lambd = latent heat of vaporization (MJ kg-1)
    q = water density (1000 kg L-1)
    ea = actual vapour pressure (kPa)
    ed = saturation vapour pressure (kPa)
    """
    tmax = tmax.to_numpy()
    tmin = tmin.to_numpy()
    rhmin = rhmin.to_numpy()
    rhmax = rhmax.to_numpy()
    u2 = u2.to_numpy()

    # Inputs
    ta = (tmax+tmin)/2
    if solar is None:
        rn = net.to_numpy()
    else:
        solar = solar.to_numpy()
        rn = Rn_calc(solar, tmax, tmin, rhmin, rhmax, elevation, latitude,
            meteoindex)
    psi = psi_calc(elevation)
    vpc = vpc_calc(ta)
    ea = ea_calc(tmin=tmin, tmax=tmax, rhmin=rhmin, rhmax=rhmax)
    es = es_calc(tmin=tmin, tmax=tmax)
    lambd = lambda_calc(ta)
    q = 1

    w = 2.6 * (1+0.536 * u2)
    num1 = vpc * rn / (lambd*(vpc+psi))
    num2 = psi*(es-ea) * w / (lambd*(vpc+psi))
    pet = (num1 + num2)
    pet = pd.DataFrame(data=pet, index=meteoindex)
    num1 = pd.DataFrame(data=num1, index=meteoindex)
    num2 = pd.DataFrame(data=num2, index=meteoindex)
    return pet, num1, num2


def penman_monteith(u2, tmax, tmin, rhmin, rhmax, elevation, latitude,
                    meteoindex, solar=None, net=None):
    """
    u2 = wind speed at 2m
    Rn = net solar radiation (MJ m-2 day-1)
    vpc = Slope of vapor pressure curve (kPa degC-1)
    psi  = psychrometric constant (kPa degC-1)
    lambd = latent heat of vaporization (MJ kg-1)
    q = water density (1000 kg L-1)
    ea = actual vapour pressure (kPa)
    ed = saturation vapour pressure (kPa)
    """

    tmax = tmax.to_numpy()
    tmin = tmin.to_numpy()
    rhmin = rhmin.to_numpy()
    rhmax = rhmax.to_numpy()
    u2 = u2.to_numpy()

    # Inputs
    ta = (tmax+tmin)/2
    if solar is None:
        rn = net.to_numpy()
    else:
        solar = solar.to_numpy()
        rn = Rn_calc(solar, tmax, tmin, rhmin, rhmax, elevation, latitude,
            meteoindex)
    psi = psi_calc(elevation)
    vpc = vpc_calc(ta)
    ea = ea_calc(tmin=tmin, tmax=tmax, rhmin=rhmin, rhmax=rhmax)
    es = es_calc(tmin=tmin, tmax=tmax)
    lambd = lambda_calc(ta)
    q = 1
    aird = 1.225
    airq = 0.001

    rs = np.full(len(lambd), 69)
    ra = 208/u2
    w = 1500/ra
    one = np.full(len(lambd), 1)
    num1 = vpc * rn / (lambd*(vpc+(psi*(1+rs/ra))))
    num2 = psi*(es-ea) * aird * airq / ra / (lambd*(vpc+(psi*(1+rs/ra))))
    pet = (num1 + num2)
    pet = pd.DataFrame(data=pet, index=meteoindex)
    return pet


def priestley_taylor(tmin, tmax, elevation, latitude, rhmin, rhmax, meteoindex,
                     solar=None, net=None):
    """
    albedo = surface albedo
    Rn = net solar radiation (MJ m-2 day-1)
    vpc = Slope of vapor pressure curve (kPa degC-1)
    psi  = psychrometric constant (kPa degC-1)
    lambd = latent heat of vaporization (MJ kg-1)
    q = water density (1000 kg L-1)
    """
    # Inputs
    ta = (tmax+tmin)/2
    if solar is None:
        rn = net
    else:
        rn = Rn_calc(solar, tmax, tmin, rhmin, rhmax, elevation, latitude,
            meteoindex)
    psi = psi_calc(elevation)
    vpc = vpc_calc(ta)
    lambd = lambda_calc(ta)
    q = 1

    albedo = 1.26
    pet = (albedo * vpc * rn)/(lambd*q*(vpc+psi))
    pet = pd.DataFrame(data=pet, index=meteoindex)
    return pet


def kimberly_penman(meteoindex, tmin, tmax, rhmin, rhmax, elevation, latitude,
                    u2, solar=None, net=None,):
    """
    Jd = Julian day
    u2 = wind speed at 2m
    Rn = net solar radiation (MJ m-2 day-1)
    vpc = Slope of vapor pressure curve (kPa degC-1)
    psi  = psychrometric constant (kPa degC-1)
    lambd = latent heat of vaporization (MJ kg-1)
    q = water density (1000 kg L-1)
    ea = actual vapour pressure (kPa)
    ed = saturation vapour pressure (kPa)
    """
    # Inputs
    ta=(tmax+tmin)/2
    if solar is None:
        rn=net
    else:
        rn = Rn_calc(solar, tmax, tmin, rhmin, rhmax, elevation, latitude,
            meteoindex)
    psi = psi_calc(elevation)
    vpc = vpc_calc(ta)
    ea = ea_calc(tmin=tmin, tmax=tmax, rhmin=rhmin, rhmax=rhmax)
    es = es_calc(tmin=tmin, tmax=tmax)
    lambd = lambda_calc(ta)
    q = 1
    jd = day_of_year(meteoindex)

    w = (0.4 + 0.14 * np.exp(-((jd-173)/58)**2)) + \
        (0.605 + 0.345 * np.exp(-((jd-243)/80)**2)) * u2
    num1 = vpc * rn / (lambd * (psi + vpc))
    num2 = psi*(es-ea)*w / (lambd * (psi + vpc))
    pet = num1 + num2
    pet = pd.DataFrame(data=pet, index=meteoindex)
    return pet, num1, num2


def hamon(tmin, tmax, meteoindex, latitude):
    """
    DL = day length (h dayK1)
    Ta = air temperature (degrees C)
    """
    # Inputs
    ta=(tmax+tmin)/2
    sunset_hangle = sunset_hangle_calc(meteoindex, latitude)
    dl = 24/np.pi * sunset_hangle # hours of daylight

    pet = (dl/12)**2 * np.exp(ta/16)
    pet = pd.DataFrame(data=pet, index=meteoindex)
    return pet


def makink(solar, tmin, tmax, elevation):
    """
    Rg = global short-wave radiation (MJ m-2 day-1)
    vpc = Slope of vapor pressure curve (kPa degC-1)
    psi  = psychrometric constant (kPa degC-1)
    lambd = latent heat of vaporization (MJ kg-1)
    q = water density (1000 kg L-1)
    """
    # Inputs
    ta = (tmax+tmin)/2
    psi = psi_calc(elevation)
    vpc = vpc_calc(ta)
    lambd = lambda_calc(ta)
    q = 1

    pet = 1/(lambd * q) *(0.63* solar * vpc/(vpc+psi))
    pet = pd.DataFrame(data=pet, index=tmin.index)
    return pet


def hargreaves(meteoindex, tmax, tmin, latitude):
    """
    Ra = extraterrestrial radiation (MJ m-2 day-1)
    lambd = latent heat of vaporization (MJ kg-1)
    Ta = air temperature (degrees C)
    q = water density (1000 kg L-1)
    """
    # Inputs
    ta=(tmax+tmin)/2
    lambd = lambda_calc(ta)
    q = 1
    ra = Ra_calc(meteoindex, latitude)

    pet = 0.0023 * ra/(q*lambd) * (tmax - tmin)**0.5*(ta + 17.8)
    pet = pd.DataFrame(data=pet, index=meteoindex)
    return pet


def blaney_criddle(tmin, tmax, meteoindex, latitude, k=0.8, tx=4):
    """
    Ta = air temperature (degrees C)
    k= (0.82)(Oudin)
    D = bright sunshine (h day-1)
    mean daily percentage of annual daytime hours (%)
    """
    # Inputs
    ta=(tmax+tmin)/2
    sunset_hangle = sunset_hangle_calc(meteoindex, latitude)
    dl = 24 / np.pi * sunset_hangle  # hours of daylight

    d = dl/24  # bright sunshine (h day-1)

    pet = k * d * (0.46 * ta + tx)
    pet = pd.DataFrame(data=pet, index=meteoindex)
    return pet


def jensen_haise(tmax, tmin, meteoindex, latitude,solar, Ct=0.02, tx=4):
    """
    Ra = extraterrestrial radiation (MJ m-2 day-1)
    lambd = latent heat of vaporization (MJ kg-1)
    Ta = air temperature (degrees C)
    q = water density (1000 kg L-1)
    """
    # Inputs
    ta = (tmax+tmin)/2
    lambd = lambda_calc(ta)
    q = 1
    ra = Ra_calc(meteoindex, latitude)

#    pet = ra * (ta) / (lambd * q*40)
    pet = Ct*(ta-tx)*solar
    pet = pd.DataFrame(data=pet, index=meteoindex)
    return pet