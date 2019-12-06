import numpy as np
import pandas as pd
from input import Rn_calc, vpc_calc, ea_calc, es_calc, psi_calc, \
    lambda_calc, day_of_year, sunset_hangle_calc, Ra_calc


def penman(tmax, tmin, rhmin, rhmax, elevation, latitude, meteoindex, u2,
           solar=None, net=None, a=2.6, b=0.536):
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

    w = a * (1 + b * u2)
    num1 = vpc * rn / (lambd*(vpc+psi))
    num2 = psi*(es-ea) * w / (lambd*(vpc+psi))
    pet = (num1 + num2)
    pet = pd.DataFrame(data=pet, index=meteoindex)
    num1 = pd.DataFrame(data=num1, index=meteoindex)
    num2 = pd.DataFrame(data=num2, index=meteoindex)
    return pet, num1, num2


def penman_monteith(u2, tmax, tmin, rhmin, rhmax, elevation, latitude,
                    meteoindex, solar=None, net=None, rs=69, ra1=208):
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
    aird = 1.225
    airq = 0.001

    ra = ra1/u2
    num1 = vpc * rn / (lambd*(vpc+(psi*(1+rs/ra))))
    num2 = psi*(es-ea) * aird * airq / ra / (lambd*(vpc+(psi*(1+rs/ra))))
    pet = (num1 + num2)
    pet = pd.DataFrame(data=pet, index=meteoindex)
    return pet

def veronika_pm(u2, t, rh, elevation, latitude,
                    meteoindex, solar=None, net=None, g=None):
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

    t = t.to_numpy()
    rh = rh.to_numpy()
    u2 = u2.to_numpy()

    # Inputs
    if solar is None:
        rn = net.to_numpy()
    else:
        solar = solar.to_numpy()
        rn = Rn_calc(solar, t, rh, elevation, latitude,
            meteoindex)
    if g is None:
        g = 0
    else:
        for i in range(0, len(g)):
            if -0.09 < g[i] < 0.27:
                pass
            else:
                g[i] = 0

    psi = 0.0615
    vpc = vpc_calc(t)
    ea = ea_calc(t=t, rh=rh)
    es = es_calc(t=t)
    lambd = lambda_calc(t)
    cn = 37
    cd = np.full(len(lambd), 1)
    for i in range(0, len(rn)):
        if rn[i] > 0:
            cd[i] = 0.34
        else:
            cd[i] = 0.96

    num1 = 0.408 * vpc * (rn-g) + psi*cn*u2*(es-ea)/(t+273)
    num2 = (vpc+(psi*(1+0.24*u2)))
    pet = (num1 / num2)
    pet = pd.DataFrame(data=pet, index=meteoindex)
    num1 = pd.DataFrame(data=num1, index=meteoindex)
    num2 = pd.DataFrame(data=num2, index=meteoindex)
    return pet

def fao_pm(tmax, tmin, rhmin, rhmax, elevation, latitude, meteoindex, u2,
           solar=None, net=None, g=None, cn=900, cd=0.34):
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

    num1 = 0.408 * vpc * (rn-g) + psi*cn*u2*(es-ea)/(ta+273)
    den2 = (vpc+(psi*(1+cd*u2)))
    pet = (num1 / den2)
    pet = pd.DataFrame(data=pet, index=meteoindex)
    return pet

def fao_pm_hourly(u2, t, rh, elevation, latitude,
                    meteoindex, solar=None, net=None, g=None):
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

    t = t.to_numpy()
    rh = rh.to_numpy()
    u2 = u2.to_numpy()

    # Inputs
    if solar is None:
        rn = net.to_numpy()
    else:
        solar = solar.to_numpy()
        rn = Rn_calc(solar, t, rh, elevation, latitude,
            meteoindex)
    if g is None:
        g = 0
    else:
        for i in range(0, len(g)):
            if -0.09 < g[i] < 0.27:
                pass
            else:
                g[i] = 0

    psi = psi_calc(elevation)
    vpc = vpc_calc(t)
    ea = ea_calc(t=t, rh=rh)
    es = es_calc(t=t)
    lambd = lambda_calc(t)
    cn = 37
    cd = np.full(len(lambd), 1)
    for i in range(0, len(rn)):
        if rn[i] > 0:
            cd[i] = 0.34
        else:
            cd[i] = 0.96

    num1 = 0.408 * vpc * (rn-g) + psi*cn*u2*(es-ea)/(t+273)
    num2 = (vpc+(psi*(1+cd*u2)))
    pet = (num1 / num2)
    pet = pd.DataFrame(data=pet, index=meteoindex)
    num1 = pd.DataFrame(data=num1, index=meteoindex)
    num2 = pd.DataFrame(data=num2, index=meteoindex)
    return pet

def priestley_taylor(tmin, tmax, elevation, latitude, rhmin, rhmax, meteoindex,
                     solar=None, net=None, albedo=1.26):
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

    pet = (albedo * vpc * rn)/(lambd*(vpc+psi))
    pet = pd.DataFrame(data=pet, index=meteoindex)
    return pet


def kimberly_penman(meteoindex, tmin, tmax, rhmin, rhmax, elevation, latitude,
                    u2, solar=None, net=None, a=1, b=2):
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

    jd = day_of_year(meteoindex)

    w = (0.4 + 0.14 * np.exp(-((jd-173)/58)**2)) + \
        (0.605 + 0.345 * np.exp(-((jd-243)/80)**2)) * u2
    num1 = vpc * rn / (lambd * (psi + vpc))
    num2 = psi*(es-ea)*w * (a+b*u2) / (lambd * (psi + vpc))
    pet = num1 + num2
    pet = pd.DataFrame(data=pet, index=meteoindex)
    return pet, num1, num2


def hamon(tmin, tmax, meteoindex, latitude):
    """
    DL = day length (h dayK1)
    Ta = air temperature (degrees C)
    """
    # Inputs
    ta = (tmax+tmin)/2
    sunset_hangle = sunset_hangle_calc(meteoindex, latitude)
    dl = 24/np.pi * sunset_hangle  # hours of daylight

    pet = (dl/12)**2 * np.exp(ta/16)
    pet = pd.DataFrame(data=pet, index=meteoindex)
    return pet


def makink(solar, tmin, tmax, elevation, f=1):
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

    pet = f * 1/(lambd * q) * (0.63 * solar * vpc/(vpc+psi) - 14)
    pet = pd.DataFrame(data=pet, index=tmin.index)
    return pet


def hargreaves(meteoindex, tmax, tmin, latitude, a=0.0023, tx=17.8):
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

    pet = a * ra/(q*lambd) * (tmax - tmin)**0.5*(ta + tx)
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


def jensen_haise(tmax, tmin, meteoindex, solar, cr=0.02, tx=4):
    """
    Ra = extraterrestrial radiation (MJ m-2 day-1)
    lambd = latent heat of vaporization (MJ kg-1)
    Ta = air temperature (degrees C)
    q = water density (1000 kg L-1)
    """
    # Inputs
    ta = (tmax+tmin)/2
    lambd = lambda_calc(ta)

    pet = 1 / lambd * cr*(ta-tx)*solar
    pet = pd.DataFrame(data=pet, index=meteoindex)
    return pet


def oudin(tmax, tmin, meteoindex, latitude, k1=5, k2=100):
    """
    Ra = extraterrestrial radiation (MJ m-2 day-1)
    lambd = latent heat of vaporization (MJ kg-1)
    Ta = air temperature (degrees C)
    q = water density (1000 kg L-1)
    """
    # Inputs
    ta = (tmax+tmin)/2
    lambd = lambda_calc(ta)
    ra = Ra_calc(meteoindex, latitude)
    q=1
    petd=[]
    for i in range(0, len(ta)):
        if (ta[i] + k1) < 0:
            petd.append(0)
        else:
            petv = ra[i] *(ta+k1) / (lambd * q * k2)
            petd.append(petv)
    pet = pd.DataFrame(data=petd, index=meteoindex)
    return pet
