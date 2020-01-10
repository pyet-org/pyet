import numpy as np
import pandas as pd
from inputs import Rn_calc, vpc_calc, ea_calc, es_calc, psi_calc, \
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
    return pet


def penman_monteith(u2, tmax, tmin, rhmin, rhmax, elevation, latitude,
                    meteoindex, solar=None, net=None, rs1=70, ra1=208, h=None,
                    lai=None, rl=100):
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
    if h is None:
        ra = ra1/u2
    else:
        ra = (np.log((2-0.67*h)/(0.123*h)))*(np.log((2-0.67*h)/(0.0123*h)))\
             /(0.41**2)/(u2)

    if lai is None:
        rs = rs1
    else:
        laief=0.5*lai
        #laief=lai/(0.3*lai+1.2)
        rs = rl/(0.5*laief)

    psi = psi_calc(elevation)
    vpc = vpc_calc(ta)
    pressure=101.3 * ((293 - 0.0065 * elevation) / 293) ** (5.26)
    qa = (pressure)/(1.01*(ta+273)*0.287)
    cp = psi* 0.622* 2.45/pressure
    ea = ea_calc(tmin=tmin, tmax=tmax, rhmin=rhmin, rhmax=rhmax)
    es = es_calc(tmin=tmin, tmax=tmax)
    lambd = lambda_calc(ta)


    num1 = vpc * rn / (lambd*(vpc+(psi*(1+rs/ra))))
    num2 = 86400 * cp * qa * (es-ea) / ra / (lambd*(vpc+(psi*(1+rs/ra))))
    pet = (num1 + num2)
    pet = pd.DataFrame(data=pet, index=meteoindex)
    pet.columns = [0]
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
            meteoindex).to_numpy()
    if g is None:
        g=0

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
    return pet


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

    pet = f * 1/(lambd) * 0.61 * solar * vpc/(vpc+psi) - 0.12
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


def blaney_criddle(tmin, tmax, meteoindex, latitude, k=0.75, tx=8):
    """
    Ta = air temperature (degrees C)
    k= (0.82)(Oudin)
    D = bright sunshine (h day-1)
    mean daily percentage of annual daytime hours (%)
    """
    # Inputs
    ta = (tmax+tmin)/2
    months = np.arange (1, 13)
    daily_p = (0.195, 0.23, 0.27, 0.305, 0.34, 0.355, 0.345, 0.32, 0.28, 0.24,
               0.25, 0.17)
    dl = pd.DataFrame (np.nan, index=meteoindex, columns=["dl"])
    for month, dyp in zip (months, daily_p):
        dl["dl"].loc[(dl.index.month == month)] = dyp
    dl=dl.to_numpy().flatten()
    pet = k * dl * (0.46 * ta + tx)
    pet = pd.DataFrame(data=pet, index=meteoindex)
    return pet, dl


def jensen_haise(tmax, tmin, meteoindex, solar, cr=0.025, tx=-3):
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
    petd=[]
    for i in range(0, len(ta)):
        if (ta[i] + k1) < 0:
            petd.append(0)
        else:
            petv = ra[i] *(ta[i]+k1) / (lambd[i] * k2)
            petd.append(petv)
    pet = pd.DataFrame(data=petd, index=meteoindex)
    return pet

def pm_yang(u2, tmax, tmin, rhmin, rhmax, elevation, latitude,
            meteoindex, solar=None, net=None, co=301, ra1=208,
            rrs=55, srs=0.0009, h=None, lai=None, rl=100):
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

    if h is None:
        ra = ra1/u2
    else:
        ra = (np.log((2-0.67*h)/(0.123*h)))*(np.log((2-0.67*h)/(0.0123*h)))\
             /(0.41**2)/(u2)

    if lai is None:
        rs = rrs * (1+srs * (co-300))
    else:
        rs = (rrs * (1+srs *(co-300)))*0.12*24 /lai
    psi = psi_calc(elevation)
    vpc = vpc_calc(ta)
    pressure = 101.3 * ((293 - 0.0065 * elevation) / 293) ** (5.26)
    qa = (pressure)/(1.01*(ta+273)*0.287)
    cp = psi* 0.622* 2.45/pressure
    ea = ea_calc(tmin=tmin, tmax=tmax, rhmin=rhmin, rhmax=rhmax)
    es = es_calc(tmin=tmin, tmax=tmax)
    lambd = lambda_calc(ta)


    num1 = vpc * rn / (lambd*(vpc+(psi*(1+rs/ra))))
    num2 = 86400 * cp * qa * (es-ea) / ra / (lambd*(vpc+(psi*(1+rs/ra))))
    pet = (num1 + num2)
    pet = pd.DataFrame(data=pet, index=meteoindex)
    pet.columns = [0]
    return pet

def fao_yang(tmax, tmin, rhmin, rhmax, elevation, latitude, meteoindex, u2,
           solar=None, net=None, g=None, cn=900, co=320):
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
            meteoindex).to_numpy()
    if g is None:
        g=0

    psi = psi_calc(elevation)
    vpc = vpc_calc(ta)
    ea = ea_calc(tmin=tmin, tmax=tmax, rhmin=rhmin, rhmax=rhmax)
    es = es_calc(tmin=tmin, tmax=tmax)

    num1 = 0.408 * vpc * (rn-g) + psi*cn*u2*(es-ea)/(ta+273)
    den2 = (vpc+(psi*(1+u2*(0.34+2.4*10**(-4)*(co-300)))))
    pet = (num1 / den2)
    pet = pd.DataFrame(data=pet, index=meteoindex)
    return pet