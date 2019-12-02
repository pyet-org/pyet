# -*- coding: utf-8 -*-
"""
Created on Mon Sep  2 08:31:08 2019

@author: matevz
"""
def et_val_hourly(lista):
    lista_copy = lista.copy()
    import numpy as np
    ps = 0.6108 * np.exp(17.27 * lista_copy["T"]/(lista_copy["T"]+237.3))
    pa = ps * lista_copy["rh"]/100
    
    kzp = 4098 * ps / pow((lista_copy["T"] + 237.3),2)
    return lista_copy["T"], ps, pa, lista_copy["wind"],lista_copy["netrad"], kzp


def et_val(lista):
    import numpy as np
    import pandas as pd
    #Naredi kopijo uporabljenih vhodnih podatkov. S tem se zavarujemo, da ne pride do sprememb v izvorni datoteki vhodnih podatkov.
    lista_copy = lista.copy()
    #Uporabljene konstante
    #------------------------------------------------
    #Izracun dolgovalovnega sevanja
    #Vrednost Steffan Boltzmanove konstante
    stef = 4.903*10**-9
    #Izračun povprečne Temperature
#    lista_copy["T_max"] = lista_copy["T_max"] + 3
#    lista_copy["T_min"] = lista_copy["T_min"] + 3
    T_pov = (lista_copy["T_max"] + lista_copy["T_min"]) / 2
    #Izračun krivulje zasičenega parnega tlaka (kzp)
    kzp = 4098 * (0.6108 * np.exp((17.27 * T_pov) / (T_pov + 237.3))) / (T_pov + 237.3)**2
    #Izracun zasicenega parnega tlaka pri maks dnevni (esmax)
    psmax = 0.6108 * np.exp(17.27 * lista_copy["T_max"] / (lista_copy["T_max"] + 237.3))
    #Izracun zasicenega parnega tlaka pri min dnevni  (pzmin)
    psmin = 0.6108 * np.exp(17.27 * lista_copy["T_min"] / (lista_copy["T_min"] + 237.3))
    #Izracun dejanskega parnega tlaka (ea)
    pa1 = lista_copy["rh_max"]* psmin / 100
    pa2 = lista_copy["rh_min"] * psmax / 100
    pa = (pa1 + pa2) / 2
    #Izracun srednje vrednosti zasičenega parnega tlaka (pzs)
    ps = (psmax + psmin) / 2
    #Korekcija vetra glede na višino merjenja
    veter = lista_copy["wind"] #* 4.87 / np.log(67.8 * 1 - 5.42)
    #---------------------------------------------------------
    #Izracun neto sevanja na površini rastline (Rn)
    #Dolocitev dneva v letu
    J = pd.to_numeric(lista_copy.index.strftime('%j'))
    #Izračun relativne inverzne razdalja (dr)
    dr = 1 + 0.33 * np.cos(2 * np.pi / 365 * J)
    #Izračun sončnega odklona SO
    so = 0.409 * np.sin(2 * np.pi / 365 * J - 1.39)
    #Izracun zemljepisne sirine pri 16.039875
    zs = np.pi / 180 * 47.49
    #Določitev kota soncnega zahoda sz
    sz = np.arccos(-np.tan(so) * np.tan(zs))
    #Ekstrateresticnega sevanje (Ra)
    Ra = 24 * 60 / np.pi * 0.082 * dr * (sz * np.sin(zs) * np.sin(so) + np.cos(zs) *
    np.cos(so) * np.sin(sz))
    #Izracun radiacija jasnega neba (rso), nadmorska visina 198
    Rso = (0.75 + (2 * 10**-5) * 198) * Ra
    #Izracun kratkovalovnega sevanja Rs - reflekcisjki koeficient/albedo je za referenčno travnato površino enak 0.23
    Rns = (1-0.23) * lista_copy["netr"]
    #Izracun neto sevanja na površini rastline (Rnl)
    Rnl = np.absolute((stef * ((lista_copy["T_max"] + 273.2)**4 + (lista_copy["T_min"] + 273.2)**4) / \
    2 * (0.34 - 0.14 * (pa**0.5)) * (1.35 * (lista_copy["netr"] / Rso) - 0.35)))

    #------------------------------Ra na novo
    sd = 0.409 * np.sin(2 * np.pi / 365 * J - 1.39)
    lat=np.pi / 180 * 47.49
    ws = np.arccos(-np.tan(lat) * np.tan(sd))
    dr = 1 + 0.33 * np.cos(2 * np.pi / 365 * J)
    #extraterrestrial radiation [MJ m-2 day-1]
    Ra = 24 * 60 / np.pi * 0.082 * dr * (ws * np.sin(lat) * np.sin(sd)
    + np.cos(lat) * np.cos(sd) * np.sin(ws))

    #Izračun neto sevanja
    Rn = Rns - Rnl

    for i in range (0,len(Rn)):
        if Rn[i] <= 0:
            Rn[i] = 0.001
        else:
            pass
    return T_pov, lista_copy["T_max"],lista_copy["T_min"],ps, pa, veter,lista_copy["netr"], lista_copy["rh_max"],lista_copy["rh_min"],kzp, Ra
