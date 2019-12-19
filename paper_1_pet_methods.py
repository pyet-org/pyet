# -*- coding: utf-8 -*-
"""
Created on Thu Dec 19 14:29:50 2019

@author: matevz
"""
import datetime
import pandas as pd
import numpy as np
import pickle
import sys
sys.path.append(r'C:\Matevz_arbeit\etmodul\etmodul')
sys.path
import matplotlib.pyplot as plt
from stats import RMSE, nanp, R2, MBE
from pet import penman, penman_monteith, priestley_taylor, fao_pm, \
kimberly_penman, hamon, makink, hargreaves, blaney_criddle, jensen_haise, oudin
from sklearn.linear_model import LinearRegression
#------------------------------------------------------------------------------
# import lysimeter data
lysih = pd.read_csv("data\\lysimeter_et.txt", sep="\t", parse_dates=True, 
                         index_col=0, dayfirst=True)
lysid = lysih.resample("d").sum()

#------------------------------------------------------------------------------
# import meteo Boku data
meteoboku = pd.read_csv("data\\meteo_hourly.txt", sep="\t", parse_dates=True, 
                         index_col=0, dayfirst=True)

meteoboku["netrad"][meteoboku["netrad"]<0] = 0
start = "2015-01-01 00:00:00"
end = "2017-12-31 23:00:00"

meteoamean = meteoboku[:][start:end].resample("d").mean()
meteomax = meteoboku[:][start:end].resample("d").max()
meteomin = meteoboku[:][start:end].resample("d").min()
meteosum = meteoboku[:][start:end].resample("d").sum() * 0.0036

tmax = meteomax["temperature"][start:end]
tmin = meteomax["temperature"][start:end]
rhmax = meteod["rh"][start:end]
rhmin = meteomin["rh"][start:end]
elevation = 700
latitude = 47.4920
meteoindex = meteoamean[start:end].index
u2 = meteoamean["wind"][start:end]
net = meteosum["netrad"][start:end]
solar = meteosum["Short_rad"][start:end]
#global_radh = radzamg[:][start:end].interpolate()
#global_rad = global_radh.resample("d").sum() * 0.0036
