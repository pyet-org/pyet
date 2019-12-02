# -*- coding: utf-8 -*-
"""
Created on Mon Dec  2 11:29:05 2019

@author: matevz
"""


#datah = pd.read_csv("hourly.txt", sep="\t", parse_dates=True, index_col=0,
#                    dayfirst=True)
#
#datah["Netradiation"][datah["Netradiation"]<0] = 0
#datah["Shortwaveradiation"][datah["Shortwaveradiation"]>1550] = np.nan
#
#
#datamean = datah.resample("d").mean()
#datamax = datah.resample("d").max()
#datamin = datah.resample("d").min()
#datasum = datah.resample("d").sum() * 0.0036
#
#
#data={"tmax": datamax["temperature"], "tmin": datamin["temperature"],
#      "rhmax": datamax["rh"], "rhmin": datamin["rh"], 
#      "solar": datasum["Shortwaveradiation"], "netrad": datasum["Netradiation"],
#      "wind": datamean["wind"], "g": datasum["soillfux"]}
#
#datad = pd.DataFrame(data=data, index = datamean.index)
#datad1 = datad[:]["2015-01-01 00:00:00":"2016-12-31 00:00:00"]
#datad1.plot(subplots=True)
#
#pickle_out = open("datad","wb")
#pickle.dump(datad1, pickle_out)

#eth = pd.read_csv("data/et_daily.txt", sep="\t", parse_dates=True, index_col=0,
#                    dayfirst=True)
#
#ethour = eth.resample("d").sum()
#ethour=ethour[:]["2015-01-01 00:00:00":"2016-12-31 00:00:00"]
#pickle_out = open("etday","wb")
#pickle.dump(ethour, pickle_out)

# example 1
import pandas as pd
import numpy as np
import pickle
import matplotlib.pyplot as plt

pickle_in = open("data/datad","rb")
datad = pickle.load(pickle_in)

pickle_in1 = open("data/etday","rb")
etday = pickle.load(pickle_in1)

from etmodul.pet import penman, penman_monteith, priestley_taylor, \
kimberly_penman, hamon, makink, hargreaves, blaney_criddle, jensen_haise
#inputs
tmax = datad["tmax"]
tmin = datad["tmin"]
rhmax = datad["rhmax"]
rhmin = datad["rhmin"]
elevation = 700
latitude = 47.4920
meteoindex = datad.index
u2 = datad["wind"]
net = datad["netrad"]
solar = datad["solar"]


pen, num1, num2= penman(tmax, tmin, rhmin, rhmax, elevation, latitude, meteoindex, u2, net=net)
pm = penman_monteith(u2, tmax, tmin, rhmin, rhmax, elevation, latitude, meteoindex, net=net)
pt = priestley_taylor(tmin, tmax, elevation, latitude, rhmin, rhmax,meteoindex, net=net)
kp, num1, num2 = kimberly_penman(meteoindex, tmin, tmax, rhmin, rhmax, elevation,latitude, u2, net=net)
ham = hamon(tmin, tmax, meteoindex, latitude)
mak = makink(solar, tmin, tmax, elevation)
har = hargreaves(meteoindex, tmax, tmin, latitude)
bc = blaney_criddle(tmin, tmax, meteoindex, latitude, k=0.7, tx=1.2)
jh = jensen_haise(tmax, tmin, meteoindex, latitude,solar, Ct=0.0068, tx=-10)

fig,ax = plt.subplots()
plt.plot(etday.index, etday)
#plt.plot(pen.index, pen)
#plt.plot(pen.index, num1)
#plt.plot(pen.index, num2)
#plt.plot(pm.index, pm)
#plt.plot(pt.index, pt)
#plt.plot(kp.index, kp)
#plt.plot(ham.index, ham)
#plt.plot(mak.index, mak)
#plt.plot(har.index, har)
#plt.plot(bc.index, bc)  # do inverse estimation
plt.plot(jh.index, jh)   # do inverse estimation
plt.legend()

# determine for bc
from lmfit import minimize, Parameters
params = Parameters()
params.add("k1", value = 0.6, min = 0.45, max = 1.2)
params.add("k2", -5., value = -5., min = -10., max = 10)

#
def get_residual(params,data,tmax, tmin, meteoindex, latitude):
    k1 = params["k1"].value
    k2 = params["k2"].value
    
    bc = jensen_haise(tmin, tmax, meteoindex, latitude, Ct=k1, tx=k2)
    bc=bc.to_numpy()
    return data - bc

out = minimize(get_residual, params, args = (etday.to_numpy(), tmax,tmin, meteoindex, latitude))

# determine for jh
params = Parameters()
params.add("k1", value = 0.002, min = 0.001, max = 0.01)
params.add("k2", value = -5., min = -10., max = 20)

#
def get_residual(params,data,tmax, tmin, meteoindex, latitude, solar):
    k1 = params["k1"].value
    k2 = params["k2"].value
    
    jh = jensen_haise(tmax, tmin, meteoindex, latitude,solar, Ct=k1, tx=k2)
    model = jh.to_numpy()
    return data - model

out = minimize(get_residual, params, args = (etday.to_numpy(), tmax,tmin, meteoindex, latitude, solar))