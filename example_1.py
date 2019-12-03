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
from etmodul.statistics import RMSE, nanp, R2, MBE

pickle_in = open("data/datad","rb")
datad = pickle.load(pickle_in)

pickle_in1 = open("data/etday","rb")
etday = pickle.load(pickle_in1)
etday[etday>10]=1.102
etday[etday==np.nan]=1

from etmodul.pet import penman, penman_monteith, priestley_taylor, \
kimberly_penman, hamon, makink, hargreaves, blaney_criddle, jensen_haise
from sklearn.linear_model import LinearRegression
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
plt.plot(pen.index, pen)
plt.plot(pm.index, pm)
plt.plot(pt.index, pt)
plt.plot(kp.index, kp)
plt.plot(ham.index, ham)
plt.plot(mak.index, mak)
plt.plot(har.index, har)
plt.plot(bc.index, bc)  # do inverse estimation
plt.plot(jh.index, jh)   # do inverse estimation
plt.legend()

shorts = (pen, pm, pt, kp, ham, mak, har, bc, jh)
names = ("Penman", "Penman-Monteith", "Priestley-Taylor", "Kimberley-Penman",
          "Hamon", "Makink", "Hargreaves", "Blaney-Criddle", "Jensen-Haise")
xaxis = (0,0,0,1,1,1,2,2,2)
yaxis = (0,1,2,0,1,2,0,1,2)
fig, ax = plt.subplots(nrows=3, ncols=3, figsize=(9,9))
for short, x, y, name in zip(shorts,xaxis,yaxis, names):
    ax[x,y].scatter(etday, short, marker=".", c="k")
    ax[x,y].plot(np.arange(0,100), np.arange(0,100), c="r")
    ax[x,y].set_xlim(0, etday.max().values)
    ax[x,y].set_ylim(0, etday.max().values)
    ax[x,y].grid(True, linestyle="-.", alpha=1)
    if y == 0:
        ax[x,y].set_ylabel("PET_obs [mm/day]")
    if x== 2:
        ax[x,y].set_xlabel("PET_sim [mm/day]")
    ax[x,y].text(0.2, 6.3, name + "\n" + \
      "y = " + str(round(model.coef_[0,0],3)) + "*k + " + \
      str(round(model.intercept_[0],3)) + "\n" + \
      "R2 = "+str(round(R2(etday.to_numpy(), short.to_numpy()),3)) + "\n" + \
      "RMSE = "+str(round(RMSE(etday.to_numpy(), short.to_numpy()),3)) + "\n" +\
      "MBE = "+str(round(MBE(etday.to_numpy(), short.to_numpy()),3)), size = 10)
    model = LinearRegression(fit_intercept=True)
    x1=etday.to_numpy()
    y1=short.to_numpy()
    model.fit(x1, y1)
    xfit = etday.to_numpy()
    yfit = model.predict(xfit)
    ax[x,y].plot(xfit, yfit)
    x2 = np.arange(0,100)
plt.tight_layout()
