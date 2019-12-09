# -*- coding: utf-8 -*-
"""
Created on Mon Dec  2 11:29:05 2019

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

meteoh = pd.read_csv("C:\\Matevz_arbeit\\My_papers\\Veronika_papers\\response_of_evapo_non-ranfall_to_climate_change\\daten\\meteo_hourly.txt", sep="\t", parse_dates=True, 
                         index_col=0, dayfirst=True)

meteoh["Short_rad"].plot()

radio = pd.read_csv("C:\\Matevz_arbeit\\My_papers\\Veronika_papers\\response_of_evapo_non-ranfall_to_climate_change\\daten\\2016_2018_radiation.txt", sep="\t", parse_dates=True, 
                         index_col=0, dayfirst=True)
radioh = radio.resample("h").sum()
radio["tPAR_Top"].plot()
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
#------------------------------------------------------------------------------
# kcb determination alfalda hay
kinit = 0.4
kmid = 1.2
kend=1.15
maxh=0.7
#maxroot=1-2
p=0.55

linit=10
ldew=30
lmid=25
llate=10

linit1=5
ldew1=10
lmid1=10
llate1=5


k=pd.DataFrame(np.nan, index=etday.index, columns=["k"])
start1 = ("2015-04-01", "2016-04-01")
end1 = ("2015-10-31", "2016-10-31")
cuts1 = ("2015-05-27","2015-07-29","2015-10-11", 
         "2016-05-30","2016-07-27","2016-10-04")

for start, end in zip(start1,end1):
    start_1=datetime.datetime.strptime(start, "%Y-%m-%d")
    end_1=datetime.datetime.strptime(end, "%Y-%m-%d")
    k["k"][start_1]=kinit
    k["k"][end_1]=kinit
    k["k"][start_1 + datetime.timedelta(days=ldew)]=kmid
        
for cut in cuts1:
    cut_1=datetime.datetime.strptime(cut, "%Y-%m-%d")
    after_cut=cut_1 + datetime.timedelta(days=1)
    k["k"][cut_1]=kmid
    k["k"][after_cut]=kinit
    if cut_1.month == 10:
        pass
    else:
        k["k"][cut_1 + datetime.timedelta(days=ldew1)]=kmid

k["k"]=k["k"].interpolate()

k["k"].plot(marker=".")
k1.plot(marker=".")
start = ("2015-04-01", "2016-04-01", "2017-04-01", "2018-04-01")
cuts = (["2015-05-27","2015-07-29","2015-10-11"],
        ["2016-05-30","2016-07-27","2016-10-04"],
        ["2017-05-29","2017-07-25","2017-10-03"],
        ["2018-05-28","2018-07-23","2018-10-01"])

start1 = ()
end = ()


k1 = etday["GS-2"]/pm[0]
k1[k1>1.6]=np.nan
k1[k1<0.35]=np.nan

#-----------------------------------------------------------------------------



pickle_in = open("data/datad","rb")
datad = pickle.load(pickle_in)

pickle_in1 = open("data/etday","rb")
etday = pickle.load(pickle_in1)
etday[etday>10]=1.102
etday[etday==np.nan]=1

tmax1=tmax["2015-05-27":"2015-05-29","2016-05-27":"2016-05-29"]
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



pt = priestley_taylor(tmin, tmax, elevation, latitude, rhmin, rhmax, meteoindex,
                     net=net)
mak = makink(solar, tmin, tmax, elevation)
bc, days = blaney_criddle(tmin, tmax, meteoindex, latitude, k=0.75)
har = hargreaves(meteoindex, tmax, tmin, latitude)
fao = fao_pm(tmax, tmin, rhmin, rhmax, elevation, latitude, meteoindex, u2,
           net=net)
kp = kimberly_penman(meteoindex, tmin, tmax, rhmin, rhmax, 
                                 elevation,latitude, u2, net=net)
jh = jensen_haise(tmax, tmin, meteoindex, solar)
pm, num1, num2 = penman_monteith(u2, tmax, tmin, rhmin, rhmax, elevation, 
                                 latitude, meteoindex, net=net, h=h,
                                 lai=laids)
pen = penman(tmax, tmin, rhmin, rhmax, elevation, latitude, meteoindex, u2,
             net=net)
oud = oudin(tmax, tmin, meteoindex, latitude, k1=5, k2=100)
ham = hamon(tmin, tmax, meteoindex, latitude)

starts=("2015-01-01", "2015-10-31", "2016-10-31" )
ends = ("2015-04-01", "2016-04-01", "2016-12-31")


for end, start in zip(ends,starts):
    for short in shorts:
        short[:][start:end]=np.nan
    



fig,ax = plt.subplots()
plt.plot(etday.index, etday)
plt.plot(ham.index, ham)
plt.plot(oud.index, oud)
plt.plot(pen.index, pen)
plt.plot(pm.index, pm)
plt.plot(jh.index, jh)
plt.plot(kp.index, kp)
plt.plot(fao.index, fao)
plt.plot(bc.index, bc)
plt.plot(mak.index, mak)
plt.plot(pt.index, pt)
plt.plot(har.index, har)



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
    model = LinearRegression(fit_intercept=True)
    x1=etday.to_numpy()
    y1=short.to_numpy()
    model.fit(x1, y1)    
    ax[x,y].text(0.2, 17.3, name + "\n" + \
      "y = " + str(round(model.coef_[0,0],3)) + "*k + " + \
      str(round(model.intercept_[0],3)) + "\n" + \
      "R2 = "+str(round(R2(etday.to_numpy(), short.to_numpy()),3)) + "\n" + \
      "RMSE = "+str(round(RMSE(etday.to_numpy(), short.to_numpy()),3)) + "\n" +\
      "MBE = "+str(round(MBE(etday.to_numpy(), short.to_numpy()),3)), size = 10)

    xfit = etday.to_numpy()
    yfit = model.predict(xfit)
    ax[x,y].plot(xfit, yfit)
    x2 = np.arange(0,100)
plt.tight_layout()

etday.to_numpy()
