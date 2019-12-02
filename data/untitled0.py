# -*- coding: utf-8 -*-
"""
Created on Mon Dec  2 11:29:05 2019

@author: matevz
"""
import pandas as pd
import numpy as np
import pickle

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

# example 1

pickle_in = open("datad","rb")
datad = pickle.load(pickle_in)

from etmodule.pet import penman