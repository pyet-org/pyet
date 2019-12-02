# -*- coding: utf-8 -*-
"""
Created on Mon Sep  2 08:24:23 2019

@author: matevz
"""

import numpy as np
import pandas as pd

def RMSE(obs,sim):
    from math import sqrt
    n = np.count_nonzero(~np.isnan(obs))
    rmse = sqrt(np.nansum(np.power((sim-obs),2)/n))
    return rmse

def nanp (obs, sim):
    from math import sqrt
    
    sim[np.isnan(obs)] = np.nan
    obs[np.isnan(sim)] = np.nan
    
    nad = obs-sim
    pod = obs-np.nanmean(obs)
    
    Ens = 1-(np.nansum(np.power(nad,2))/(np.nansum(np.power(pod,2))))
    return Ens

def R2 (obs, sim):
    from math import sqrt
    
    sim[np.isnan(obs)] = np.nan
    obs[np.isnan(sim)] = np.nan
    
    nad = obs-sim
    obsmean = obs-np.nanmean(obs)
    simmean = sim-np.nanmean(sim)
    
    
    R2 = 1-(np.power((np.nansum(obsmean*simmean)),2))/(np.nansum(np.power(obsmean,2))*np.nansum(np.power(simmean,2)))

    return R2