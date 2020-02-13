# -*- coding: utf-8 -*-
"""
Created on Mon Sep  2 08:24:23 2019

@author: matevz
"""

import numpy as np
import pandas as pd

def rmse(observed, simulated):
    sim, obs = np.array(simulated), np.array(observed)
    mse = np.nanmean((obs - sim) ** 2)
    rmse = np.sqrt(mse)
    return rmse

def MBE(obs,sim):
    n = np.count_nonzero(~np.isnan(obs))
    mbe = (np.nansum((sim-obs)/n))
    return mbe

def nanp (obs, sim):
    from math import sqrt
    
    sim[np.isnan(obs)] = np.nan
    obs[np.isnan(sim)] = np.nan
    
    nad = obs-sim
    pod = obs-np.nanmean(obs)
    
    Ens = 1-(np.nansum(np.power(nad,2))/(np.nansum(np.power(pod,2))))
    return Ens

def R2 (observed, simulated):
    sim, obs = np.array(simulated), np.array(observed)
    mean_obs = np.nanmean(obs)
    squared_error = np.sum((sim - obs) * (sim - obs))
    squared_error_mean = np.sum((sim - mean_obs) * (sim - mean_obs))
    return 1 - (squared_error/squared_error_mean)

def nash(observed, simulated):
    sim, obs = np.array(simulated), np.array(observed)
    mean_obs = np.nanmean(obs)
    num = np.nansum((obs-sim) ** 2)
    denom = np.nansum((obs-mean_obs)**2)
    return 1 - (num/denom)

def t_stats(observed, simulated):
    n = np.count_nonzero(~np.isnan(observed))
    mbe1 = (np.nansum((simulated - observed) / n))
    rmse1 = rmse(observed, simulated)
    return (((n-1)*(mbe1**2))/(rmse1**2 - mbe1**2))**(0.5)