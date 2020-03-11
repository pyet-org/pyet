import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import pyet as et

meteo = pd.read_csv("data/hydrus_meteo.txt", index_col=0,
                    delim_whitespace=True)
meteo.index = pd.date_range("2015-01-01 00:00:00", "2016-01-01 00:00:00")
wind = meteo["wind"] * 0.01157  # km/d to m/s
tmax = meteo["tmax"]
tmin = meteo["tmin"]
ta = (tmax + tmin) / 2
rh = meteo["rhmean"]
solar = meteo["rad"]  # [MJ/m2/d]
elevation = 145.93
latitude = -17.94 * 3.141592654 / 180

hydrus_pm = pd.read_csv("data/hydrus_pm.txt", index_col=0,
                        delim_whitespace=True)
hydrus_har = pd.read_csv("data/hydrus_har.txt", index_col=0,
                         delim_whitespace=True)

# Penman_monteith from Hydrus
pm = et.pm_hydrus(wind, elevation, latitude, rs=solar, tmax=tmax,
                  tmin=tmin, rh=rh, croph=0.6)
plt.plot(np.linspace(0, np.max(pm), 120), np.linspace(0, np.max(pm), 120))
plt.scatter(pm, hydrus_pm["ET"])

ham = et.hamon(ta, latitude)
plt.scatter(pm, ham)
