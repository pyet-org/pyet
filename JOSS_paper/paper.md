---
title: 'PyEt: A Python package for estimating reference and potential evaporation'
tags:
  - Python
  - Evaporation
  - Water balance
  - Hydrology
  - Evapotranspiration
authors:
  - name: Matevž Vremec
    orcid: 0000-0002-2730-494X
    affiliation: 1 
  - name: Raoul Collenteur
    orcid: 0000-0001-7843-1538
    affiliation: 1
affiliations:
 - name: Institute of Earth Sciences, NAWI Graz Geocenter Austria, University of Graz, Austria
   index: 1
date: 2021
bibliography: paper.bib
---

# Summary - write

The evaporation (ET) of water from land surfaces to the atmosphere is a major component of the 
water cycle and accurate estimates of the flux are essential to the water and agricultural sector.  
As direct observations of ET still present a diffulty, ET is often estimated from available 
meteorological data using empirical formulas. There exist a variety of such formulas, where most 
tend to estimate the potential evaporation (PET), which is the maximum amount of water that would 
be evaporated if enough water were available. Consistently implementing these methods is difficult 
due to the significant diversity in the level of input data, process representation, assumptions 
and data requirements. The goal of `PyEt` is to provide a Python package that includes several methods 
for the estimation of PET, is well documented and simple to use. `PyEt` currently includes eighteen 
different methods to estimate PET and various methods to estimate surface and aerodynamic resistance. 
The package allows users to compute and compare potential evaporation estimates using different 
approaches with minimum effort. The structure of the package allows the implementation of the methods 
in other hydrological models or sensitivity analysis.

# Methods

`PyEt` includes 18 well-known PET methods, decribed in ([@xu2000evaluation{, [@xu2001evaluation], 
[@oudin2005potential], [@jensen1990evapotranspiration]), to estimate potential evaporation at a single 
location using input data at daily resolution.

| Classification | Common method name        | Data needed | PyEt Method        | Reference                      |
|----------------|---------------------------|-------------|--------------------|--------------------------------|
| Combination    | Penman                    | RH, T, U, D |`penman`            |[@penman1948natural]            |
|                | Penman-Monteith           | RH, T, U, D |`pm`                |[@monteith1965evaporation]      |
|                | FAO-56                    | RH, T, U, D |`pm_fao56`          |[@allen1998crop]                |
|                | Priestley-Taylor          | T, D        |`priestley_taylor`  |[@priestley1972assessment]      |
|                | Kimberly-Penman           | RH, T, U, D |`kimberly_penman`   |[@wright1982new]                |
|                | Thom-Oliver               | RH, T, U, D |`thom_oliver`       |[@thom1977penman]               |
| Temperature    | Blaney–Criddle            | T, D        |`blaney_criddle`    |[@blaney1952determining]        |
|                | Hamon                     | T           |`hamon`             |[@hamon1963estimating]          |
|                | Romanenko                 | RH, T       |`romanenko`         |[@xu2001evaluation]             |
|                | Linacre                   | T           |`linacre`           |[@linacre1977simple]            |
| Radiation      | Turc                      | T, D        |`turc`              |[@xu2001evaluation]             |
|                | Jensen–Haise              | T, D        |`jensen_haise`      |[@jensen1963estimating]         |
|                | McGuinness–Bordne         | T, D        |`mcguinness_bordne` |[@mcguinness1972comparison]     |
|                | Hargreaves                | T           |`hargreaves`        |[@hargreaves1982estimating]     |
|                | Doorenbos–Pruitt (FAO-24) | RH, T, U, D |`fao_24`            |[@jensen1990evapotranspiration] |
|                | Abtew                     | T, D        |`abtew`             |[@abtew1996evapotranspiration]  |
|                | Makkink                   | T, D        |`makkink`           |[@makkink1957testing]           |
|                | Oudin                     | T           |`oudin`             |[@oudin2005potential]           |
T, Temperature; U, Wind Speed; D, Radiation; RH, Relative Humidity. Adapted from [@oudin2005potential].

# Example application

This example shows how `PyEt` can be used to calculate potential evaporation using various PET methods. 
The PET in the example is estimated for the city of Maribor (Slovenia), using online available data from the
Slovenian Environmental Agency (ARSO).  
+

``` python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import pyet

# Import meteorological data 
meteo = pd.read_csv("data/meteo_maribor_2017.csv", parse_dates=True, index_col="Date")
tmax, tmin, rh, wind, rs = [meteo[col] for col in meteo.columns]
elevation = 279
lat = 46.5678 * np.pi / 180

et_penman = pyet.pm(wind, Rs=Rs, elevation=elevation, lat=lat, tmax=tmax, tmin=tmin, rh=rh)
et_pt = pyet.priestley_taylor(wind, rs=rs, elevation=elevation, lat=lat, tmax=tmax, tmin=tmin, rh=rh)
et_mak = pyet.makkink(Rs, tmax=tmax, tmin=tmin, elevation=elevation)
et_hamon = pyet.hamon(tmean.index, tmean, lat)
```
, where $R_s$ is the incoming solar radiation [MJ m-2 d-1], $elevation$ the site elevation [m], 
$lat$ the the site latitude [rad], $tmax$ and $tmin$ the maximum and minimum temperature [°C] and 
$rh$ the mean relative humidity [%].

![Daily potential evaporation for Maribor (Slovenia) estimated according to [@monteith1965evaporation], 
[@priestley1972assessment], [@makkink1957testing] and [@hamon1963estimating]](Figure1.png)

# Concluding remarks

This paper presents `PyEt`, a Python package to estimate potential evaporation (PET) from available 
meteorological data. Using `PyEt`, users can estimate PET using 18 different methods with only one line
 of Python-code. At this stage, the methods are implemented only for 1D data (e.g. time series data), 
 while further developments will focus on enabling 2D and 3D data input (Numpy Arrays, Array, and NetCDF files). 
 The authors believe that `PyEt` makes a great contribution to the hydrology, meteorology and agricultural 
 community, enabling a simple PET estimation and comparison between different PET estimates. `PyEt` methods 
 and estimates can be further used in hydrological models and sensitivity analyses. The `PyEt` design enables
 simple extension of the software with new capabilities and new model options, the authors therefore welcome 
 new code contribution and suggestions from the community to improve the overall package performance.  `PyEt`
 is free and open-source software under the LGPL-3.0 license and is available at
 http://www.github.com/phydrus/PyEt.  Full documentation is available on ReadTheDocs
 (https://pyet.readthedocs.io). The notebook and data necessary to reproduce the ﬁgure in this manuscript 
 are available through PyEt GitHub page.

# Acknowledgements

# References