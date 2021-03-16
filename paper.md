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

# Summary

The evaporation (ET) of water from land surfaces to the atmosphere is a major component of the water 
cycle and accurate estimates of the flux are essential to the water and agricultural sector. As 
direct measurement of ET is difficult, the evaporation flux is often estimated from available 
meteorological data using empirical formulas. There exist a variety of such formulas, where most 
tend to estimate the potential evaporation (PET), which is the maximum amount of water that would 
be evaporated if enough water were available. Consistently implementing these methods is difficult 
due to the significant diversity in the level of input data, process representation, assumptions 
and data requirements. The goal of `PyEt` is to provide a Python package that includes a wide range 
of methods for the estimation of PET, is fully documented and easy to use. Currently, `PyEt` includes 
eighteen different methods to estimate PET and various methods to estimate surface and aerodynamic 
resistance (see Table below). The package allows users to compute and compare potential evaporation 
estimates using different approaches with minimum effort. The structure of the package allows the 
implementation of the methods in other hydrological models or sensitivity analysis.

# Statement of Need
Provide the statement of need here.

| Common method name  | Data needed | PyEt Method        | Reference                      |
|---------------------|-------------|--------------------|--------------------------------|
| Penman                    | RH, T, U, D |`penman`            |[@penman1948natural]            |
| Penman-Monteith           | RH, T, U, D |`pm`                |[@monteith1965evaporation]      |
| FAO-56                    | RH, T, U, D |`pm_fao56`          |[@allen1998crop]                |
| Priestley-Taylor          | T, D        |`priestley_taylor`  |[@priestley1972assessment]      |
| Kimberly-Penman           | RH, T, U, D |`kimberly_penman`   |[@wright1982new]                |
| Thom-Oliver               | RH, T, U, D |`thom_oliver`       |[@thom1977penman]               |
| Blaney–Criddle            | T, D        |`blaney_criddle`    |[@blaney1952determining]        |
| Hamon                     | T           |`hamon`             |[@hamon1963estimating]          |
| Romanenko                 | RH, T       |`romanenko`         |[@xu2001evaluation]             |
| Linacre                   | T           |`linacre`           |[@linacre1977simple]            |
| Turc                      | T, D        |`turc`              |[@xu2001evaluation]             |
| Jensen–Haise              | T, D        |`jensen_haise`      |[@jensen1963estimating]         |
| McGuinness–Bordne         | T, D        |`mcguinness_bordne` |[@mcguinness1972comparison]     |
| Hargreaves                | T           |`hargreaves`        |[@hargreaves1982estimating]     |
| Doorenbos–Pruitt (FAO-24) | RH, T, U, D |`fao_24`            |[@jensen1990evapotranspiration] |
| Abtew                     | T, D        |`abtew`             |[@abtew1996evapotranspiration]  |
| Makkink                   | T, D        |`makkink`           |[@makkink1957testing]           |
| Oudin                     | T           |`oudin`             |[@oudin2005potential]           |

T, Temperature; U, Wind Speed; D, Radiation; RH, Relative Humidity. Adapted from [@oudin2005potential].

# Example application

This example shows how `PyEt` can be used to compute potential evaporation using different methods. 
The potential evaporation is estimated for the city of Maribor (Slovenia), using online available 
data from the Slovenian Environmental Agency (ARSO). The potential evaporation is computed using 
four different methods (Penman, Priestley-Taylor, Makking, Hammon) and plotted for comparison.

``` python
import pandas as pd
import pyet

# Import meteorological data 
meteo = pd.read_csv("data/meteo_maribor_2017.csv", parse_dates=True, 
                    index_col="Date")
tmax, tmin, rh, wind, rs = [meteo[col] for col in meteo.columns]
lat = 46.5678 * np.pi / 180 

et_penman = pyet.pm(wind, rs=rs, elevation=279, lat=lat, tmax=tmax, 
					tmin=tmin, rh=rh)
et_pt = pyet.priestley_taylor(wind, rs=rs, elevation=279, lat=lat, 
							  tmax=tmax, tmin=tmin, rh=rh)
et_mak = pyet.makkink(rs, tmax=tmax, tmin=tmin, elevation=279)
et_hamon = pyet.hamon(tmean.index, tmean, lat)
```
In the code above $R_s$ is the incoming solar radiation [MJ m-2 d-1], $elevation$ the site elevation [m], 
$lat$ the the site latitude [rad], $tmax$ and $tmin$ the maximum and minimum temperature [°C] and 
$rh$ the mean relative humidity [%].

![Daily potential evaporation for Maribor (Slovenia) estimated according to [@monteith1965evaporation], 
[@priestley1972assessment], [@makkink1957testing] and [@hamon1963estimating]](Figure1.png)

# Concluding remarks

This paper presents `PyEt`, a Python package to estimate potential evaporation from available 
meteorological data. Using `PyEt`, users can estimate PET using 18 different methods with only a few lines
of Python-code. At this stage (PyEt v1.0), the methods are implemented for 1D data (e.g. time series data).
Further developments will focus on enabling 2D and 3D data input (Numpy Arrays, Array, and NetCDF files).
The authors believe that `PyEt` will a valuable contribution to the hydrology, meteorology and agricultural 
communities, enabling a simple PET estimation and comparison between different PET estimates. `PyEt` methods 
and estimates can be further used in hydrological models and sensitivity analyses. The `PyEt` design enables
simple extension of the software with new capabilities and new model options. The authors warmly welcome 
code contributions, bug reports, and feedback from the community to improve further improve the package.
`PyEt` is free and open-source software under the LGPL-3.0 license and is available at
http://www.github.com/phydrus/PyEt. Full documentation is available on ReadTheDocs (https://pyet.readthedocs.io). 
The notebook and data necessary to reproduce the Figure in this manuscript are available through PyEt GitHub page.

# Acknowledgements
The first author is funded by XX. The second authors is funded by the Austrian Science Fund (FWF) under Research 
Grant W1256 (Doctoral Programme Climate 525 Change: Uncertainties, Thresholds and Coping Strategies).

# References