---
title: '`pyet`: A Python package for estimating evaporation'
tags:
  - Python
  - Evaporation
  - Water balance
  - Hydrology
authors:
  - name: Matevž Vremec
    orcid: 0000-0002-2730-494X
    affiliation: 1 
  - name: Raoul A. Collenteur
    orcid: 0000-0001-7843-1538
    affiliation: 1
affiliations:
 - name: Institute of Earth Sciences, NAWI Graz Geocenter Austria, University of Graz, Austria
   index: 1
date: 2021
bibliography: paper.bib
---

# Summary

The evaporation of water from land surfaces to the atmosphere is a major component of the hydrological cycle. Accurate estimates of this flux are essential to the water and agricultural sectors. As direct measurement of evaporation is difficult, the evaporation flux is commonly estimated from more easily obtained meteorological data using empirical formulas. There is a wide variety of such formulas, ranging from complex physically based methods (e.g., Penman-Monteith), to simpler temperature-based methods (e.g., Hamon). Most of the formulas estimate the potential evaporation (PE), the maximum amount of water that can be evaporated if enough water and energy are available. Water resources experts rely on these estimates to get a first impression of the hydrological water balance. Due to the diversity in the level of input data, process representations, and assumptions, the technical implementation of these methods can be challenging and prone to errors. The goal of `pyet` is to provide a Python package that includes a wide range of methods for the estimation of PE, is fully documented, tested, and easy to use. Currently, `pyet` includes eighteen different methods to estimate daily PE, and various methods to estimate surface and aerodynamic resistance and other missing meteorological data. 

# Statement of Need

A number of standardized calculation methods of potential and reference evaporation exist, from which the benchmark reports of the American Society of Civil Engineers [@walter2000asce] and the Food and Agriculture Organization [@allen1998crop] are among the most widely adopted. There are several Python packages that implemented one of the estimation methods (@pyet0_2019, @eto_2019, @refet_2020 and @evap_2020). None of these packages provide implementations of alternative methods with a common API, and often lack testing and documentation. Outside of Python, the R package `Evapotranspiration` (@Danlu2016r) is the most closely related to `pyet` project and provides similar functionalities to the R community. The aim of `pyet` is to provide a simple-to-use and well documented Python package that includes a wide variety of evaporation formulas for the Python community. The implemented methods in `pyet` are benchmarked against literature values and tested with continuous integration to ensure the correctness of the implementation. Allowing multiple levels of input data, `pyet` is also applicable in regions with sparsely distributed measurement stations, where standard meteorological data (e.g., wind, relative humidity) are often unavailable.

# Estimation of Evaporation 

Evaporation is the process where a substance (here, water) is converted from a liquid into a vapor phase. In this paper, the term evaporation is used to refer to the total evaporation flux from a land surface, consisting of transpiration (evaporation of water by vegetation), soil evaporation, and interception evaporation [@miralles2020]. The potential evaporation flux is determined by meteorological conditions, whereas water availability determines if actual evaporation occurs at its potential rate. Measurements of actual evaporation are limited to point measurements and require expensive instrumentation, which is often unavailable for practical applications. A common method of obtaining actual evaporation rates is thus by using hydrological models, which compute actual evaporation by reducing PE due to water limitations. These models require accurate input of PE, which is often calibrated to match regional and plot characteristics [@OUDIN2005290].

A widely used formula to determine the potential evaporation is using the Penman-Monteith formulation [@allen1998crop, @walter2000asce], which requires measurements of air temperature at 2 m above the ground, relative humidity, wind speed, and radiation. These requirements limit the applicability of the Penman-Monteith, as not all input variables are as widely available. Alternative methods that require less input data are therefore applied in more data scarce regions or in climate scenarios, which often only offer projections of radiation and temperature. For an introduction to the importance of alternative methods to estimate PE, we refer to @OUDIN2005290.

# Functionality

`Pyet` currently contains 18 different PE methods (summarized in Table 1). The package also provides utility functions to estimate missing meteorological data (e.g., solar and net radiation) and the methods are currently implemented for 1D data (time series data). In future work we hope to expand the functionality to 2D and 3D data (Numpy array, XArray, and NetCDF). The implemented methods are benchmarked against open-source data to ensure proper functioning of the methods [@allen1998crop] using continous integration.

Table 1: PE and surface/aerodynamic resistance methods included in `pyet`. T, Temperature; U, Wind Speed; 
D, Radiation; RH, Relative Humidity; $h_{crop}$, crop height; LAI, Leaf area index; $[CO_2]$ - atmospheric
$CO_2$ concentration. Adapted from @OUDIN2005290.

| Method            | Data needed     | PyEt Method       | Reference                   |
|-------------------|-----------------|-------------------|-----------------------------|
| Penman            | RH, T, U, D     |`penman`           |@penman1948natural           |
| Penman-Monteith   | RH, T, U, D     |`pm`               |@monteith1965evaporation     |
|                   |                 |                   |@schymanski_2017             |
| FAO-56            | RH, T, U, D     |`pm_fao56`         |@allen1998crop               |
| Priestley-Taylor  | T, D            |`priestley_taylor` |@priestley1972assessment     |
| Kimberly-Penman   | RH, T, U, D     |`kimberly_penman`  |@wright1982new               |
| Thom-Oliver       | RH, T, U, D     |`thom_oliver`      |@thom1977penman              |
| Blaney–Criddle    | T, D            |`blaney_criddle`   |@blaney1952determining       |
| Hamon             | T               |`hamon`            |@hamon1963estimating         |
| Romanenko         | RH, T           |`romanenko`        |@xu2001evaluation            |
| Linacre           | T               |`linacre`          |@linacre1977simple           |
| Turc              | T, D            |`turc`             |@xu2001evaluation            |
| Jensen–Haise      | T, D            |`jensen_haise`     |@jensen1963estimating        |
| McGuinness–Bordne | T, D            |`mcguinness_bordne`|McGuinness & Bordne          |
|                   |                 |                   |(1972)                       |
| Hargreaves        | T               |`hargreaves`       |Hargreaves & Samani          |
|                   |                 |                   |(1982)                       |
| Doorenbos–Pruitt  | RH, T, U, D     |`fao_24`           |@jensen1990evapotranspiration|
|(FAO-24)           |                 |                   |                             |
| Abtew             | T, D            |`abtew`            |@abtew1996evapotranspiration |
| Makkink           | T, D            |`makkink`          |@makkink1957testing          |
| Oudin             | T               |`oudin`            |@OUDIN2005290                |
| Aerodynamic       | U,($h_{crop}$)  |`calc_res_aero`    |@allen1998crop               |
| resistance        |                 |                   |                             |
| Surface           | (LAI),($[CO_2]$)|`calc_res_surf`    |@allen1998crop               |
| resistance        |                 |                   |@Yang2018HydrologicIO        |

# Example application

This example shows how `pyet` can be used to estimate potential evaporation using different methods. The potential evaporation is estimated for the town of De Bilt (The Netherlands), using data provided by the The Royal Netherlands Meteorological Institute (KNMI). The potential evaporation is computed using three different methods (Penman, Priestley-Taylor, and Makkink) and plotted for comparison.

``` python
# Import needed packages
import pandas as pd
import pyet

# Import meteorological data 
data = pd.read_csv("data/etmgeg_260.txt", skiprows=46, delimiter=",", 
                   skipinitialspace=True, index_col="YYYYMMDD", 
				   parse_dates=True)
tmean, tmax, tmin, rh, wind, rs = (data.TG/10, data.TX/10, data.TN/10,
                                   data.UG, data.FG/10, data.Q/100)

# Compute Evaporation
et_penman = pyet.pm(tmean, wind, rs=rs, elevation=4, lat=0.91, 
					tmax=tmax, tmin=tmin, rh=rh)
et_pt = pyet.priestley_taylor(tmean, wind, rs=rs, elevation=4, 
							  lat=0.91, tmax=tmax, tmin=tmin, rh=rh)
et_makkink = pyet.makkink(tmean, rs, elevation=4)
```
In the code above `rs` is the incoming solar radiation [MJ m-2 d-1], `elevation` the altitude above sea level [m], `lat` the site latitude [rad], `tmean` the mean temperature, `tmax` and `tmin` the maximum and minimum temperature [°C] and `rh` the mean relative humidity [%]. The results show that PE differs based on the estimation method that is used, thus demonstrating the uncertainty associated with the PE estimation methods.

![Daily and cumulative potential evaporation for De Bilt estimated according to @monteith1965evaporation, 
@priestley1972assessment, @makkink1957testing and @OUDIN2005290.](Figure1.png)


# Concluding remarks

In this paper we presented `pyet`, a Python package to estimate daily potential evaporation using available meteorological data. With `pyet`, users can estimate PE using 18 different methods with only a few lines of Python code. At this stage (pyet v1.0), the methods are implemented for 1D data (e.g. time series data).
Further developments will focus on enabling 2D and 3D data input (Numpy Arrays, Array, and NetCDF files). The `pyet` design allows a simple extension of the software with new capabilities and new model options. The authors welcome code contributions, bug reports, and feedback from the community to further improve the package. `pyet` is free and open-source software under the MIT license and is available at http://www.github.com/phydrus/pyet. Full documentation is available on ReadTheDocs (https://pyet.readthedocs.io). The notebook and data necessary to reproduce the figures in this manuscript are available through `pyet` GitHub page.

# Acknowledgements

The first author is partly funded by the Earth System Sciences programme of the Austrian Academy of 
Sciences (project ClimGrassHydro). The second author is funded by the Austrian Science Fund (FWF) under Research 
Grant W1256 (Doctoral Programme Climate 525 Change: Uncertainties, Thresholds and Coping Strategies).

# References
