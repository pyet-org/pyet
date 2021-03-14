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
date: 12 March 2021
bibliography: paper.bib

# Summary

Evaporation (ET) is a major compenent of the hydrological cycle and therefore of high importance for the agricultural and energy sector. 
In this manuscript, the term evaporation is used instead of the commonly used potential evapotranspiration.  

# Introduction

`PyEt` is an python package for calculating reference and potential evaporation (ET). 
PyEt currently includes eighteen methods for calculating ET. 

# Evaporation models

| Classification | Common method name        | Data needed | Reference                      |
|----------------|---------------------------|-------------|--------------------------------|
| Combination    | Penman                    | RH, T, U, D |[@penman1948natural]            |
|                | Penman-Monteith           | RH, T, U, D |[@monteith1965evaporation]      |
|                | Priestley-Taylor          | T, D        |[@priestley1972assessment]      |
|                | Kimberly-Penman           | RH, T, U, D |[@wright1982new]                |
|                | Thom-Oliver               | RH, T, U, D |[@thom1977penman]               |
| Temperature    | Thornthwaite              | T           |[@thornthwaite1948approach]     |
|                | Blaney–Criddle            | T, D        |[@blaney1952determining]        |
|                | Hamon                     | T           |[@hamon1963estimating]          |
|                | Romanenko                 | RH, T       |[@xu2001evaluation]             |
|                | Linacre                   | T           |[@linacre1977simple]            |
| Radiation      | Turc                      | T, D        |[@xu2001evaluation]             |
|                | Jensen–Haise              | T, D        |[@jensen1963estimating]         |
|                | McGuinness–Bordne         | T, D        |[@mcguinness1972comparison]     |
|                | Hargreaves                | T           |[@hargreaves1982estimating]     |
|                | Doorenbos–Pruitt (FAO-24) | RH, T, U, D |[@jensen1990evapotranspiration] |
|                | Abtew                     | T, D        |[@abtew1996evapotranspiration]  |
|                | Makkink                   | T, D        |[@makkink1957testing]           |
|                | Oudin                     | T           |[@oudin2005potential]           |

# Acknowledgements

# References