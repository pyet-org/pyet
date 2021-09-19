<img src=https://raw.githubusercontent.com/phydrus/pyet/d7fdd87719588c00326e692f3b1a47b32161e533/docs/_static/logo.png width=120, align=left>

# pyet: Estimation of Potential Evaporation

<a href="http://www.gnu.org/licenses/gpl-3.0.txt"><img src=https://img.shields.io/github/license/phydrus/pyet> </a>
[![Build Status](https://travis-ci.org/phydrus/PyEt.svg?branch=master)](https://travis-ci.org/github/phydrus/PyEt)
<a href="https://pypi.python.org/pypi/pyet"> <img src=https://img.shields.io/pypi/v/pyet.svg> </a>
[![Documentation Status](https://readthedocs.org/projects/pyet/badge/?version=latest)](https://pyet.readthedocs.io/en/latest/?badge=latest)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/7ed73a2e80784ccf90317c1af8c0cc17)](https://app.codacy.com/gh/phydrus/pyet?utm_source=github.com&utm_medium=referral&utm_content=phydrus/pyet&utm_campaign=Badge_Grade_Dashboard)
[![Codacy Badge](https://app.codacy.com/project/badge/Coverage/e49f23e356f441688422ec32cfcf6aaa)](https://www.codacy.com/gh/phydrus/pyet/dashboard?utm_source=github.com&utm_medium=referral&utm_content=phydrus/pyet&utm_campaign=Badge_Coverage)

pyet is an open source python package for calculating reference and potential 
evaporation (PE). Currently eighteen methods for calculating PE are 
implemented:

| Classification | Common method name        | Data needed | pyet Method        | Reference                   |
|----------------|---------------------------|-------------|--------------------|-----------------------------|
| Combination    | Penman                    | RH, T, U, D |`penman`            |Penman (1948)                |
|                | Penman-Monteith           | RH, T, U, D |`pm`                |Monteith (1965)              |
|                | FAO-56                    | RH, T, U, D |`pm_fao56`          |Allen et al. (1998)          |
|                | Priestley-Taylor          | T, D        |`priestley_taylor`  |Priestley and Taylor (1972)  |
|                | Kimberly-Penman           | RH, T, U, D |`kimberly_penman`   |Wright (1982)                |
|                | Thom-Oliver               | RH, T, U, D |`thom_oliver`       |Thom and Oliver (1977)       |
| Temperature    | Blaney–Criddle            | T, D        |`blaney_criddle`    |Blaney and Criddle (1952)    |
|                | Hamon                     | T           |`hamon`             |Hamon (1963)                 |
|                | Romanenko                 | RH, T       |`romanenko`         |Xu and Singh (2001)          |
|                | Linacre                   | T           |`linacre`           |Linacre (1977)               |
| Radiation      | Turc                      | T, D        |`turc`              |Xu and Singh (2001)          |
|                | Jensen–Haise              | T, D        |`jensen_haise`      |Jensen (1963)                |
|                | McGuinness–Bordne         | T, D        |`mcguinness_bordne` |McGuinness (1972)            |
|                | Hargreaves                | T           |`hargreaves`        |Hargreaves and Samani (1982) |
|                | Doorenbos–Pruitt (FAO-24) | RH, T, U, D |`fao_24`            |Jensen et al. (1990)         |
|                | Abtew                     | T, D        |`abtew`             |Abtew (1996)                 |
|                | Makkink                   | T, D        |`makkink`           |Makkink (1957)               |
|                | Oudin                     | T           |`oudin`             |Oudin (2005)                 |

T, Temperature; U, Wind Speed; D, Radiation; RH, Relative Humidity. Adapted from [@oudin2005potential].

## Examples and Documentation

Examples of using pyet can be found in the example folder. This folder also 
contains a number of Jupyter Notebooks that thoroughly explain the use of the 
software. Documentation is hosted on [ReadTheDocs](https://pyet.readthedocs.io).

After defining the input data, evaporation is estimated using only one 
line of python code:

`>>> pyet.pm_fao56(tmean, wind, rn=rn, tmax=tmax, tmin=tmin, rh=rh, elevation=elevation)`

## Bug reports and Questions

pyet is in active development, and bug reports are welcome as [GitHub 
Issues](https://github.com/phydrus/pyet/issues).
General questions or discussions are possible through 
[GitHub Discussions](https://github.com/phydrus/pyet/discussions).

## Installation
The pyet package is available from the Pypi package index and can be installed 
as follows:

`>>> pip install pyet`

To install in developer mode, use the following syntax:

`>>> pip install -e .`

## Citing
To be added...
