<img src=https://raw.githubusercontent.com/phydrus/pyet/d7fdd87719588c00326e692f3b1a47b32161e533/docs/_static/logo.png width=120, align=left>

# pyet: Estimation of Potential Evaporation

<a href="https://travis-ci.org/github/phydrus/PyEt"><img src="https://api.travis-ci.org/phydrus/PyEt.svg?branch=master"><a>
<a href="https://mit-license.org/"><img src=https://img.shields.io/pypi/v/pyet.svg><a>
<a href="https://travis-ci.org/github/phydrus/PyEt"><img src=https://img.shields.io/pypi/l/pyet.svg><a>
[![Codacy Badge](https://app.codacy.com/project/badge/Grade/e49f23e356f441688422ec32cfcf6aaa)](https://www.codacy.com/gh/phydrus/pyet/dashboard?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=phydrus/pyet&amp;utm_campaign=Badge_Grade)
[![Codacy Badge](https://app.codacy.com/project/badge/Coverage/e49f23e356f441688422ec32cfcf6aaa)](https://www.codacy.com/gh/phydrus/pyet/dashboard?utm_source=github.com&utm_medium=referral&utm_content=phydrus/pyet&utm_campaign=Badge_Coverage)
<a href="https://pyet.readthedocs.io/en/latest/?badge=latest"><img src="https://readthedocs.org/projects/pyet/badge/?version=latest"><a>   
<a href="https://doi.org/10.5281/zenodo.5896800"><img src=https://zenodo.org/badge/DOI/10.5281/zenodo.5896800.svg><a>

pyet is an open source python package for calculating reference and potential evaporation (PE) for 1D (pandas.Series)
and 3D (xarray.DataArrray) data. Currently, eighteen methods for calculating daily PE are implemented:

| Classification | Common method name        | Data needed | pyet Method         | Reference                    |
|----------------|---------------------------|-------------|---------------------|------------------------------|
| Combination    | Penman                    | RH, T, U, R | `penman`            | Penman (1948)                |
|                | Penman-Monteith           | RH, T, U, R | `pm`                | Monteith (1965)              |
|                | Penman-Monteith ASCE      | RH, T, U, R | `pm`                | ASCE (2005)                  |
|                | FAO-56                    | RH, T, U, R | `pm_fao56`          | Allen et al. (1998)          |
|                | Priestley-Taylor          | RH, T, R    | `priestley_taylor`  | Priestley and Taylor (1972)  |
|                | Kimberly-Penman           | RH, T, U, R | `kimberly_penman`   | Wright (1982)                |
|                | Thom-Oliver               | RH, T, U, R | `thom_oliver`       | Thom and Oliver (1977)       |
| Temperature    | Blaney–Criddle            | T           | `blaney_criddle`    | Blaney and Criddle (1952)    |
|                | Hamon                     | T           | `hamon`             | Hamon (1963), 3 options      |
|                | Romanenko                 | RH, T       | `romanenko`         | Xu and Singh (2001)          |
|                | Linacre                   | T           | `linacre`           | Linacre (1977)               |
|                | Haude                     | RH, T       | `haude`             | Haude (1955)                 |
| Radiation      | Turc                      | RH, T, R    | `turc`              | Xu and Singh (2001)          |
|                | Jensen–Haise              | T, (R)      | `jensen_haise`      | Jensen (1963)                |
|                | McGuinness–Bordne         | T           | `mcguinness_bordne` | McGuinness (1972)            |
|                | Hargreaves                | T           | `hargreaves`        | Hargreaves and Samani (1982) |
|                | Doorenbos–Pruitt (FAO-24) | RH, T, U, R | `fao_24`            | Jensen et al. (1990)         |
|                | Abtew                     | T, R        | `abtew`             | Abtew (1996)                 |
|                | Makkink                   | T, R        | `makkink`           | Makkink (1957)               |
|                | Oudin                     | T           | `oudin`             | Oudin (2005)                 |

T, Temperature; U, Wind Speed; D, Radiation; RH, Relative Humidity. Adapted from Oudin et al. (2005).

## Examples and Documentation

Examples of using pyet can be found in the example folder:

*   [Example 1](/examples/01_example_zamg.ipynb): Estimating PE using pandas.Series

*   [Example 2](/examples/02_example_zamg_netcdf.ipynb): Estimating PE using xarray.DataArray

*   [Example 3](/examples/03_example_knmi.ipynb): Benchmarking Makkink
  against [KNMI data](https://www.knmi.nl/over-het-knmi/about)

*   [Example 4](/examples/04_example_coagmet.ipynb): Benchmarking FAO56
  against [CoAgMET data](https://coagmet.colostate.edu/)

*   [Example 5](/examples/05_example_calibration.ipynb): Calibrating the Romanenko and Abtew method against the PM-FAO56

*   [Example 6](/examples/06_worked_examples_McMahon_etal_2013.ipynb): Worked examples for estimating meteorological
  variables and potential evaporation after McMahon et al., 2013

Documentation is hosted on [ReadTheDocs](https://pyet.readthedocs.io).

After defining the input data, evaporation is estimated using only one
line of python code:

`>>> pyet.pm_fao56(tmean, wind, rn=rn, tmax=tmax, tmin=tmin, rh=rh, elevation=elevation)`

We support Python >= 3.8.

## Benchmarking

Most of the methods implemented in pyet are benchmarked against literature values from the [FAO Irrigation and
drainage paper 56](https://www.fao.org/3/x0490e/x0490e00.htm), [McMahon et al., 2013 (supplementary)](https://hess.copernicus.org/articles/17/4865/2013/) and [Schrödter, 1985](https://link.springer.com/book/10.1007/978-3-642-70434-5). In addition, two comparative analysis between daily PE estimated with pyet and other organizations is
made:

*   `pyet.pm_fao56` against daily PE estimated with ASCE Penman-Monteith from [CoAgMET](https://coagmet.colostate.edu/) (
  Colorado State University),

*   `pyet.makkink` against daily PE estimated with Makkink from The Royal Netherlands Meteorological
  Institute ([KNMI](https://www.knmi.nl/over-het-knmi/about)).

## Dimensions

As of version v1.2.0, pyet is compatible with both Pandas.Series and xarray.DataArray, which means you can now estimate
potential evaporation for both point and gridded data.

## Bug reports and Questions

pyet is in active development, and bug reports are welcome as [GitHub
Issues](https://github.com/phydrus/pyet/issues).
General questions or discussions are possible through
[GitHub Discussions](https://github.com/phydrus/pyet/discussions).

## Installation

The pyet package is available from the Pypi package index and can be installed
as follows::

`>>> pip install pyet`

To install in developer mode, use the following syntax:

`>>> pip install -e .`

## Citing

If you use PyEt in one of your studies, please cite the PyEt EGU abstract:

*   Vremec, M. and Collenteur, R.: PyEt - a Python package to estimate potential and reference evapotranspiration, EGU
  General Assembly 2021, online, 19–30 Apr 2021, EGU21-15008, https://doi.org/10.5194/egusphere-egu21-15008, 2021.

```Reference
@inproceedings{vremec2021pyet,
  title={PyEt-a Python package to estimate potential and reference evapotranspiration},
  author={Vremec, Matev{\v{z}} and Collenteur, Raoul},
  booktitle={EGU General Assembly Conference Abstracts},
  pages={EGU21--15008},
  year={2021},
  doi={https://doi.org/10.5194/egusphere-egu21-15008}
}
```
