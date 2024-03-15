<img src=https://raw.githubusercontent.com/phydrus/pyet/d7fdd87719588c00326e692f3b1a47b32161e533/docs/_static/logo.png width=120, align=left>

# pyet: Estimation of Potential Evapotranspiration

[![codacy-coverage-reporter](https://github.com/pyet-org/pyet/actions/workflows/ci.yml/badge.svg)](https://github.com/pyet-org/pyet/actions/workflows/ci.yml)
<a href="https://pypi.org/project/pyet/"><img src=https://img.shields.io/pypi/v/pyet.svg><a>
<a href="https://mit-license.org/"><img src=https://img.shields.io/pypi/l/pyet.svg><a>
[![Codacy Badge](https://app.codacy.com/project/badge/Grade/e49f23e356f441688422ec32cfcf6aaa)](https://www.codacy.com/gh/phydrus/pyet/dashboard?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=phydrus/pyet&amp;utm_campaign=Badge_Grade)
[![Codacy Badge](https://app.codacy.com/project/badge/Coverage/e49f23e356f441688422ec32cfcf6aaa)](https://www.codacy.com/gh/phydrus/pyet/dashboard?utm_source=github.com&utm_medium=referral&utm_content=phydrus/pyet&utm_campaign=Badge_Coverage)
<a href="https://pyet.readthedocs.io/en/latest/?badge=latest"><img src="https://readthedocs.org/projects/pyet/badge/?version=latest"><a>
<a href="https://doi.org/10.5281/zenodo.5896800"><img src=https://zenodo.org/badge/DOI/10.5281/zenodo.5896800.svg><a>

pyet is an open source python package for calculating reference and potential Evapotranspiration (PET) for 1D (pandas.Series)
and 3D (xarray.DataArrray) data. Currently, eighteen methods for calculating daily PET are implemented:

| Method name       | pyet function     | T      | RH         | R      | u2     | Lat.   | El.    | Benchmarked?      |
|:------------------|:------------------|:-------|:-----------|:-------|:-------|:-------|:-------|:------------------|
| Penman            | penman            | &check; $^a$ | &check; $^{b,c}$ | &check; $^d$ | &check;      | &check; $^d$ | &check; $^e$ | &check;           |
| Penman-Monteith   | pm                | &check; $^a$ | &check; $^{b,c}$ | &check; $^d$ | &check;      | &check; $^d$ | &check; $^e$ | &check;           |
| ASCE-PM           | pm_asce           | &check; $^a$ | &check; $^{b,c}$ | &check; $^d$ | &check;      | &check; $^d$ | &check; $^e$ | &check;           |
| FAO-56            | pm_fao56          | &check; $^a$ | &check; $^{b,c}$ | &check; $^d$ | &check;      | &check; $^d$ | &check; $^e$ | &check;           |
| Priestley-Taylor  | priestley_taylor  | &check;      | &check; $^h$     | &check; $^h$ | -      | &check; $^h$ | &check; $^e$ | &check;           |
| Kimberly-Penman   | kimberly_penman   | &check; $^a$ | &check; $^{b,c}$ | &check; $^d$ | &check;      | &check; $^d$ | &check; $^e$ | -                 |
| Thom-Oliver       | thom_oliver       | &check; $^a$ | &check; $^{b,c}$ | &check; $^d$ | &check;      | &check; $^d$ | &check; $^e$ | -                 |
| Blaney-Criddle    | blaney_criddle    | &check;      | - $^i$     | - $^i$ | - $^i$ | &check;      | -      | &check;           |
| Hamon             | hamon             | &check;      | -          | -      | -      | &check;      | -      | &check;           |
| Romanenko         | romanenko         | &check;      | &check;          | -      | -      | -      | -      | &check;           |
| Linacre           | linacre           | &check; $^j$ | -          | -      | -      | -      | &check;      | &check;           |
| Haude             | haude             | &check;      | &check; $^k$     | -      | -      | -      | -      | &check;           |
| Turc              | turc              | &check;      | &check;          | &check;      | -      | -      | -      | &check;           |
| Jensen-Haise      | jensen_haise      | &check;      | -          | &check; $^l$ | -      | &check; $^l$ | -      | &check;           |
| McGuinness-Bordne | mcguinness_bordne | &check;      | -          | -      | -      | &check;      | -      | &check;           |
| Hargreaves        | hargreaves        | &check; $^m$ | -          | -      | -      | &check;      | -      | &check;           |
| FAO-24 radiation  | fao_24            | &check;      | &check;          | &check;      | &check;      | -      | &check; $^e$ | -                 |
| Abtew             | abtew             | &check;      | -          | &check;      | -      | -      | -      | &check;         |
| Makkink           | makkink           | &check;      | -          | &check;      | -      | -      | &check; $^e$ | &check;           |
| Oudin             | oudin             | &check;      | -          | -      | -      | &check;      | -      | -                 |

$^a$ $T_{max}$ and $T_{min}$ can also be provided. $^b$ $RH_{max}$ and $RH_{min}$ can also be provided. $^c$ If actual vapor pressure is provided, RH is not needed.  $^d$ Input for radiation can be (1) Net radiation, (2) solar radiation or (3) sunshine hours. If (1), then latitude is not needed. If (1, 3) latitude and elevation is needed. $^e$ One must provide either the atmospheric pressure or elevation. $^f$ The PM method can be used to estimate potential crop evapotranspiration, if leaf area index or crop height data is available. $^g$ The effect of $CO_2$ on stomatal resistance can be included using the formulation of Yang et al. 2019.  $^h$ If net radiation is provided, RH and Lat are not needed. $^i$ If method==2, $u_2$, $RH_{min}$ and sunshine hours are required. $^j$ Additional input of $T_{max}$ and $T_{min}$, or $T_{dew}$. $^k$ Input can be $RH$ or actual vapor pressure. $^l$ If method==1, latitude is needed instead of $R_s$. $^m$ $T_{max}$ and $T_{min}$ also needed.

## Examples and Documentation

Examples of using *pyet* can be found in the example folder:

*   [Example 1](<https://pyet.readthedocs.io/en/dev/examples/01_example_zamg.html#>): Estimating PET using pandas.Series

*   [Example 2](<https://pyet.readthedocs.io/en/dev/examples/02_example_zamg_netcdf.html>): Estimating PET using xarray.DataArray

*   [Example 3](<https://pyet.readthedocs.io/en/dev/examples/03_example_knmi.html>): Benchmarking Makkink
  against [KNMI data](https://www.knmi.nl/over-het-knmi/about)

*   [Example 4](<https://pyet.readthedocs.io/en/dev/examples/04_example_coagmet.html>): Benchmarking FAO56
  against [CoAgMET data](https://coagmet.colostate.edu/)

*   [Example 5](<https://pyet.readthedocs.io/en/dev/examples/05_example_calibration.html>): Calibrating the Romanenko and Abtew method against the PM-FAO56

*   [Example 6](<https://pyet.readthedocs.io/en/dev/examples/06_worked_examples_McMahon_etal_2013.html>): Worked examples for estimating meteorological
  variables and potential evapotranspiration after McMahon et al., 2013

*   [Example 7](<https://pyet.readthedocs.io/en/dev/examples/07_example_climate_change.html>): Example for estimating potential evapotranspiration under
  warming and elevated $CO_2$ concentrations following Yang et al., (2019)

*   [Example 8](<https://pyet.readthedocs.io/en/dev/examples/08_crop_coefficient.html>): Determining the crop coefficient function with Python

*   [Example 9](<https://pyet.readthedocs.io/en/dev/examples/09_CMIP6_data.html>): Estimating PET using CMIP data

*   [Example 10](<https://pyet.readthedocs.io/en/dev/examples/10_example_paper.html>): Notebook supporting PyEt manuscript

Documentation is hosted on [ReadTheDocs](https://pyet.readthedocs.io).

After defining the input data, evapotranspiration is estimated using only one
line of python code:

`>>> pyet.pm_fao56(tmean, wind, rn=rn, tmax=tmax, tmin=tmin, rh=rh, elevation=elevation)`

We support Python >= 3.8.

## Benchmarking

Most of the methods implemented in *pyet* are benchmarked against literature values from the [FAO Irrigation and
drainage paper 56](https://www.fao.org/3/x0490e/x0490e00.htm), [McMahon et al., 2013 (supplementary)](https://hess.copernicus.org/articles/17/4865/2013/) and [Schrödter, 1985](https://link.springer.com/book/10.1007/978-3-642-70434-5). In addition, two comparative analysis between daily PE estimated with *pyet* and other organizations is
made:

*   `pyet.pm_fao56` against daily PET estimated with ASCE Penman-Monteith from [CoAgMET](https://coagmet.colostate.edu/) (
  Colorado State University),

*   `pyet.makkink` against daily PET estimated with Makkink from The Royal Netherlands Meteorological
  Institute ([KNMI](https://www.knmi.nl/over-het-knmi/about)).

## Data dimensions

As of version v1.2., *pyet* is compatible with both Pandas.Series and xarray.DataArray, which means you can now estimate
potential evapotranspiration for both point and gridded data.

## Bug reports and Questions

pyet is in active development, and bug reports are welcome as [GitHub
Issues](https://github.com/phydrus/pyet/issues).
General questions or discussions are possible through
[GitHub Discussions](https://github.com/phydrus/pyet/discussions).

## Installation

The *pyet* package is available from the Pypi package index and can be installed
as follows::

`>>> pip install pyet`

To install in developer mode, use the following syntax:

`>>> pip install -e .`

## Citing

If you use *pyet* in one of your studies, please cite the *pyet* EGU abstract:

* Vremec, M., Collenteur, R. A., and Birk, S.: Technical note: Improved handling of potential evapotranspiration in
hydrological studies with PyEt, Hydrol. Earth Syst. Sci. Discuss. [preprint], https://doi.org/10.5194/hess-2022-417,
in review, 2023.

```Reference
@Article{hess-2022-417,
  AUTHOR = {Vremec, M. and Collenteur, R. A. and Birk, S.},
  TITLE = {Technical note: Improved handling of potential evapotranspiration in hydrological studies with \textit{PyEt}},
  JOURNAL = {Hydrology and Earth System Sciences Discussions},
  VOLUME = {2023},
  YEAR = {2023},
  PAGES = {1--23},
  URL = {https://hess.copernicus.org/preprints/hess-2022-417/},
  DOI = {10.5194/hess-2022-417}
}
```
