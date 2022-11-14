<img src=https://raw.githubusercontent.com/phydrus/pyet/d7fdd87719588c00326e692f3b1a47b32161e533/docs/_static/logo.png width=120, align=left>

# pyet: Estimation of Potential Evapotranspiration

<a href="https://travis-ci.org/github/phydrus/PyEt"><img src="https://api.travis-ci.org/phydrus/PyEt.svg?branch=master"><a>
<a href="https://mit-license.org/"><img src=https://img.shields.io/pypi/v/pyet.svg><a>
<a href="https://travis-ci.org/github/phydrus/PyEt"><img src=https://img.shields.io/pypi/l/pyet.svg><a>
[![Codacy Badge](https://app.codacy.com/project/badge/Grade/e49f23e356f441688422ec32cfcf6aaa)](https://www.codacy.com/gh/phydrus/pyet/dashboard?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=phydrus/pyet&amp;utm_campaign=Badge_Grade)
[![Codacy Badge](https://app.codacy.com/project/badge/Coverage/e49f23e356f441688422ec32cfcf6aaa)](https://www.codacy.com/gh/phydrus/pyet/dashboard?utm_source=github.com&utm_medium=referral&utm_content=phydrus/pyet&utm_campaign=Badge_Coverage)
<a href="https://pyet.readthedocs.io/en/latest/?badge=latest"><img src="https://readthedocs.org/projects/pyet/badge/?version=latest"><a>   
<a href="https://doi.org/10.5281/zenodo.5896800"><img src=https://zenodo.org/badge/DOI/10.5281/zenodo.5896800.svg><a>

pyet is an open source python package for calculating reference and potential Evapotranspiration (PET) for 1D (pandas.Series)
and 3D (xarray.DataArrray) data. Currently, eighteen methods for calculating daily PET are implemented:

| Method name       | pyet function     | T      | RH         | R      | u2     | Lat.   | El.    | Benchmarked? |
|:------------------|:------------------|:-------|:-----------|:-------|:-------|:-------|:-------|:---------|
| Penman            | penman            | &check; $^a$ | &check; $^{b,c}$ | &check; $^d$ | &check;      | &check; $^d$ | &check; $^e$ | -        |
| Penman-Monteith   | pm                | &check; $^a$ | &check; $^{b,c}$ | &check; $^d$ | &check;      | &check; $^d$ | &check; $^e$ | &check;        |
| ASCE-PM           | pm_asce           | &check; $^a$ | &check; $^{b,c}$ | &check; $^d$ | &check;      | &check; $^d$ | &check; $^e$ | &check;        |
| FAO-56            | pm_fao56          | &check; $^a$ | &check; $^{b,c}$ | &check; $^d$ | &check;      | &check; $^d$ | &check; $^e$ | &check;        |
| Priestley-Taylor  | priestley_taylor  | &check;      | &check; $^h$     | &check; $^h$ | -      | &check; $^h$ | &check; $^e$ | &check;        |
| Kimberly-Penman   | kimberly_penman   | &check; $^a$ | &check; $^{b,c}$ | &check; $^d$ | &check;      | &check; $^d$ | &check; $^e$ | -        |
| Thom-Oliver       | thom_oliver       | &check; $^a$ | &check; $^{b,c}$ | &check; $^d$ | &check;      | &check; $^d$ | &check; $^e$ | -        |
| Blaney-Criddle    | blaney_criddle    | &check;      | - $^i$     | - $^i$ | - $^i$ | &check;      | -      | &check;        |
| Hamon             | hamon             | &check;      | -          | -      | -      | &check;      | -      | -        |
| Romanenko         | romanenko         | &check;      | &check;          | -      | -      | -      | -      | -        |
| Linacre           | linacre           | &check; $^j$ | -          | -      | -      | -      | &check;      | -        |
| Haude             | haude             | &check;      | &check; $^k$     | -      | -      | -      | -      | &check;        |
| Turc              | turc              | &check;      | &check;          | &check;      | -      | -      | -      | &check;        |
| Jensen-Haise      | jensen_haise      | &check;      | -          | &check; $^l$ | -      | &check; $^l$ | -      | -        |
| McGuinness-Bordne | mcguinness_bordne | &check;      | -          | -      | -      | &check;      | -      | -        |
| Hargreaves        | hargreaves        | &check; $^m$ | -          | -      | -      | &check;      | -      | &check;        |
| FAO-24            | fao_24            | &check;      | &check;          | &check;      | &check;      | -      | &check; $^e$ | -        |
| Abtew             | abtew             | &check;      | -          | &check;      | -      | -      | -      | -        |
| Makkink           | makkink           | &check;      | -          | &check;      | -      | -      | &check; $^e$ | &check;        |
| Oudin             | oudin             | &check;      | -          | -      | -      | &check;      | -      | -        |

$^a$ $T_{max}$ and $T_{min}$ can also be provided. $^b$ $RH_{max}$ and $RH_{min}$ can also be provided. $^c$ If actual vapor pressure is provided, RH is not needed.  $^d$ Input for radiation can be (1) Net radiation, (2) solar radiation or (3) sunshine hours. If (1), then latitude is not needed. If (1, 3) latitude and elevation is needed. $^e$ One must provide either the atmospheric pressure or elevation. $^f$ The PM method can be used to estimate potential crop evapotranspiration, if leaf area index or crop height data is available. $^g$ The effect of $CO_2$ on stomatal resistance can be included using the formulation of Yang et al. 2019.  $^h$ If net radiation is provided, RH and Lat are not needed. $^i$ If method==2, $u_2$, $RH_{min}$ and sunshine hours are required. $^j$ Additional input of $T_{max}$ and $T_{min}$, or $T_{dew}$. $^k$ Input can be $RH$ or actual vapor pressure. $^l$ If method==1, latitude is needed instead of $R_s$. $^m$ $T_{max}$ and $T_{min}$ also needed.
 
## Examples and Documentation

Examples of using *pyet* can be found in the example folder:

*   [Example 1](/examples/01_example_zamg.ipynb): Estimating PET using pandas.Series

*   [Example 2](/examples/02_example_zamg_netcdf.ipynb): Estimating PET using xarray.DataArray

*   [Example 3](/examples/03_example_knmi.ipynb): Benchmarking Makkink
  against [KNMI data](https://www.knmi.nl/over-het-knmi/about)

*   [Example 4](/examples/04_example_coagmet.ipynb): Benchmarking FAO56
  against [CoAgMET data](https://coagmet.colostate.edu/)

*   [Example 5](/examples/05_example_calibration.ipynb): Calibrating the Romanenko and Abtew method against the PM-FAO56

*   [Example 6](/examples/06_worked_examples_McMahon_etal_2013.ipynb): Worked examples for estimating meteorological
  variables and potential evapotranspiration after McMahon et al., 2013

*   [Example 7](/examples/07_example_climate_change.ipynb): Example for estimating potential evapotranspiration under 
  warming and elevated $CO_2$ concentrations following Yang et al., (2019) 

*   [Example 8](/examples/08_crop_coefficient.ipynb): Determining the crop coefficient function with Python 

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

## Dimensions

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

*   Vremec, M. and Collenteur, R.: *pyet* - a Python package to estimate potential and reference evapotranspiration, EGU
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
