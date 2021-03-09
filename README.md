# PyEt: Computation of Evapotranspiration

<a href="http://www.gnu.org/licenses/gpl-3.0.txt"><img src=https://img.shields.io/github/license/phydrus/pyet> </a>
[![Build Status](https://travis-ci.org/phydrus/PyEt.svg?branch=master)](https://travis-ci.org/github/phydrus/PyEt)
<a href="https://pypi.python.org/pypi/pyet"> <img src=https://img.shields.io/pypi/v/pyet.svg> </a>
[![Documentation Status](https://readthedocs.org/projects/pyet/badge/?version=latest)](https://pyet.readthedocs.io/en/latest/?badge=latest)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/7ed73a2e80784ccf90317c1af8c0cc17)](https://app.codacy.com/gh/phydrus/PyEt?utm_source=github.com&utm_medium=referral&utm_content=phydrus/PyEt&utm_campaign=Badge_Grade_Dashboard)


PyEt is an open source python package for calculating reference and potential 
evapotranspiration(ET). Currently nine methods for calculating ET are 
implemented:

* FAO Penman-Monteith (Monteith, 1965; FAO, 1990)
* FAO-56 Penman-Monteith (Monteith, 1965; Allen et al, 1998)
* Penman (1948)
* Hargreaves (Hargreaves and Samani, 1982; 1985)
* Hamon (1961)
* Priestley and Taylor (1972)
* Jensen and Haise (1963)
* Makkink (1957)
* Upscaled corrected Penman-Monteith by Schymanski (2017)

## Examples and Documentation

Examples of using PyEt can be found in the example folder. This folder also 
contains a number of Jupyter Notebooks that thoroughly explain the use of the 
software. Documentation is hosted on [ReadTheDocs](pyet.readthedocs.io.).

After defining the input data, evapotranspiration is estimated using only one 
line of python code:

`>>> pyet.pm_fao56(wind, elevation, latitude, solar=solar, tmax=tmax, tmin=tmin, rh=rh)`

## Bug reports and Questions

PyEt is in active development, and bug reports are welcome as [GitHub 
Issues](https://github.com/phydrus/PyEt/issues).
General questions or discussions are possible through 
[GitHub Discussions](https://github.com/phydrus/PyEt/discussions).

## Installation
The PyEt package is available from the Pypi package index and can be installed 
as follows:

`>>> pip install pyet`

To install in developer mode, use the following syntax:

`>>> pip install -e .`

## Citing
To be added...
