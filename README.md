# PyEt

<a href="http://www.gnu.org/licenses/gpl-3.0.txt"><img src=https://img.shields.io/github/license/phydrus/pyet> </a>
[![Build Status](https://travis-ci.org/phydrus/PyEt.svg?branch=master)](https://travis-ci.org/github/phydrus/PyEt)
<a href="https://pypi.python.org/pypi/pyet"> <img src=https://img.shields.io/pypi/v/pyet.svg> </a>
[![Documentation Status](https://readthedocs.org/projects/pyet/badge/?version=latest)](https://pyet.readthedocs.io/en/latest/?badge=latest)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/7ed73a2e80784ccf90317c1af8c0cc17)](https://app.codacy.com/gh/phydrus/PyEt?utm_source=github.com&utm_medium=referral&utm_content=phydrus/PyEt&utm_campaign=Badge_Grade_Dashboard)


PyEt is an open source python package for calculating both reference and potential evapotranspiration(ET).

In the current state, 8 methods for calculating ET are implemented:
* FAO Penman-Monteith (Monteith, 1965; FAO, 1990)
* FAO-56 Penman-Monteith (Monteith, 1965; Allen et al, 1998)
* Penman (1948)
* Hargreaves (Hargreaves and Samani, 1982; 1985)
* Hamon (1961)
* Priestley and Taylor (1972)
* Jensen and Haise (1963)
* Makkink (1957)

## Examples and Documentation
Examples of using PyEt can be found in the examples folders.

## Installation

To install in developer mode, use the following syntax:

`>>> pip install -e .`
