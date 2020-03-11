﻿# PyEt
PyEt is an open source python package for calculating both reference and potential evapotranspiration(ET).

In the current state, 10 methods for calculating ET are implemented:
* FAO Penman-Monteith (Monteith, 1965; FAO, 1990)
* FAO-56 Penman-Monteith (Monteith, 1965; Allen et al, 1998)
* Penman (1948)
* Hargreaves (Hargreaves and Samani, 1982; 1985)
* Hamon (1961)
* Priestley and Taylor (1972)
* Kimberly–Penman (Wright, 1982)
* Jensen and Haise (1963)
* Makkink (1957)

## Examples and Documentation
Examples of using PyEt can be found in the examples folders. This folder also containts a Jupyter Notebook that compares PyEt output of Penman-Monteith and Hargreaves method with HYDRUS-1D calculations.

## Installation

To install in developer mode, use the following syntax:

`>>> pip install -e .`
