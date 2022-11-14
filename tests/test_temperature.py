"""This file tests for the temperature methods for a minimal functioning."""

import numpy as np
import pandas as pd

import pyet as et

tmean = pd.Series(data=20 * np.sin(np.linspace(0, 1, 365) * 2 * np.pi),
                  index=pd.date_range("2001-01-01", "2001-12-31", freq="D"))


def test_blaney_criddle():
    et.blaney_criddle(tmean, lat=0.9)
    return True


def test_haude():
    et.haude(tmean, rh=0.5)
    return True


def test_hamon():
    et.hamon(tmean, lat=0.9)
    return True


def test_romanenko():
    et.romanenko(tmean, rh=0.5)
    return True


def test_linacre():
    et.linacre(tmean, 0, lat=0.9)
    return True
