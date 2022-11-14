"""This file tests all methods for a minimal functioning."""
import unittest

import numpy as np
import pandas as pd

import pyet as et

tmean = pd.Series(data=20 * np.sin(np.linspace(0, 1, 365) * 2 * np.pi),
                  index=pd.date_range("2001-01-01", "2001-12-31", freq="D"))
tmax = tmean + 2
tmin = tmean - 2
rh = pd.Series(data=60 * np.sin(np.linspace(0, 1, 365) * np.pi),
               index=pd.date_range("2001-01-01", "2001-12-31", freq="D"))
rs = pd.Series(data=10 * np.sin(np.linspace(0, 1, 365) * np.pi),
               index=pd.date_range("2001-01-01", "2001-12-31", freq="D"))
wind = pd.Series(data=5 * np.sin(np.linspace(0, 1, 365) * np.pi),
                 index=pd.date_range("2001-01-01", "2001-12-31", freq="D"))
lat = 0.9
elevation = 20


class Testall(unittest.TestCase):
    def test_calculate_all(self):
        et_df = et.calculate_all(tmean, wind, rs, elevation, lat, tmax, tmin,
                                 rh)
        self.assertIsInstance(et_df, pd.DataFrame)
