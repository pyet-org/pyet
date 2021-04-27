import unittest

import numpy as np
import pandas as pd

import pyet as et


class Testhargreaves(unittest.TestCase):
    def test_hargreaves(self):
        # Based on example 18, p 78 TestFAO56
        tmax = pd.Series([26.6], index=pd.DatetimeIndex(["2015-07-15"]))
        tmin = pd.Series([14.8], index=pd.DatetimeIndex(["2015-07-15"]))
        tmean = (tmax + tmin) / 2
        lat = 45.72 * np.pi / 180
        har = et.hargreaves(tmean, tmax, tmin, lat)
        self.assertAlmostEqual(float(har), 5.0, 1)
