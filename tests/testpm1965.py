import unittest
import numpy as np
import pandas as pd

import pyet as et


class TestFAO56(unittest.TestCase):
    def test_et_fao56(self):
        # Based on Example 18, p. 74 FAO.
        wind = pd.Series([2.078], index=pd.DatetimeIndex(["2015-07-06"]))
        tmax = pd.Series([21.5], index=pd.DatetimeIndex(["2015-07-06"]))
        tmin = pd.Series([12.3], index=pd.DatetimeIndex(["2015-07-06"]))
        rhmax = pd.Series([84], index=pd.DatetimeIndex(["2015-07-06"]))
        rhmin = pd.Series([63], index=pd.DatetimeIndex(["2015-07-06"]))
        rs = pd.Series([22.07], index=pd.DatetimeIndex(["2015-07-06"]))
        n = 9.25
        nn = 16.1
        elevation = 100
        latitude = 50.80 * np.pi / 180
        pm_1965 = round(et.pm1965(wind, elevation, latitude, solar=rs,
                                  tmax=tmax, tmin=tmin, rhmax=rhmax,
                                  rhmin=rhmin, n=n, nn=nn), 1)
        self.assertAlmostEqual(pm_1965.values, 3.9, 1)
