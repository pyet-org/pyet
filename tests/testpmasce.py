import unittest
import numpy as np
import pandas as pd

import pyet as et


class TestASCE(unittest.TestCase):
    def test_et_asce(self):
        # Based on Example 18, p. 74 FAO.
        dindex = pd.DatetimeIndex(["2015-07-06"])
        wind = pd.Series([2.078], index=dindex)
        tmax = pd.Series([21.5], index=dindex)
        tmin = pd.Series([12.3], index=dindex)
        rhmax = pd.Series([84], index=dindex)
        rhmin = pd.Series([63], index=dindex)
        rs = pd.Series([22.07], index=dindex)
        n = 9.25
        nn = 16.1
        elevation = 100
        latitude = 50.80 * np.pi / 180
        pm_asce = round(et.pm_asce(wind, elevation, latitude, solar=rs,
                                   tmax=tmax, tmin=tmin, rhmax=rhmax,
                                   rhmin=rhmin, n=n, nn=nn), 1)
        self.assertAlmostEqual(pm_asce.values, 3.9, 1)
