import unittest

import numpy as np
import pandas as pd

import pyet as et


class Testalternative(unittest.TestCase):
    def test_hargreaves(self):
        # Based on example S19.46,p 78 TestFAO56
        tmax = pd.Series([26.6], index=pd.DatetimeIndex(["2015-07-15"]))
        tmin = pd.Series([14.8], index=pd.DatetimeIndex(["2015-07-15"]))
        tmean = (tmax + tmin) / 2
        lat = 45.72 * np.pi / 180
        har = et.hargreaves(tmean, tmax, tmin, lat)
        self.assertAlmostEqual(float(har), 5.0, 1)

    def test_makkink(self):
        # Based on example S19.91, p 80 McMahon_etal_2013
        tmean = pd.Series([11.5], index=pd.DatetimeIndex(["1980-07-20"]))
        rs = pd.Series([17.194], index=pd.DatetimeIndex(["1980-07-20"]))
        elevation = 546
        lambda0 = 2.45  # Lambda in McMahon
        lambda1 = et.calc_lambda(tmean)  # Lambda in pyet
        lambda_corr = lambda1 / lambda0
        n = -0.12  # correction for different formula in McMahon_2013
        mak = (et.makkink(tmean, rs, elevation=elevation,
                          k=0.61) * lambda_corr) + n
        self.assertAlmostEqual(float(mak), 2.3928, 2)

    def test_blaney_criddle(self):
        # Based on example S19.93, McMahon_etal_2013
        tmean = pd.Series([11.5], index=pd.DatetimeIndex(["1980-07-20"]))
        rhmin = pd.Series([25], index=pd.DatetimeIndex(["1980-07-20"]))
        wind = pd.Series([0.5903], index=pd.DatetimeIndex(["1980-07-20"]))
        n = pd.Series([10.7], index=pd.DatetimeIndex(["1980-07-20"]))
        py = pd.Series([0.2436], index=pd.DatetimeIndex(["1980-07-20"]))
        lat = -23.7951 * np.pi / 180
        bc = et.blaney_criddle(tmean, lat, wind=wind, n=n, rhmin=rhmin,
                               py=py, method=2)
        self.assertAlmostEqual(float(bc), 3.1426, 2)

    def test_turc(self):
        # Based on example S19.99, McMahon_etal_2013
        tmean = pd.Series([11.5], index=pd.DatetimeIndex(["1980-07-20"]))
        rhmean = pd.Series([48], index=pd.DatetimeIndex(["1980-07-20"]))
        rs = pd.Series([17.194], index=pd.DatetimeIndex(["1980-07-20"]))
        turc = et.turc(tmean, rs, rhmean, k=0.32)
        self.assertAlmostEqual(float(turc), 2.6727, 1)

    def test_hargreaves_samani(self):
        # Based on example S19.101, McMahon_etal_2013
        tmean = pd.Series([11.5], index=pd.DatetimeIndex(["1980-07-20"]))
        tmax = pd.Series([21], index=pd.DatetimeIndex(["1980-07-20"]))
        tmin = pd.Series([2], index=pd.DatetimeIndex(["1980-07-20"]))
        lambda0 = 2.45  # Lambda in McMahon
        lambda1 = et.calc_lambda(tmean)  # Lambda in pyet
        lambda_corr = lambda1 / lambda0
        lat = -23.7951 * np.pi / 180
        har = et.hargreaves(tmean, tmax, tmin, lat, method=1)
        self.assertAlmostEqual(float(har * lambda_corr), 4.1129, 2)

    def test_priestley_taylor(self):
        # Based on example S19.109, McMahon_etal_2013
        tmean = pd.Series([11.5], index=pd.DatetimeIndex(["1980-07-20"]))
        rn = pd.Series([8.6401], index=pd.DatetimeIndex(["1980-07-20"]))
        elevation = 546
        lambda0 = 2.45  # Lambda in McMahon
        lambda1 = et.calc_lambda(tmean)  # Lambda in pyet
        lambda_corr = lambda1 / lambda0
        pt = et.priestley_taylor(tmean, rn=rn, elevation=elevation, alpha=1.26)
        self.assertAlmostEqual(float(pt * lambda_corr), 2.6083, 2)

    def test_blaney_criddle1(self):
        # Based on example 5.1, P89 Schrodter 1985
        tmean = pd.Series([17.3], index=pd.DatetimeIndex(["1980-07-20"]))
        lat = 50 * np.pi / 180
        bc = et.blaney_criddle(tmean, lat, method=0)
        self.assertAlmostEqual(float(bc), 3.9, 2)

    def test_haude(self):
        # Based on example 5.2, P95 Schrodter 1985
        tmean = pd.Series([21.5], index=pd.DatetimeIndex(["1980-07-20"]))
        ea = pd.Series([1.19], index=pd.DatetimeIndex(["1980-07-20"]))
        k = 0.26 / 0.35  # 0.35 value is taken in pyet
        e0 = et.calc_e0(tmean)
        rh = ea / e0 * 100
        haude = et.haude(tmean, rh=rh, k=k)
        self.assertAlmostEqual(float(haude), 3.6, 1)
