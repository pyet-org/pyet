import unittest

import numpy as np
import pandas as pd

import pyet as et


class Testr(unittest.TestCase):
    lat = -0.6095
    elevation = 48
    tmax = pd.Series(
        [28.8, 27.4, 29., 26.3, 32.7, 35.7, 33.8, 33.8, 35.6, 27.5],
        index=pd.date_range("2001-03-01", "2001-03-10", freq="d"))
    tmin = pd.Series(
        [15.1, 14., 16.3, 16.2, 17.8, 26.7, 20.6, 22.3, 22.7, 17.1],
        index=pd.date_range("2001-03-01", "2001-03-10", freq="d"))
    tmean = (tmax + tmin) / 2
    rhmax = pd.Series(
        [68, 77, 69, 70, 73, 36, 51, 46, 34, 67],
        index=pd.date_range("2001-03-01", "2001-03-10", freq="d"))
    rhmin = pd.Series(
        [30, 25, 30, 34, 29, 20, 18, 18, 10, 30],
        index=pd.date_range("2001-03-01", "2001-03-10", freq="d"))
    rh = (rhmax + rhmin) / 2
    rs = pd.Series(
        [20.44477761, 20.34991941, 20.25373437, 20.1562439, 20.05747057,
         19.95743816, 19.85617161, 19.75369701, 19.65004162, 19.54523386],
        index=pd.date_range("2001-03-01", "2001-03-10", freq="d"))
    uz = pd.Series(
        [2.65625, 2.78472222, 2.49305556, 3.73611111, 3.30902778,
         4.93402778, 2.90972222, 4.48958333, 3.99305556, 4.31597222],
        index=pd.date_range("2001-03-01", "2001-03-10", freq="d"))
    n = pd.Series(
        [8.6, 8.6, 8.6, 8.6, 8.6, 8.6, 8.6, 8.6, 8.6, 8.6],
        index=pd.date_range("2001-03-01", "2001-03-10", freq="d"))
    lambda0 = 2.45  # Lambda in Danlu R Evapotranspiration package
    lambda1 = et.calc_lambda(tmean)  # Lambda in pyet
    lambda_corr = lambda1 / lambda0

    def test_penman(self):
        wind_penman = Testr.uz * np.log(2 / 0.001) / np.log(10 / 0.001)
        pyet_penman = (et.penman(Testr.tmean, wind_penman, rs=Testr.rs,
                                 elevation=Testr.elevation, lat=Testr.lat,
                                 tmax=Testr.tmax, tmin=Testr.tmin, rh=Testr.rh,
                                 rhmax=Testr.rhmax, rhmin=Testr.rhmin,
                                 aw=2.626, bw=1.381,
                                 albedo=0.08) * Testr.lambda_corr).round(
            1).tolist()
        r_penman = [6.8, 6.7, 6.7, 6.9, 7.6, 10.1, 7.9, 9.2, 9.4, 7.3]
        self.assertEqual(r_penman, pyet_penman, 1)

    def test_fao56(self):
        wind = Testr.uz * 4.87 / np.log(67.8 * 10 - 5.42)
        pyet_fao56 = et.pm_fao56(Testr.tmean, wind, rs=Testr.rs,
                                 elevation=Testr.elevation, lat=Testr.lat,
                                 tmax=Testr.tmax, tmin=Testr.tmin,
                                 rh=Testr.rh, rhmax=Testr.rhmax,
                                 rhmin=Testr.rhmin
                                 ).round(1).tolist()
        r_fao56 = [5.1, 5., 5., 5.2, 5.9, 8.8, 6.3, 7.8, 8., 5.7]
        self.assertEqual(r_fao56, pyet_fao56, 1)

    def test_pm(self):
        wind = Testr.uz * 4.87 / np.log(67.8 * 10 - 5.42)
        pyet_pm = (et.pm(Testr.tmean, wind, rs=Testr.rs,
                         elevation=Testr.elevation, lat=Testr.lat,
                         tmax=Testr.tmax, tmin=Testr.tmin, rh=Testr.rh,
                         rhmax=Testr.rhmax,
                         rhmin=Testr.rhmin) * Testr.lambda_corr).round(
            1).tolist()
        r_fao56 = [5.1, 5., 5., 5.2, 5.9, 8.8, 6.3, 7.8, 8., 5.7]
        self.assertEqual(r_fao56, pyet_pm, 1)

    def test_asce_pm(self):
        wind = Testr.uz * 4.87 / np.log(67.8 * 10 - 5.42)
        pyet_pm = et.pm_asce(Testr.tmean, wind, rs=Testr.rs,
                             elevation=Testr.elevation, lat=Testr.lat,
                             tmax=Testr.tmax, tmin=Testr.tmin, rh=Testr.rh,
                             rhmax=Testr.rhmax, rhmin=Testr.rhmin).round(
            1).tolist()
        r_fao56 = [5.1, 5., 5., 5.2, 5.9, 8.8, 6.3, 7.8, 8., 5.7]
        self.assertEqual(r_fao56, pyet_pm, 1)

    def test_makkink(self):
        pyet_makk = (et.makkink(Testr.tmean, rs=Testr.rs,
                                elevation=Testr.elevation,
                                k=0.61) * Testr.lambda_corr - 0.12).round(
            1).tolist()
        r_makkink = [3.5, 3.4, 3.5, 3.4, 3.6, 3.8, 3.6, 3.7, 3.7, 3.3]
        self.assertEqual(r_makkink, pyet_makk, 1)

    def test_makkink_knmi(self):
        pyet_makk_knmi = et.makkink_knmi(Testr.tmean, rs=Testr.rs).round(
            1).tolist()
        r_makkink = [3.8, 3.8, 3.9, 3.8, 4.0, 4.3, 4.0, 4.1, 4.1, 3.7]
        self.assertEqual(r_makkink, pyet_makk_knmi, 1)

    def test_pt(self):
        pyet_pt = (et.priestley_taylor(Testr.tmean, rs=Testr.rs, lat=Testr.lat,
                                       elevation=Testr.elevation,
                                       tmax=Testr.tmax, tmin=Testr.tmin,
                                       rhmax=Testr.rhmax, rhmin=Testr.rhmin,
                                       rh=Testr.rh) * Testr.lambda_corr).round(
            1).tolist()
        r_pt = [4., 3.9, 4., 3.9, 4.2, 4.1, 3.9, 3.9, 3.6, 3.7]
        self.assertEqual(r_pt, pyet_pt, 1)

    def test_hargreaves(self):
        pyet_har = (
            et.hargreaves(Testr.tmean, Testr.tmax, Testr.tmin, Testr.lat,
                          method=1)).round(1).tolist()
        r_har = [4.6, 4.3, 4.3, 3.7, 5.4, 4.6, 4.8, 4.4, 4.9, 3.7]
        self.assertEqual(r_har, pyet_har, 1)

    def test_hamon(self):
        pyet_pm = et.hamon(Testr.tmean, lat=Testr.lat, n=Testr.n,
                           tmax=Testr.tmax, tmin=Testr.tmin, method=3).round(
            1).tolist()
        r_hamon = [1.5, 1.4, 1.5, 1.4, 1.8, 2.4, 2., 2.1, 2.2, 1.5]
        self.assertEqual(r_hamon, pyet_pm, 1)

    def test_blaney_criddle(self):
        wind = Testr.uz * 4.87 / np.log(67.8 * 10 - 5.42)
        pyet_bc = (et.blaney_criddle(Testr.tmean, Testr.lat, rhmin=Testr.rhmin,
                                     wind=wind, method=2, n=Testr.n,
                                     clip_zero=False)).round(1).tolist()
        r_bc = [3., 3., 3.1, 2.9, 3.6, 5., 4.1, 4.5, 4.9, 3.3]
        self.assertEqual(r_bc, pyet_bc, 1)

    def test_romanenko(self):
        pyet_romanenko = (
            et.romanenko(Testr.tmean, rh=Testr.rh, tmax=Testr.tmax,
                         tmin=Testr.tmin, rhmax=Testr.rhmax,
                         rhmin=Testr.rhmin)).round(1).tolist()
        r_romanenko = [9.3, 8.9, 9.4, 8.2, 10.6, 16.8, 14., 14.7, 17.4, 9.2]
        self.assertEqual(r_romanenko, pyet_romanenko, 1)

    def test_mcguinnessbordne(self):
        pyet_mgb = (
                et.mcguinness_bordne(Testr.tmean,
                                     lat=Testr.lat) * Testr.lambda_corr).round(
            1).tolist()
        r_mgb = [5.8, 5.5, 5.9, 5.6, 6.4, 7.6, 6.7, 6.8, 7., 5.6]
        self.assertEqual(r_mgb, pyet_mgb, 1)

    def test_jensen_haise(self):
        pyet_jh = (
                et.jensen_haise(Testr.tmean,
                                Testr.rs) * Testr.lambda_corr).round(
            1).tolist()
        r_jh = [5.2, 4.9, 5.3, 5., 5.8, 7., 6.1, 6.3, 6.4, 5.]
        self.assertEqual(r_jh, pyet_jh, 1)

    def test_linacre(self):
        tdew = [10.2375, 8.8375, 11.5125, 10.2875, 11.9875, 9.475, 7.7625,
                7.2375, 3.0125, 11.325]
        pyet_linacre = (
            et.linacre(Testr.tmean, Testr.elevation, Testr.lat,
                       tdew=tdew)).round(1).tolist()
        r_linacre = [4.4, 4.3, 4.4, 4.2, 5.4, 9.1, 7.5, 8., 9.9, 4.3]
        self.assertEqual(r_linacre, pyet_linacre, 1)

    def test_abtew(self):
        pyet_abtew = (
                et.abtew(Testr.tmean, Testr.rs,
                         k=0.52) * Testr.lambda_corr).round(1).tolist()
        r_abtew = [4.3, 4.3, 4.3, 4.3, 4.3, 4.2, 4.2, 4.2, 4.2, 4.1]
        self.assertEqual(r_abtew, pyet_abtew, 1)

    def test_turc(self):
        pyet_turc = (
            et.turc(Testr.tmean, Testr.rs, Testr.rh)).round(1).tolist()
        # There is a mistake in the R-Package
        # The R-Package does not consider RH < 50
        # if pyet.turc also does nonsider this, the result is the same
        r_turc = [4.2, 4., 4.2, 4., 4.3, 6.1, 5.4, 5.5, 6.2, 4.1]
        self.assertEqual(r_turc, pyet_turc, 1)
