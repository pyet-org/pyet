import unittest
import numpy as np
import pandas as pd

import pyet as et


class TestFAO56(unittest.TestCase):
    def test_press_calc(self):
        # Based on Example 2 and 4, p. 32 FAO.
        p1 = et.calc_press(1800)
        p2 = et.calc_press(1200)
        self.assertAlmostEqual(p1, 81.8, 1)
        self.assertAlmostEqual(p2, 87.9, 1)
        # Basen on Table C-2 in ASCE(2001).
        p3 = et.calc_press(1462.4)
        self.assertAlmostEqual(p3, 85.17, 1)

    def test_vpc(self):
        # Based on ASCE Table C-3
        tmean = np.array([21.65, 22.9, 23.7, 22.8, 24.3,
                          26.0, 26.1, 26.4, 23.9, 24.2])
        vpc1 = et.calc_vpc(tmean).round(3).tolist()
        vpcr = np.array([0.1585, 0.1692, 0.1762, 0.1684, 0.1820,
                         0.199, 0.1996, 0.2027, 0.1781, 0.1809]).round(
            3).tolist()
        self.assertAlmostEqual(vpc1, vpcr, 1)

    def test_psy_calc(self):
        # Based on Example 2, p. 32 FAO.
        psy = et.calc_psy(81.8)
        self.assertAlmostEqual(psy, 0.054, 3)
        # Basen on Table C-2 in ASCE(2001).
        psy1 = et.calc_psy(85.17)
        self.assertAlmostEqual(psy1, 0.0566, 3)

    def test_e0_calc(self):
        # Based on Example 3 and 4, p. 36 FAO.
        e0_1 = et.calc_e0(24.5)
        e0_2 = et.calc_e0(15)
        e0_3 = et.calc_e0(19.75)
        e0_4 = et.calc_e0(19.5)
        self.assertAlmostEqual(e0_1, 3.075, 2)
        self.assertAlmostEqual(e0_2, 1.705, 2)
        self.assertAlmostEqual(e0_3, 2.3, 2)
        self.assertAlmostEqual(e0_4, 2.267, 2)

    def test_es_calc(self):
        # Based on Example 3 and 6, p. 36 FAO.
        es = et.calc_es(tmax=24.5, tmin=15)
        self.assertAlmostEqual(es, 2.39, 3)

    def test_ea_calc(self):
        # Based on Example 5, p. 39 FAO.
        ea1 = et.calc_ea(tmax=25, tmin=18, rhmax=82, rhmin=54)
        rhmean = (82 + 54) / 2
        ea2 = et.calc_ea(tmax=25, tmin=18, rh=rhmean)
        self.assertAlmostEqual(ea1, 1.70, 2)
        self.assertAlmostEqual(ea2, 1.78, 2)

    def test_relative_distance(self):
        # Based on Example 8, p. 47 FAO.
        rd = et.relative_distance(246)
        self.assertAlmostEqual(rd, 0.985, 3)

    def test_solar_declination(self):
        # Based on Example 8, p. 47 FAO.
        sd = et.solar_declination(246)
        self.assertAlmostEqual(sd, 0.12, 2)

    def test_sunset_angle(self):
        # Based on Example 8, p. 47 FAO.
        sangle = et.sunset_angle(0.12, -0.35)
        self.assertAlmostEqual(sangle, 1.527, 3)

    def test_day_of_year(self):
        # Based on Example 8, p. 47 FAO.
        doy = et.day_of_year(pd.to_datetime("2015-09-03"))
        self.assertAlmostEqual(doy, 246, 1)
        # Based on ASCD Table C-3
        dindex = pd.date_range("2020-07-01", "2020-07-10")
        doy1 = et.day_of_year(dindex).tolist()
        doyr = [183, 184, 185, 186, 187,
                188, 189, 190, 191, 192]
        self.assertEqual(doy1, doyr, 1)

    def test_extraterrestrial_r(self):
        # Based on Example 8, p. 47 FAO.
        extrar = et.extraterrestrial_r(pd.to_datetime("2015-09-03"), -0.35)
        self.assertAlmostEqual(extrar, 32.2, 1)

    def test_daylight_hours(self):
        # Based on Example 9, p. 47 FAO.
        dayhours = et.daylight_hours(pd.to_datetime("2015-09-03"), -0.35)
        self.assertAlmostEqual(dayhours, 11.7, 1)

    def test_calc_rad_long(self):
        # Based on Example 10, p. 52 FAO.
        Rs = pd.Series([14.5], index=pd.DatetimeIndex(["2015-05-15"]))
        Rnl = et.calc_rad_long(Rs, tmax=25.1, tmin=19, ea=2.1, rso=18.8)
        self.assertAlmostEqual(float(Rnl), 3.6, 1)

    def test_et_fao56(self):
        # Based on Example 18, p. 72 FAO.
        wind = pd.Series([2.078], index=pd.DatetimeIndex(["2015-07-06"]))
        tmax = pd.Series([21.5], index=pd.DatetimeIndex(["2015-07-06"]))
        tmin = pd.Series([12.3], index=pd.DatetimeIndex(["2015-07-06"]))
        tmean = (tmax + tmin) / 2
        rhmax = pd.Series([84], index=pd.DatetimeIndex(["2015-07-06"]))
        rhmin = pd.Series([63], index=pd.DatetimeIndex(["2015-07-06"]))
        Rs = pd.Series([22.07], index=pd.DatetimeIndex(["2015-07-06"]))
        n = 9.25
        nn = 16.1
        elevation = 100
        lat = 50.80 * np.pi / 180
        et56 = et.pm_fao56(tmean, wind, elevation=elevation, lat=lat, rs=Rs,
                           tmax=tmax, tmin=tmin, rhmax=rhmax,
                           rhmin=rhmin, n=n, nn=nn)
        self.assertAlmostEqual(float(et56), 3.9, 1)
