import unittest
import numpy as np
import pandas as pd
import pyet as et

# Data from "THE ASCE STANDARDIZED REFERENCE EVAPOTRANSPIRATION EQUATION"
# Appendix C page C-7
LAT = 0.7053  # in radians
ELEV = 1462.4
LZ = 105  # deg. west of greenwich
LM = 105  # deg. west of greenwich

tmean = np.array(
    [30.9, 31.2, 29.1, 28.3, 26, 22.9, 20.1, 19.9, 18.4, 16.5, 15.4, 15.5,
     13.5, 13.2, 16.2, 20, 22.9, 26.4, 28.2, 29.8, 30.9, 31.8, 32.5, 32.9,
     32.4, 30.2, 30.6, 28.3, 25.9, 23.9])
ea = [1.09, 1.15, 1.21, 1.21, 1.13, 1.2, 1.35, 1.35, 1.32, 1.26, 1.34, 1.31,
      1.26, 1.24, 1.31, 1.36, 1.39, 1.25, 1.17, 1.03, 1.02, 0.98, 0.87, 0.86,
      0.93, 1.14, 1.27, 1.27, 1.17, 1.2]
rs = [2.24, 1.65, 0.34, 0.32, 0.08, 0, 0, 0, 0, 0, 0, 0, 0, 0.03, 0.46, 1.09,
      1.74, 2.34, 2.84, 3.25, 3.21, 3.34, 2.96, 2.25, 1.35, 0.88, 0.79, 0.27,
      0.03, 0]
wind = [3.74, 3.3, 1.06, 2.8, 2.03, 0.96, 0.53, 0.87, 0.28, 0.46, 0.92, 0.63,
        0.63, 0.27, 1.15, 1.18, 0.81, 0.67, 1.4, 1.81, 1.9, 2.54, 2.67,
        2.85, 2.55, 3.14, 2.56, 2.72, 3.01, 2.64]
hindex = pd.date_range("2000-7-1 16:00", "2000-7-2 21:00", freq="H")
rs_series = pd.Series(rs, hindex)
tmean_series = pd.Series(tmean, hindex)


class TestASCE(unittest.TestCase):
    def test_day_of_year_hour(self):
        doy = et.day_of_year(hindex).tolist()
        doyr = [183, 183, 183, 183, 183, 183, 183, 183, 184, 184, 184, 184,
                184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184, 184,
                184, 184, 184, 184, 184, 184]
        self.assertEqual(doy, doyr, 1)

    def test_vpc_hour(self):
        vpc1 = et.calc_vpc(tmean).round(2).tolist()
        vpcr = np.array(
            [0.2548, 0.2586, 0.2329, 0.2237, 0.1989, 0.1692, 0.1457, 0.1441,
             0.1328, 0.1195, 0.1124, 0.113, 0.1009, 0.0992, 0.1175, 0.1449,
             0.1692, 0.203, 0.2226, 0.2412, 0.2548, 0.2664, 0.2757, 0.2811,
             0.2743, 0.2461, 0.251, 0.2237, 0.1979, 0.1783]).round(
            2).tolist()
        self.assertEqual(vpc1, vpcr, 1)

    def test_e0_calc1(self):
        # Based on Example 3 and 4, p. 36 FAO.
        e0 = et.calc_e0(tmean).round(1).tolist()
        e0_ref = np.array(
            [4.467, 4.544, 4.029, 3.846, 3.361, 2.792, 2.353, 2.324,
             2.116, 1.877, 1.74, 1.761, 1.547, 1.517, 1.842, 2.338, 2.792,
             3.442, 3.824, 4.195, 4.467, 4.701, 4.891, 5.002, 4.863,
             4.292, 4.391, 3.846, 3.342, 2.966]).round(1).tolist()
        self.assertEqual(e0, e0_ref, 1)

    def test_relative_distance_hour(self):
        rd = et.relative_distance(184)
        self.assertAlmostEqual(rd, 0.967, 3)

    def test_solar_declination(self):
        sd = et.solar_declination(184)
        self.assertAlmostEqual(sd, 0.4017, 2)

    def test_omega1_omega2(self):
        omega1, omega2 = et.sunset_angle_hour(hindex, LAT, LZ, LM)
        omega1_ref = np.array(
            [0.773, 1.035, 1.297, 1.558, 1.82, 1.941, 1.941, 1.941, 1.941,
             -1.939, -1.939, -1.939, -1.939, -1.939, -1.846, -1.584,
             -1.322, -1.06, -0.799, - 0.537, - 0.275, - 0.013, 0.249,
             0.51, 0.772, 1.034, 1.296, 1.558, 1.819, 1.939]).round(1).tolist()
        omega2_ref = np.array(
            [1.035, 1.297, 1.558, 1.82, 1.941, 1.941, 1.941, 1.941, 1.941,
             -1.939, -1.939, -1.939, -1.939, -1.846, -1.584, -1.322,
             -1.06, -0.799, - 0.537, - 0.275, - 0.013, 0.249, 0.51,
             0.772, 1.034, 1.296, 1.558, 1.819, 1.939, 1.939]).round(
            1).tolist()
        self.assertEqual(omega1.round(1).tolist(), omega1_ref, 3)
        self.assertEqual(omega2.round(1).tolist(), omega2_ref, 3)

    def test_extra_rad_hour(self):
        ra_hour = et.extraterrestrial_r_hour(hindex, LAT, LZ, LM)
        ra_ref = np.array(
            [3.26, 2.52, 1.68, 0.81, 0.09, 0, 0, 0, 0, 0, 0, 0, 0, 0.05, 0.72,
             1.59, 2.43, 3.19, 3.81, 4.26, 4.49, 4.51, 4.29, 3.87, 3.26, 2.52,
             1.68, 0.81, 0.09, 0]).round(1).tolist()
        self.assertEqual(ra_hour.round(1).tolist(), ra_ref, 3)

    def test_rad_rso_hour(self):
        ra_hour = et.extraterrestrial_r_hour(hindex, LAT, LZ, LM)
        rad_rso = et.calc_rso(ra_hour, ELEV, kab=0.779)
        rso_ref = np.array(
            [2.54, 1.96, 1.31, 0.63, 0.07, 0, 0, 0, 0, 0, 0, 0, 0, 0.04, 0.56,
             1.24, 1.9, 2.49, 2.97, 3.32, 3.5, 3.51, 3.34, 3.01, 2.54, 1.96,
             1.31, 0.63, 0.07, 0]).round(1).tolist()
        self.assertEqual(rad_rso.round(1).tolist(), rso_ref, 3)

#    def test_rad_rnl_hour(self):
#        rad_long = et.calc_rad_long(rs_series, tmean=tmean, ea=ea,
#                                    elevation=ELEV, lat=LAT, lz=LZ, lon=LM,
#                                    kab=0.779, freq="H").round(1)
#        rlong_ref = np.array(
#            [0.284, 0.262, 0, 0.104, 0.313, 0.23, 0.211, 0.21, 0.208, 0.207,
#             0.198, 0.201, 0.198, 0.161, 0.193, 0.224, 0.245, 0.278, 0.299,
#             0.331, 0.308, 0.333, 0.317, 0.248, 0.134, 0.084, 0.147,
#             0.07, 0.077, 0.076]).round(1).tolist()
#        self.assertEqual(rad_long.round(1).tolist(), rlong_ref, 3)
# Check
#    def test_rad_rn_hour(self):
#        rad_long = et.calc_rad_long(rs_series, tmean=tmean, ea=ea,
#                                    elevation=ELEV, lat=LAT, lz=LZ, lon=LM,
#                                    kab=0.779, freq="H")
#        rad_rn = (1 - 0.23) * rs_series - rad_long
#        rneto_ref = np.array(
#            [1.44, 1.009, 0.262, 0.142, -0.251, -0.23, -0.211, -0.21, -0.208,
#             -0.207, -0.198, -0.201, -0.198, -0.138, 0.161, 0.616, 1.095,
#             1.524, 1.888, 2.171, 2.164, 2.239, 1.962, 1.485, 0.905, 0.593,
#             0.461, 0.138, -0.054, -0.076]).round(1).tolist()
#        self.assertEqual(rad_rn.round(1).tolist(), rneto_ref, 3)

    def test_etos_hour(self):
        etos = et.pm_asce(tmean_series, wind, rs=rs_series, elevation=ELEV,
                          lat=LAT, lz=LZ, lon=LM, kab=0.779, ea=ea, freq="H")
        etos_ref = np.array(
            [0.61, 0.48, 0.14, 0.2, 0.06, 0.01, -0.01, 0, -0.02, -0.01, -0.01,
             -0.01, -0.01, -0.01, 0.06, 0.19, 0.32, 0.46, 0.6, 0.72, 0.73,
             0.79, 0.74, 0.62, 0.44, 0.35, 0.29, 0.19, 0.1, 0.08,
             ]).round(1).tolist()
        self.assertEqual(etos.round(1).tolist(), etos_ref, 3)

    def test_etrs_hour(self):
        etrs = et.pm_asce(tmean_series, wind, rs=rs_series, elevation=ELEV,
                          lat=LAT, lz=LZ, lon=LM, kab=0.779, ea=ea, freq="H",
                          etype="rs")
        etrs_ref = np.array(
            [0.82, 0.66, 0.2, 0.33, 0.09, 0.02, -0.01, 0, -0.03, -0.02,
             -0.01, -0.02, -0.02, -0.02, 0.08, 0.23, 0.37, 0.52, 0.7, 0.86,
             0.88, 0.97, 0.93, 0.81, 0.6, 0.52, 0.42, 0.31, 0.14, 0.11]).round(
            1).tolist()
        self.assertEqual(etrs.round(1).tolist(), etrs_ref, 3)
