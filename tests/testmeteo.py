import unittest

import numpy as np
import pandas as pd

import pyet as et


class Testmeteo(unittest.TestCase):
    def test_calc_rad_sol_in(self):
        # Based on example 10, p 50 TestFAO56
        tindex = pd.to_datetime("2021-5-15")
        lat = -22.9 * np.pi / 180
        n = 7.1
        rad_sol_in = et.calc_rad_sol_in(tindex, lat, n)
        self.assertAlmostEqual(float(rad_sol_in), 14.5, 1)
