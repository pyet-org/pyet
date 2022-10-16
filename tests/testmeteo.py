import unittest

import numpy as np
import pandas as pd

import pyet as et


class Testmeteo(unittest.TestCase):
    def test_calc_rad_sol_in(self):
        # Based on example 10, p 50 TestFAO56
        lat = -22.9 * np.pi / 180
        n = pd.Series(7.1, pd.DatetimeIndex(["2021-5-15"]))
        rad_sol_in = et.calc_rad_sol_in(n, lat)
        self.assertAlmostEqual(float(rad_sol_in.values), 14.5, 1)
