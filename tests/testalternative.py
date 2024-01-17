import unittest

from numpy import full, pi
from pandas import Series, DatetimeIndex, date_range
from xarray import testing, DataArray

import pyet as et


class Testalternative(unittest.TestCase):
    def test_hargreaves(self):
        # Based on example S19.46,p 78 TestFAO56
        tmax = Series([26.6], index=DatetimeIndex(["2015-07-15"]))
        tmin = Series([14.8], index=DatetimeIndex(["2015-07-15"]))
        tmean = (tmax + tmin) / 2
        lat = 45.72 * pi / 180
        har = et.hargreaves(tmean, tmax, tmin, lat)
        self.assertAlmostEqual(float(har.iloc[0]), 5.0, 1)

        tmax_xr = DataArray(
            full((2, 3, 3), 26.6),
            coords=[
                ("time", date_range(start="2015-07-15", periods=2)),
                ("y", [1, 2, 3]),
                ("x", [1, 2, 3]),
            ],
        )
        tmin_xr = DataArray(
            full((2, 3, 3), 14.8),
            coords=[
                ("time", date_range(start="2015-07-15", periods=2)),
                ("y", [1, 2, 3]),
                ("x", [1, 2, 3]),
            ],
        )
        tmean_xr = (tmax_xr + tmin_xr) / 2
        calculated_har = et.hargreaves(tmean_xr, tmax_xr, tmin_xr, lat)
        expected_har = DataArray(
            full((2, 3, 3), 5.0),
            coords=[
                ("time", date_range(start="2015-07-15", periods=2)),
                ("y", [1, 2, 3]),
                ("x", [1, 2, 3]),
            ],
        )
        testing.assert_allclose(calculated_har.round(1), expected_har.round(1))
        lat_xr = DataArray(
            full((3, 3), lat),
            coords=[("y", [1, 2, 3]), ("x", [1, 2, 3])],
        )
        calculated_har1 = et.hargreaves(tmean_xr, tmax_xr, tmin_xr, lat_xr)
        testing.assert_allclose(calculated_har1.round(1), expected_har.round(1))

    def test_makkink(self):
        # Based on example S19.91, p 80 McMahon_etal_2013
        tmean = Series([11.5], index=DatetimeIndex(["1980-07-20"]))
        rs = Series([17.194], index=DatetimeIndex(["1980-07-20"]))
        elevation = 546
        lambda0 = 2.45  # Lambda in McMahon
        lambda1 = et.calc_lambda(tmean)  # Lambda in pyet
        lambda_corr = lambda1 / lambda0
        n = -0.12  # correction for different formula in McMahon_2013
        mak = (et.makkink(tmean, rs, elevation=elevation, k=0.61) * lambda_corr) + n
        self.assertAlmostEqual(float(mak), 2.3928, 2)

        tmean_xr = DataArray(
            full((2, 3, 3), 11.5),
            coords=[
                ("time", date_range(start="2015-07-15", periods=2)),
                ("y", [1, 2, 3]),
                ("x", [1, 2, 3]),
            ],
        )
        rs_xr = DataArray(
            full((2, 3, 3), 17.194),
            coords=[
                ("time", date_range(start="2015-07-15", periods=2)),
                ("y", [1, 2, 3]),
                ("x", [1, 2, 3]),
            ],
        )
        calculated_mak = (
            et.makkink(tmean_xr, rs_xr, elevation=elevation, k=0.61)
            * lambda_corr.values
        ) + n

        expected_mak = DataArray(
            full((2, 3, 3), 2.3928),
            coords=[
                ("time", date_range(start="2015-07-15", periods=2)),
                ("y", [1, 2, 3]),
                ("x", [1, 2, 3]),
            ],
        )
        testing.assert_allclose(calculated_mak.round(2), expected_mak.round(2))
        elevation_xr = DataArray(
            full((3, 3), elevation),
            coords=[("y", [1, 2, 3]), ("x", [1, 2, 3])],
        )
        calculated_mak1 = (
            et.makkink(tmean_xr, rs_xr, elevation=elevation_xr, k=0.61)
            * lambda_corr.values
        ) + n
        testing.assert_allclose(calculated_mak1.round(2), expected_mak.round(2))

    def test_makkink_knmi(self):
        # Based on knmi station 260 (de Bilt)
        tmean = Series([0.9], index=DatetimeIndex(["1980-01-01"]))
        rs = Series([2.53], index=DatetimeIndex(["1980-01-01"]))
        mak = et.makkink_knmi(tmean, rs).round(1)
        self.assertAlmostEqual(float(mak.iloc[0]), 0.3, 2)

        tmean_xr = DataArray(
            full((2, 3, 3), 0.9),
            coords=[
                ("time", date_range(start="1980-01-01", periods=2)),
                ("y", [1, 2, 3]),
                ("x", [1, 2, 3]),
            ],
        )
        rs_xr = DataArray(
            full((2, 3, 3), 2.53),
            coords=[
                ("time", date_range(start="1980-01-01", periods=2)),
                ("y", [1, 2, 3]),
                ("x", [1, 2, 3]),
            ],
        )
        calculated_mak1 = et.makkink_knmi(tmean_xr, rs_xr)

        expected_mak1 = DataArray(
            full((2, 3, 3), 0.3),
            coords=[
                ("time", date_range(start="1980-01-01", periods=2)),
                ("y", [1, 2, 3]),
                ("x", [1, 2, 3]),
            ],
        )
        testing.assert_allclose(calculated_mak1.round(1), expected_mak1.round(1))

    def test_blaney_criddle(self):
        # Based on example S19.93, McMahon_etal_2013
        tmean = Series([11.5], index=DatetimeIndex(["1980-07-20"]))
        rhmin = Series([25], index=DatetimeIndex(["1980-07-20"]))
        wind = Series([0.5903], index=DatetimeIndex(["1980-07-20"]))
        n = Series([10.7], index=DatetimeIndex(["1980-07-20"]))
        py = Series([0.2436], index=DatetimeIndex(["1980-07-20"]))
        lat = -23.7951 * pi / 180
        bc = et.blaney_criddle(tmean, lat, wind=wind, n=n, rhmin=rhmin, py=py, method=2)
        self.assertAlmostEqual(float(bc.iloc[0]), 3.1426, 2)

        tmean_xr = DataArray(
            full((2, 3, 3), 11.5),
            coords=[
                ("time", date_range(start="1980-07-20", periods=2)),
                ("y", [1, 2, 3]),
                ("x", [1, 2, 3]),
            ],
        )
        rhmin_xr = DataArray(
            full((2, 3, 3), 25),
            coords=[
                ("time", date_range(start="1980-07-20", periods=2)),
                ("y", [1, 2, 3]),
                ("x", [1, 2, 3]),
            ],
        )
        wind_xr = DataArray(
            full((2, 3, 3), 0.5903),
            coords=[
                ("time", date_range(start="1980-07-20", periods=2)),
                ("y", [1, 2, 3]),
                ("x", [1, 2, 3]),
            ],
        )
        n_xr = DataArray(
            full((2, 3, 3), 10.7),
            coords=[
                ("time", date_range(start="1980-07-20", periods=2)),
                ("y", [1, 2, 3]),
                ("x", [1, 2, 3]),
            ],
        )
        py_xr = DataArray(
            full((2, 3, 3), 0.2436),
            coords=[
                ("time", date_range(start="1980-07-20", periods=2)),
                ("y", [1, 2, 3]),
                ("x", [1, 2, 3]),
            ],
        )
        calculated_bc = et.blaney_criddle(
            tmean_xr, lat, wind=wind_xr, n=n_xr, rhmin=rhmin_xr, py=py_xr, method=2
        )
        expected_bc = DataArray(
            full((2, 3, 3), 3.1426),
            coords=[
                ("time", date_range(start="1980-07-20", periods=2)),
                ("y", [1, 2, 3]),
                ("x", [1, 2, 3]),
            ],
        )
        testing.assert_allclose(calculated_bc.round(2), expected_bc.round(2))
        lat_xr = DataArray(
            full((3, 3), lat),
            coords=[("y", [1, 2, 3]), ("x", [1, 2, 3])],
        )
        calculated_bc1 = et.blaney_criddle(
            tmean_xr, lat_xr, wind=wind_xr, n=n_xr, rhmin=rhmin_xr, py=py_xr, method=2
        )
        testing.assert_allclose(calculated_bc1.round(2), expected_bc.round(2))

    def test_turc(self):
        # Based on example S19.99, McMahon_etal_2013
        tmean = Series([11.5], index=DatetimeIndex(["1980-07-20"]))
        rhmean = Series([48], index=DatetimeIndex(["1980-07-20"]))
        rs = Series([17.194], index=DatetimeIndex(["1980-07-20"]))
        turc = et.turc(tmean, rs, rhmean)
        self.assertAlmostEqual(float(turc.iloc[0]), 2.6727, 1)

        tmean_xr = DataArray(
            full((2, 3, 3), 11.5),
            coords=[
                ("time", date_range(start="1980-07-20", periods=2)),
                ("y", [1, 2, 3]),
                ("x", [1, 2, 3]),
            ],
        )
        rh_xr = DataArray(
            full((2, 3, 3), 48),
            coords=[
                ("time", date_range(start="1980-07-20", periods=2)),
                ("y", [1, 2, 3]),
                ("x", [1, 2, 3]),
            ],
        )
        rs_xr = DataArray(
            full((2, 3, 3), 17.194),
            coords=[
                ("time", date_range(start="1980-07-20", periods=2)),
                ("y", [1, 2, 3]),
                ("x", [1, 2, 3]),
            ],
        )
        calculated_turc1 = et.turc(tmean_xr, rs_xr, rh_xr)

        expected_turc1 = DataArray(
            full((2, 3, 3), 2.6727),
            coords=[
                ("time", date_range(start="1980-07-20", periods=2)),
                ("y", [1, 2, 3]),
                ("x", [1, 2, 3]),
            ],
        )
        testing.assert_allclose(calculated_turc1.round(1), expected_turc1.round(1))

    def test_hargreaves_samani(self):
        # Based on example S19.101, McMahon_etal_2013
        tmean = Series([11.5], index=DatetimeIndex(["1980-07-20"]))
        tmax = Series([21], index=DatetimeIndex(["1980-07-20"]))
        tmin = Series([2], index=DatetimeIndex(["1980-07-20"]))
        lambda0 = 2.45  # Lambda in McMahon
        lambda1 = et.calc_lambda(tmean)  # Lambda in pyet
        lambda_corr = lambda1 / lambda0
        lat = -23.7951 * pi / 180
        har = et.hargreaves(tmean, tmax, tmin, lat, method=1)
        self.assertAlmostEqual(float(har.iloc[0] * lambda_corr), 4.1129, 2)

        tmean_xr = DataArray(
            full((2, 3, 3), 11.5),
            coords=[
                ("time", date_range(start="1980-07-20", periods=2)),
                ("y", [1, 2, 3]),
                ("x", [1, 2, 3]),
            ],
        )
        tmax_xr = DataArray(
            full((2, 3, 3), 21),
            coords=[
                ("time", date_range(start="1980-07-20", periods=2)),
                ("y", [1, 2, 3]),
                ("x", [1, 2, 3]),
            ],
        )
        tmin_xr = DataArray(
            full((2, 3, 3), 2),
            coords=[
                ("time", date_range(start="1980-07-20", periods=2)),
                ("y", [1, 2, 3]),
                ("x", [1, 2, 3]),
            ],
        )
        calculated_har = (
            et.hargreaves(tmean_xr, tmax_xr, tmin_xr, lat, method=1)
            * lambda_corr.values
        )
        expected_har1 = DataArray(
            full((2, 3, 3), 4.1129),
            coords=[
                ("time", date_range(start="1980-07-20", periods=2)),
                ("y", [1, 2, 3]),
                ("x", [1, 2, 3]),
            ],
        )
        expected_har1.loc[{"time": "1980-07-21"}] = 4.13
        testing.assert_allclose(calculated_har.round(2), expected_har1.round(2))
        lat_xr = DataArray(
            full((3, 3), lat),
            coords=[("y", [1, 2, 3]), ("x", [1, 2, 3])],
        )
        calculated_har1 = (
            et.hargreaves(tmean_xr, tmax_xr, tmin_xr, lat_xr, method=1)
            * lambda_corr.values
        )
        testing.assert_allclose(calculated_har1.round(2), expected_har1.round(2))

    def test_priestley_taylor(self):
        # Based on example S19.109, McMahon_etal_2013
        tmean = Series([11.5], index=DatetimeIndex(["1980-07-20"]))
        rn = Series([8.6401], index=DatetimeIndex(["1980-07-20"]))
        elevation = 546
        lambda0 = 2.45  # Lambda in McMahon
        lambda1 = et.calc_lambda(tmean)  # Lambda in pyet
        lambda_corr = lambda1 / lambda0
        pt = et.priestley_taylor(tmean, rn=rn, elevation=elevation, alpha=1.26)
        self.assertAlmostEqual(float(pt.iloc[0] * lambda_corr), 2.6083, 2)

        tmean_xr = DataArray(
            full((2, 3, 3), 11.5),
            coords=[
                ("time", date_range(start="1980-07-20", periods=2)),
                ("y", [1, 2, 3]),
                ("x", [1, 2, 3]),
            ],
        )
        rn_xr = DataArray(
            full((2, 3, 3), 8.6401),
            coords=[
                ("time", date_range(start="1980-07-20", periods=2)),
                ("y", [1, 2, 3]),
                ("x", [1, 2, 3]),
            ],
        )
        calculated_pt = (
            et.priestley_taylor(tmean_xr, rn=rn_xr, elevation=elevation, alpha=1.26)
            * lambda_corr.values
        )

        expected_pt1 = DataArray(
            full((2, 3, 3), 2.6083),
            coords=[
                ("time", date_range(start="1980-07-20", periods=2)),
                ("y", [1, 2, 3]),
                ("x", [1, 2, 3]),
            ],
        )
        testing.assert_allclose(calculated_pt.round(2), expected_pt1.round(2))
        elevation_xr = DataArray(
            full((3, 3), elevation),
            coords=[("y", [1, 2, 3]), ("x", [1, 2, 3])],
        )
        calculated_pt1 = (
            et.priestley_taylor(tmean_xr, rn=rn_xr, elevation=elevation_xr, alpha=1.26)
            * lambda_corr.values
        )
        testing.assert_allclose(calculated_pt1.round(2), expected_pt1.round(2))

    def test_blaney_criddle1(self):
        # Based on example 5.1, P89 Schrodter 1985
        tmean = Series([17.3], index=DatetimeIndex(["1980-07-20"]))
        lat = 50 * pi / 180
        bc = et.blaney_criddle(tmean, lat, method=0)
        self.assertAlmostEqual(float(bc.iloc[0]), 3.9, 1)

        tmean_xr = DataArray(
            full((2, 3, 3), 17.3),
            coords=[
                ("time", date_range(start="1980-07-20", periods=2)),
                ("y", [1, 2, 3]),
                ("x", [1, 2, 3]),
            ],
        )
        calculated_bc1 = et.blaney_criddle(tmean_xr, lat, method=0)

        expected_bc1 = DataArray(
            full((2, 3, 3), 3.9),
            coords=[
                ("time", date_range(start="1980-07-20", periods=2)),
                ("y", [1, 2, 3]),
                ("x", [1, 2, 3]),
            ],
        )
        testing.assert_allclose(calculated_bc1.round(1), expected_bc1.round(1))
        lat_xr = DataArray(
            full((3, 3), lat),
            coords=[("y", [1, 2, 3]), ("x", [1, 2, 3])],
        )
        calculated_bc2 = et.blaney_criddle(tmean_xr, lat_xr, method=0)
        testing.assert_allclose(calculated_bc2.round(1), expected_bc1.round(1))

    def test_haude(self):
        # Based on example 5.2, P95 Schrodter 1985
        tmean = Series([21.5], index=DatetimeIndex(["1980-07-20"]))
        ea = Series([1.19], index=DatetimeIndex(["1980-07-20"]))
        k = 0.26 / 0.35  # 0.35 value is taken in pyet
        e0 = et.calc_e0(tmean)
        rh = ea / e0 * 100
        haude = et.haude(tmean, rh=rh, k=k)
        self.assertAlmostEqual(float(haude.iloc[0]), 3.6, 1)

        tmean_xr = DataArray(
            full((2, 3, 3), 21.5),
            coords=[
                ("time", date_range(start="1980-07-20", periods=2)),
                ("y", [1, 2, 3]),
                ("x", [1, 2, 3]),
            ],
        )
        rh_xr = DataArray(
            full((2, 3, 3), rh.iloc[0]),
            coords=[
                ("time", date_range(start="1980-07-20", periods=2)),
                ("y", [1, 2, 3]),
                ("x", [1, 2, 3]),
            ],
        )
        calculated_haude1 = et.haude(tmean_xr, rh=rh_xr, k=k)

        expected_haude1 = DataArray(
            full((2, 3, 3), 3.6),
            coords=[
                ("time", date_range(start="1980-07-20", periods=2)),
                ("y", [1, 2, 3]),
                ("x", [1, 2, 3]),
            ],
        )
        testing.assert_allclose(calculated_haude1.round(1), expected_haude1.round(1))
