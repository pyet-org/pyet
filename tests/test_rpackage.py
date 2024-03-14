import unittest

from numpy import full, log, tile, newaxis, asarray
import pandas as pd
from xarray import DataArray, testing

import pyet as et


class Testr(unittest.TestCase):
    lat = -0.6095
    elevation = 48

    tmax = pd.Series(
        [28.8, 27.4, 29.0, 26.3, 32.7, 35.7, 33.8, 33.8, 35.6, 27.5],
        index=pd.date_range("2001-03-01", "2001-03-10", freq="d"),
    )
    tmin = pd.Series(
        [15.1, 14.0, 16.3, 16.2, 17.8, 26.7, 20.6, 22.3, 22.7, 17.1],
        index=pd.date_range("2001-03-01", "2001-03-10", freq="d"),
    )
    tmean = (tmax + tmin) / 2
    rhmax = pd.Series(
        [68, 77, 69, 70, 73, 36, 51, 46, 34, 67],
        index=pd.date_range("2001-03-01", "2001-03-10", freq="d"),
    )
    rhmin = pd.Series(
        [30, 25, 30, 34, 29, 20, 18, 18, 10, 30],
        index=pd.date_range("2001-03-01", "2001-03-10", freq="d"),
    )
    rh = (rhmax + rhmin) / 2
    rs = pd.Series(
        [
            20.44477761,
            20.34991941,
            20.25373437,
            20.1562439,
            20.05747057,
            19.95743816,
            19.85617161,
            19.75369701,
            19.65004162,
            19.54523386,
        ],
        index=pd.date_range("2001-03-01", "2001-03-10", freq="d"),
    )
    uz = pd.Series(
        [
            2.65625,
            2.78472222,
            2.49305556,
            3.73611111,
            3.30902778,
            4.93402778,
            2.90972222,
            4.48958333,
            3.99305556,
            4.31597222,
        ],
        index=pd.date_range("2001-03-01", "2001-03-10", freq="d"),
    )
    n = pd.Series(
        [8.6, 8.6, 8.6, 8.6, 8.6, 8.6, 8.6, 8.6, 8.6, 8.6],
        index=pd.date_range("2001-03-01", "2001-03-10", freq="d"),
    )
    lambda0 = 2.45  # Lambda in Danlu R Evapotranspiration package
    lambda1 = et.calc_lambda(tmean)  # Lambda in pyet
    lambda_corr = lambda1 / lambda0

    # Number of x and y dimensions
    x_dim, y_dim = 3, 3
    lat_xr = DataArray(
        full((3, 3), lat),
        coords=[("y", range(y_dim)), ("x", range(x_dim))],
    )
    elevation_xr = DataArray(
        full((3, 3), elevation),
        coords=[("y", range(y_dim)), ("x", range(x_dim))],
    )
    # Covert pandas Series to xarray DataArray
    tmax_xr = tile(tmax.values[:, newaxis, newaxis], (1, x_dim, y_dim))
    tmax_xr = DataArray(
        tmax_xr,
        coords={"time": tmax.index, "y": range(y_dim), "x": range(x_dim)},
        dims=["time", "y", "x"],
    )
    # Covert pandas Series to xarray DataArray
    tmin_xr = tile(tmin.values[:, newaxis, newaxis], (1, x_dim, y_dim))
    tmin_xr = DataArray(
        tmin_xr,
        coords={"time": tmin.index, "y": range(y_dim), "x": range(x_dim)},
        dims=["time", "y", "x"],
    )
    tmean_xr = (tmax_xr + tmin_xr) / 2
    rhmax_xr = tile(rhmax.values[:, newaxis, newaxis], (1, x_dim, y_dim))
    rhmax_xr = DataArray(
        rhmax_xr,
        coords={"time": rhmax.index, "y": range(y_dim), "x": range(x_dim)},
        dims=["time", "y", "x"],
    )
    rhmin_xr = tile(rhmin.values[:, newaxis, newaxis], (1, x_dim, y_dim))
    rhmin_xr = DataArray(
        rhmin_xr,
        coords={"time": rhmin.index, "y": range(y_dim), "x": range(x_dim)},
        dims=["time", "y", "x"],
    )
    rh_xr = (rhmax_xr + rhmin_xr) / 2
    rs_xr = tile(rs.values[:, newaxis, newaxis], (1, x_dim, y_dim))
    rs_xr = DataArray(
        rs_xr,
        coords={"time": rs.index, "y": range(y_dim), "x": range(x_dim)},
        dims=["time", "y", "x"],
    )
    uz_xr = tile(uz.values[:, newaxis, newaxis], (1, x_dim, y_dim))
    uz_xr = DataArray(
        uz_xr,
        coords={"time": uz.index, "y": range(y_dim), "x": range(x_dim)},
        dims=["time", "y", "x"],
    )
    n_xr = tile(n.values[:, newaxis, newaxis], (1, x_dim, y_dim))
    n_xr = DataArray(
        n_xr,
        coords={"time": n.index, "y": range(y_dim), "x": range(x_dim)},
        dims=["time", "y", "x"],
    )

    def test_penman(self):
        wind_penman = Testr.uz * log(2 / 0.001) / log(10 / 0.001)
        pyet_penman = (
            (
                et.penman(
                    Testr.tmean,
                    wind_penman,
                    rs=Testr.rs,
                    elevation=Testr.elevation,
                    lat=Testr.lat,
                    tmax=Testr.tmax,
                    tmin=Testr.tmin,
                    rh=Testr.rh,
                    rhmax=Testr.rhmax,
                    rhmin=Testr.rhmin,
                    aw=2.626,
                    bw=1.381,
                    albedo=0.08,
                )
                * Testr.lambda_corr
            )
            .round(1)
            .tolist()
        )
        r_penman = [6.8, 6.7, 6.7, 6.9, 7.6, 10.1, 7.9, 9.2, 9.4, 7.3]
        self.assertEqual(r_penman, pyet_penman, 1)

        expected_penman = tile(asarray(r_penman)[:, newaxis, newaxis], (1, 3, 3))
        expected_penman = DataArray(
            expected_penman,
            coords={"time": Testr.tmean.index, "y": range(3), "x": range(3)},
            dims=["time", "y", "x"],
        )
        wind_penman_xr = Testr.uz_xr * log(2 / 0.001) / log(10 / 0.001)
        penman_xr = (
            et.penman(
                Testr.tmean_xr,
                wind_penman_xr,
                rs=Testr.rs_xr,
                elevation=Testr.elevation,
                lat=Testr.lat,
                tmax=Testr.tmax_xr,
                tmin=Testr.tmin_xr,
                rh=Testr.rh_xr,
                rhmax=Testr.rhmax_xr,
                rhmin=Testr.rhmin_xr,
                aw=2.626,
                bw=1.381,
                albedo=0.08,
            )
            * Testr.lambda_corr.values[:, None, None]
        )
        testing.assert_allclose(penman_xr.round(1), expected_penman.round(1))
        penman_xr1 = (
            et.penman(
                Testr.tmean_xr,
                wind_penman_xr,
                rs=Testr.rs_xr,
                elevation=Testr.elevation_xr,
                lat=Testr.lat_xr,
                tmax=Testr.tmax_xr,
                tmin=Testr.tmin_xr,
                rh=Testr.rh_xr,
                rhmax=Testr.rhmax_xr,
                rhmin=Testr.rhmin_xr,
                aw=2.626,
                bw=1.381,
                albedo=0.08,
            )
            * Testr.lambda_corr.values[:, None, None]
        )
        testing.assert_allclose(penman_xr1.round(1), expected_penman.round(1))

    def test_fao56(self):
        wind = Testr.uz * 4.87 / log(67.8 * 10 - 5.42)
        pyet_fao56 = (
            et.pm_fao56(
                Testr.tmean,
                wind,
                rs=Testr.rs,
                elevation=Testr.elevation,
                lat=Testr.lat,
                tmax=Testr.tmax,
                tmin=Testr.tmin,
                rh=Testr.rh,
                rhmax=Testr.rhmax,
                rhmin=Testr.rhmin,
            )
            .round(1)
            .tolist()
        )
        r_fao56 = [5.1, 5.0, 5.0, 5.2, 5.9, 8.8, 6.3, 7.8, 8.0, 5.7]
        self.assertEqual(r_fao56, pyet_fao56, 1)

        expected_fao56 = tile(asarray(r_fao56)[:, newaxis, newaxis], (1, 3, 3))
        expected_fao56 = DataArray(
            expected_fao56,
            coords={"time": Testr.tmean.index, "y": range(3), "x": range(3)},
            dims=["time", "y", "x"],
        )
        wind_xr = Testr.uz_xr * 4.87 / log(67.8 * 10 - 5.42)
        calculated_fao56 = et.pm_fao56(
            Testr.tmean_xr,
            wind_xr,
            rs=Testr.rs_xr,
            elevation=Testr.elevation,
            lat=Testr.lat,
            tmax=Testr.tmax_xr,
            tmin=Testr.tmin_xr,
            rh=Testr.rh_xr,
            rhmax=Testr.rhmax_xr,
            rhmin=Testr.rhmin_xr,
        )
        testing.assert_allclose(calculated_fao56.round(1), expected_fao56.round(1))
        calculated_fao561 = et.pm_fao56(
            Testr.tmean_xr,
            wind_xr,
            rs=Testr.rs_xr,
            elevation=Testr.elevation_xr,
            lat=Testr.lat_xr,
            tmax=Testr.tmax_xr,
            tmin=Testr.tmin_xr,
            rh=Testr.rh_xr,
            rhmax=Testr.rhmax_xr,
            rhmin=Testr.rhmin_xr,
        )
        testing.assert_allclose(calculated_fao561.round(1), expected_fao56.round(1))

    def test_pm(self):
        wind = Testr.uz * 4.87 / log(67.8 * 10 - 5.42)
        pyet_pm = (
            (
                et.pm(
                    Testr.tmean,
                    wind,
                    rs=Testr.rs,
                    elevation=Testr.elevation,
                    lat=Testr.lat,
                    tmax=Testr.tmax,
                    tmin=Testr.tmin,
                    rh=Testr.rh,
                    rhmax=Testr.rhmax,
                    rhmin=Testr.rhmin,
                )
                * Testr.lambda_corr
            )
            .round(1)
            .tolist()
        )
        r_fao56 = [5.1, 5.0, 5.0, 5.2, 5.9, 8.8, 6.3, 7.8, 8.0, 5.7]
        self.assertEqual(r_fao56, pyet_pm, 1)

        expected_pm = tile(asarray(r_fao56)[:, newaxis, newaxis], (1, 3, 3))
        expected_pm = DataArray(
            expected_pm,
            coords={"time": Testr.tmean.index, "y": range(3), "x": range(3)},
            dims=["time", "y", "x"],
        )
        wind_xr = Testr.uz_xr * 4.87 / log(67.8 * 10 - 5.42)
        calculated_pm = (
            et.pm(
                Testr.tmean_xr,
                wind_xr,
                rs=Testr.rs_xr,
                elevation=Testr.elevation,
                lat=Testr.lat,
                tmax=Testr.tmax_xr,
                tmin=Testr.tmin_xr,
                rh=Testr.rh_xr,
                rhmax=Testr.rhmax_xr,
                rhmin=Testr.rhmin_xr,
            )
            * Testr.lambda_corr.values[:, None, None]
        )
        testing.assert_allclose(calculated_pm.round(1), expected_pm.round(1))
        calculated_pm1 = (
            et.pm(
                Testr.tmean_xr,
                wind_xr,
                rs=Testr.rs_xr,
                elevation=Testr.elevation_xr,
                lat=Testr.lat_xr,
                tmax=Testr.tmax_xr,
                tmin=Testr.tmin_xr,
                rh=Testr.rh_xr,
                rhmax=Testr.rhmax_xr,
                rhmin=Testr.rhmin_xr,
            )
            * Testr.lambda_corr.values[:, None, None]
        )
        testing.assert_allclose(calculated_pm1.round(1), expected_pm.round(1))

    def test_asce_pm(self):
        wind = Testr.uz * 4.87 / log(67.8 * 10 - 5.42)
        pyet_pm = (
            et.pm_asce(
                Testr.tmean,
                wind,
                rs=Testr.rs,
                elevation=Testr.elevation,
                lat=Testr.lat,
                tmax=Testr.tmax,
                tmin=Testr.tmin,
                rh=Testr.rh,
                rhmax=Testr.rhmax,
                rhmin=Testr.rhmin,
            )
            .round(1)
            .tolist()
        )
        r_fao56 = [5.1, 5.0, 5.0, 5.2, 5.9, 8.8, 6.3, 7.8, 8.0, 5.7]
        self.assertEqual(r_fao56, pyet_pm, 1)

        expected_pm = tile(asarray(r_fao56)[:, newaxis, newaxis], (1, 3, 3))
        expected_pm = DataArray(
            expected_pm,
            coords={"time": Testr.tmean.index, "y": range(3), "x": range(3)},
            dims=["time", "y", "x"],
        )
        wind_xr = Testr.uz_xr * 4.87 / log(67.8 * 10 - 5.42)
        calculated_pm = et.pm_asce(
            Testr.tmean_xr,
            wind_xr,
            rs=Testr.rs_xr,
            elevation=Testr.elevation,
            lat=Testr.lat,
            tmax=Testr.tmax_xr,
            tmin=Testr.tmin_xr,
            rh=Testr.rh_xr,
            rhmax=Testr.rhmax_xr,
            rhmin=Testr.rhmin_xr,
        )
        testing.assert_allclose(calculated_pm.round(1), expected_pm.round(1))
        calculated_pm1 = et.pm_asce(
            Testr.tmean_xr,
            wind_xr,
            rs=Testr.rs_xr,
            elevation=Testr.elevation_xr,
            lat=Testr.lat_xr,
            tmax=Testr.tmax_xr,
            tmin=Testr.tmin_xr,
            rh=Testr.rh_xr,
            rhmax=Testr.rhmax_xr,
            rhmin=Testr.rhmin_xr,
        )
        testing.assert_allclose(calculated_pm1.round(1), expected_pm.round(1))

    def test_makkink(self):
        pyet_makk = (
            (
                et.makkink(Testr.tmean, rs=Testr.rs, elevation=Testr.elevation, k=0.61)
                * Testr.lambda_corr
                - 0.12
            )
            .round(1)
            .tolist()
        )
        r_makkink = [3.5, 3.4, 3.5, 3.4, 3.6, 3.8, 3.6, 3.7, 3.7, 3.3]
        self.assertEqual(r_makkink, pyet_makk, 1)

        expected_makkink = tile(asarray(r_makkink)[:, newaxis, newaxis], (1, 3, 3))
        expected_makkink = DataArray(
            expected_makkink,
            coords={"time": Testr.tmean.index, "y": range(3), "x": range(3)},
            dims=["time", "y", "x"],
        )
        calculated_makk = (
            et.makkink(
                Testr.tmean_xr, rs=Testr.rs_xr, elevation=Testr.elevation, k=0.61
            )
            * Testr.lambda_corr.values[:, None, None]
            - 0.12
        )
        testing.assert_allclose(calculated_makk.round(1), expected_makkink.round(1))
        calculated_makk1 = (
            et.makkink(
                Testr.tmean_xr, rs=Testr.rs_xr, elevation=Testr.elevation_xr, k=0.61
            )
            * Testr.lambda_corr.values[:, None, None]
            - 0.12
        )
        testing.assert_allclose(calculated_makk1.round(1), expected_makkink.round(1))

    def test_makkink_knmi(self):
        pyet_makk_knmi = et.makkink_knmi(Testr.tmean, rs=Testr.rs).round(1).tolist()
        r_makkink = [3.8, 3.8, 3.9, 3.8, 4.0, 4.3, 4.0, 4.1, 4.1, 3.7]
        self.assertEqual(r_makkink, pyet_makk_knmi, 1)

        expected_makkink = tile(asarray(r_makkink)[:, newaxis, newaxis], (1, 3, 3))
        expected_makkink = DataArray(
            expected_makkink,
            coords={"time": Testr.tmean.index, "y": range(3), "x": range(3)},
            dims=["time", "y", "x"],
        )
        makk_knmi = et.makkink_knmi(Testr.tmean_xr, rs=Testr.rs_xr)
        testing.assert_allclose(makk_knmi.round(1), expected_makkink.round(1))

    def test_pt(self):
        pyet_pt = (
            (
                et.priestley_taylor(
                    Testr.tmean,
                    rs=Testr.rs,
                    lat=Testr.lat,
                    elevation=Testr.elevation,
                    tmax=Testr.tmax,
                    tmin=Testr.tmin,
                    rhmax=Testr.rhmax,
                    rhmin=Testr.rhmin,
                    rh=Testr.rh,
                )
                * Testr.lambda_corr
            )
            .round(1)
            .tolist()
        )
        r_pt = [4.0, 3.9, 4.0, 3.9, 4.2, 4.1, 3.9, 3.9, 3.6, 3.7]
        self.assertEqual(r_pt, pyet_pt, 1)

        expected_pt = tile(asarray(r_pt)[:, newaxis, newaxis], (1, 3, 3))
        expected_pt = DataArray(
            expected_pt,
            coords={"time": Testr.tmean.index, "y": range(3), "x": range(3)},
            dims=["time", "y", "x"],
        )
        calculated_pt = (
            et.priestley_taylor(
                Testr.tmean_xr,
                rs=Testr.rs_xr,
                lat=Testr.lat,
                elevation=Testr.elevation,
                tmax=Testr.tmax_xr,
                tmin=Testr.tmin_xr,
                rhmax=Testr.rhmax_xr,
                rhmin=Testr.rhmin_xr,
                rh=Testr.rh_xr,
            )
            * Testr.lambda_corr.values[:, None, None]
        )
        testing.assert_allclose(calculated_pt.round(1), expected_pt.round(1))
        calculated_pt1 = (
            et.priestley_taylor(
                Testr.tmean_xr,
                rs=Testr.rs_xr,
                lat=Testr.lat_xr,
                elevation=Testr.elevation_xr,
                tmax=Testr.tmax_xr,
                tmin=Testr.tmin_xr,
                rhmax=Testr.rhmax_xr,
                rhmin=Testr.rhmin_xr,
                rh=Testr.rh_xr,
            )
            * Testr.lambda_corr.values[:, None, None]
        )
        testing.assert_allclose(calculated_pt1.round(1), expected_pt.round(1))

    def test_hargreaves(self):
        pyet_har = (
            (et.hargreaves(Testr.tmean, Testr.tmax, Testr.tmin, Testr.lat, method=1))
            .round(1)
            .tolist()
        )
        r_har = [4.6, 4.3, 4.3, 3.7, 5.4, 4.6, 4.8, 4.4, 4.9, 3.7]
        self.assertEqual(r_har, pyet_har, 1)

        expected_har = tile(asarray(r_har)[:, newaxis, newaxis], (1, 3, 3))
        expected_har = DataArray(
            expected_har,
            coords={"time": Testr.tmean.index, "y": range(3), "x": range(3)},
            dims=["time", "y", "x"],
        )
        calculated_har = et.hargreaves(
            Testr.tmean_xr, Testr.tmax_xr, Testr.tmin_xr, Testr.lat, method=1
        )
        testing.assert_allclose(calculated_har.round(1), expected_har.round(1))
        calculated_har1 = et.hargreaves(
            Testr.tmean_xr, Testr.tmax_xr, Testr.tmin_xr, Testr.lat_xr, method=1
        )
        testing.assert_allclose(calculated_har1.round(1), expected_har.round(1))

    def test_hamon(self):
        pyet_pm = (
            et.hamon(
                Testr.tmean,
                lat=Testr.lat,
                n=Testr.n,
                tmax=Testr.tmax,
                tmin=Testr.tmin,
                method=3,
            )
            .round(1)
            .tolist()
        )
        r_hamon = [1.5, 1.4, 1.5, 1.4, 1.8, 2.4, 2.0, 2.1, 2.2, 1.5]
        self.assertEqual(r_hamon, pyet_pm, 1)

        expected_ham = tile(asarray(r_hamon)[:, newaxis, newaxis], (1, 3, 3))
        expected_ham = DataArray(
            expected_ham,
            coords={"time": Testr.tmean.index, "y": range(3), "x": range(3)},
            dims=["time", "y", "x"],
        )
        calculated_ham = et.hamon(
            Testr.tmean_xr,
            lat=Testr.lat,
            n=Testr.n_xr,
            tmax=Testr.tmax_xr,
            tmin=Testr.tmin_xr,
            method=3,
        )
        testing.assert_allclose(calculated_ham.round(1), expected_ham.round(1))
        calculated_ham1 = et.hamon(
            Testr.tmean_xr,
            lat=Testr.lat_xr,
            n=Testr.n_xr,
            tmax=Testr.tmax_xr,
            tmin=Testr.tmin_xr,
            method=3,
        )
        testing.assert_allclose(calculated_ham1.round(1), expected_ham.round(1))

    def test_blaney_criddle(self):
        wind = Testr.uz * 4.87 / log(67.8 * 10 - 5.42)
        pyet_bc = (
            (
                et.blaney_criddle(
                    Testr.tmean,
                    Testr.lat,
                    rhmin=Testr.rhmin,
                    wind=wind,
                    method=2,
                    n=Testr.n,
                    clip_zero=False,
                )
            )
            .round(1)
            .tolist()
        )
        r_bc = [3.0, 3.0, 3.1, 2.9, 3.6, 5.0, 4.1, 4.5, 4.9, 3.3]
        self.assertEqual(r_bc, pyet_bc, 1)

        expected_bc = tile(asarray(r_bc)[:, newaxis, newaxis], (1, 3, 3))
        expected_bc = DataArray(
            expected_bc,
            coords={"time": Testr.tmean.index, "y": range(3), "x": range(3)},
            dims=["time", "y", "x"],
        )
        wind_xr = Testr.uz_xr * 4.87 / log(67.8 * 10 - 5.42)
        calculated_bc = et.blaney_criddle(
            Testr.tmean_xr,
            Testr.lat,
            rhmin=Testr.rhmin_xr,
            wind=wind_xr,
            method=2,
            n=Testr.n_xr,
            clip_zero=False,
        )
        testing.assert_allclose(calculated_bc.round(1), expected_bc.round(1))
        calculated_bc1 = et.blaney_criddle(
            Testr.tmean_xr,
            Testr.lat_xr,
            rhmin=Testr.rhmin_xr,
            wind=wind_xr,
            method=2,
            n=Testr.n_xr,
            clip_zero=False,
        )
        testing.assert_allclose(calculated_bc1.round(1), expected_bc.round(1))

    def test_romanenko(self):
        pyet_romanenko = (
            (
                et.romanenko(
                    Testr.tmean,
                    rh=Testr.rh,
                    tmax=Testr.tmax,
                    tmin=Testr.tmin,
                    rhmax=Testr.rhmax,
                    rhmin=Testr.rhmin,
                )
            )
            .round(1)
            .tolist()
        )
        r_romanenko = [9.3, 8.9, 9.4, 8.2, 10.6, 16.8, 14.0, 14.7, 17.4, 9.2]
        self.assertEqual(r_romanenko, pyet_romanenko, 1)

        expected_rm = tile(asarray(r_romanenko)[:, newaxis, newaxis], (1, 3, 3))
        expected_rm = DataArray(
            expected_rm,
            coords={"time": Testr.tmean.index, "y": range(3), "x": range(3)},
            dims=["time", "y", "x"],
        )
        calculated_rm = et.romanenko(
            Testr.tmean_xr,
            rh=Testr.rh_xr,
            tmax=Testr.tmax_xr,
            tmin=Testr.tmin_xr,
            rhmax=Testr.rhmax_xr,
            rhmin=Testr.rhmin_xr,
        )
        testing.assert_allclose(calculated_rm.round(1), expected_rm.round(1))

    def test_mcguinnessbordne(self):
        pyet_mgb = (
            (et.mcguinness_bordne(Testr.tmean, lat=Testr.lat) * Testr.lambda_corr)
            .round(1)
            .tolist()
        )
        r_mgb = [5.8, 5.5, 5.9, 5.6, 6.4, 7.6, 6.7, 6.8, 7.0, 5.6]
        self.assertEqual(r_mgb, pyet_mgb, 1)

        expected_mb = tile(asarray(r_mgb)[:, newaxis, newaxis], (1, 3, 3))
        expected_mb = DataArray(
            expected_mb,
            coords={"time": Testr.tmean.index, "y": range(3), "x": range(3)},
            dims=["time", "y", "x"],
        )
        calc_mb = (
            et.mcguinness_bordne(Testr.tmean_xr, lat=Testr.lat)
            * Testr.lambda_corr.values[:, None, None]
        )
        testing.assert_allclose(calc_mb.round(1), expected_mb.round(1))
        calc_mb1 = (
            et.mcguinness_bordne(Testr.tmean_xr, lat=Testr.lat_xr)
            * Testr.lambda_corr.values[:, None, None]
        )
        testing.assert_allclose(calc_mb1.round(1), expected_mb.round(1))

    def test_jensen_haise(self):
        pyet_jh = (
            (et.jensen_haise(Testr.tmean, Testr.rs) * Testr.lambda_corr)
            .round(1)
            .tolist()
        )
        r_jh = [5.2, 4.9, 5.3, 5.0, 5.8, 7.0, 6.1, 6.3, 6.4, 5.0]
        self.assertEqual(r_jh, pyet_jh, 1)

        expected_jh = tile(asarray(r_jh)[:, newaxis, newaxis], (1, 3, 3))
        expected_jh = DataArray(
            expected_jh,
            coords={"time": Testr.tmean.index, "y": range(3), "x": range(3)},
            dims=["time", "y", "x"],
        )
        calc_jh = (
            et.jensen_haise(Testr.tmean_xr, Testr.rs_xr)
            * Testr.lambda_corr.values[:, None, None]
        )
        testing.assert_allclose(calc_jh.round(1), expected_jh.round(1))

    def test_linacre(self):
        tdew = [
            10.2375,
            8.8375,
            11.5125,
            10.2875,
            11.9875,
            9.475,
            7.7625,
            7.2375,
            3.0125,
            11.325,
        ]
        pyet_linacre = (
            (et.linacre(Testr.tmean, Testr.elevation, Testr.lat, tdew=tdew))
            .round(1)
            .tolist()
        )
        r_linacre = [4.4, 4.3, 4.4, 4.2, 5.4, 9.1, 7.5, 8.0, 9.9, 4.3]
        self.assertEqual(r_linacre, pyet_linacre, 1)

        expected_l = tile(asarray(r_linacre)[:, newaxis, newaxis], (1, 3, 3))
        expected_l = DataArray(
            expected_l,
            coords={"time": Testr.tmean.index, "y": range(3), "x": range(3)},
            dims=["time", "y", "x"],
        )
        tdew_xr = tile(asarray(tdew)[:, newaxis, newaxis], (1, 3, 3))
        tdew_xr = DataArray(
            tdew_xr,
            coords={"time": Testr.tmean.index, "y": range(3), "x": range(3)},
            dims=["time", "y", "x"],
        )
        calc_l = et.linacre(Testr.tmean_xr, Testr.elevation, Testr.lat, tdew=tdew_xr)
        testing.assert_allclose(calc_l.round(1), expected_l.round(1))
        calc_l1 = et.linacre(
            Testr.tmean_xr, Testr.elevation_xr, Testr.lat_xr, tdew=tdew_xr
        )
        testing.assert_allclose(calc_l1.round(1), expected_l.round(1))

    def test_abtew(self):
        pyet_abtew = (
            (et.abtew(Testr.tmean, Testr.rs, k=0.52) * Testr.lambda_corr)
            .round(1)
            .tolist()
        )
        r_abtew = [4.3, 4.3, 4.3, 4.3, 4.3, 4.2, 4.2, 4.2, 4.2, 4.1]
        self.assertEqual(r_abtew, pyet_abtew, 1)

        expected_a = tile(asarray(r_abtew)[:, newaxis, newaxis], (1, 3, 3))
        expected_a = DataArray(
            expected_a,
            coords={"time": Testr.tmean.index, "y": range(3), "x": range(3)},
            dims=["time", "y", "x"],
        )
        calc_a = (
            et.abtew(Testr.tmean_xr, Testr.rs_xr, k=0.52)
            * Testr.lambda_corr.values[:, None, None]
        )
        testing.assert_allclose(calc_a.round(1), expected_a.round(1))

    def test_turc(self):
        pyet_turc = (et.turc(Testr.tmean, Testr.rs, Testr.rh)).round(1).tolist()
        # There is a mistake in the R-Package
        # The R-Package does not consider RH < 50
        # if pyet.turc also does nonsider this, the result is the same
        r_turc = [4.2, 4.0, 4.2, 4.0, 4.3, 6.1, 5.4, 5.6, 6.2, 4.1]
        self.assertEqual(r_turc, pyet_turc, 1)

        expected_t = tile(asarray(r_turc)[:, newaxis, newaxis], (1, 3, 3))
        expected_t = DataArray(
            expected_t,
            coords={"time": Testr.tmean.index, "y": range(3), "x": range(3)},
            dims=["time", "y", "x"],
        )
        calc_turc = et.turc(Testr.tmean_xr, Testr.rs_xr, Testr.rh_xr)
        testing.assert_allclose(calc_turc.round(1), expected_t.round(1))
