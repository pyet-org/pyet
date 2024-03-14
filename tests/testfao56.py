import unittest
from numpy import pi, array, full, testing
from pandas import date_range, DatetimeIndex, Series
from xarray import DataArray, Dataset
import pyet as et


class TestFAO56(unittest.TestCase):
    def test_press_calc(self):
        # Based on Example 2 and 4, p. 32 FAO.
        elevation_s = Series(
            data=[1800, 1200, 1462.4], index=date_range(start="2020-1-1", periods=3)
        )
        calculated_pressures0 = et.calc_press(elevation_s)
        expected_pressures = Series(
            data=[81.8, 87.9, 85.17], index=date_range(start="2020-1-1", periods=3)
        )
        testing.assert_allclose(
            expected_pressures.round(1), calculated_pressures0.round(1)
        )

        # Create a 3D DataArray with dimensions 'x', 'y', and 'time'
        elevation_xr = DataArray(
            full((2, 3, 3), 1800),
            coords=[
                ("time", date_range(start="1/1/2020", periods=2)),
                ("y", [1, 2, 3]),
                ("x", [1, 2, 3]),
            ],
        )
        elevation_xr.loc[{"time": "2020-01-02"}] = 1200
        # Apply the et.calc_press function to each element of the DataArray
        calculated_pressures = et.calc_press(elevation_xr)
        # Create expected results as DataArrays
        expected_pressures = DataArray(
            full((2, 3, 3), 81.8),
            coords=[
                ("time", date_range(start="1/1/2020", periods=2)),
                ("y", [1, 2, 3]),
                ("x", [1, 2, 3]),
            ],
        )
        expected_pressures.loc[{"time": "2020-01-02"}] = 87.9
        # Check that the results are as expected
        testing.assert_allclose(calculated_pressures.round(1), expected_pressures)

    def test_vpc(self):
        # Based on ASCE Table C-3
        tmean = [21.65, 22.9, 23.7, 22.8, 24.3, 26.0, 26.1, 26.4, 23.9, 24.2]
        tmean_s = Series(tmean, index=date_range(start="2020-1-1", periods=10))
        calculated_vpc1 = et.calc_vpc(tmean_s)
        vpc1 = [
            0.1585,
            0.1692,
            0.1762,
            0.1684,
            0.1820,
            0.199,
            0.1996,
            0.2027,
            0.1781,
            0.1809,
        ]
        expected_vpc1 = Series(vpc1, index=date_range(start="2020-1-1", periods=10))
        testing.assert_allclose(expected_vpc1.round(2), calculated_vpc1.round(2))

        # Create a 3D DataArray with dimensions 'x', 'y', and 'time'
        tmean_xr = DataArray(
            full((4, 3, 3), 21.65),
            coords=[
                ("time", date_range(start="2020-1-1", periods=4)),
                ("y", [1, 2, 3]),
                ("x", [1, 2, 3]),
            ],
        )
        tmean_xr.loc[{"time": ["2020-01-02", "2020-01-03", "2020-01-04"]}] = [
            22.9,
            23.7,
            22.8,
        ]
        # Apply the et.calc_press function to each element of the DataArray
        calculated_vpc = et.calc_vpc(tmean_xr)
        # Create expected results as DataArrays
        expected_vpc = DataArray(
            full((4, 3, 3), 0.158),
            coords=[
                ("time", date_range(start="2020-1-1", periods=4)),
                ("y", [1, 2, 3]),
                ("x", [1, 2, 3]),
            ],
        )
        expected_vpc.loc[{"time": ["2020-01-02", "2020-01-03", "2020-01-04"]}] = [
            0.1692,
            0.1762,
            0.1684,
        ]
        # Check that the results are as expected
        testing.assert_allclose(calculated_vpc.round(3), expected_vpc.round(3))

    def test_psy_calc(self):
        # Based on Example 2, p. 32 FAO. and Table C-2 in ASCE(2001)
        p0 = Series([81.8, 85.17], index=date_range(start="2020-1-1", periods=2))
        calculated_psy = et.calc_psy(p0)
        expected_psy = Series(
            [0.054, 0.0566], index=date_range(start="2020-1-1", periods=2)
        )
        testing.assert_allclose(expected_psy.round(3), calculated_psy.round(3))
        # Create a 3D DataArray with dimensions 'x', 'y', and 'time'
        press_xr = DataArray(
            full((2, 3, 3), 81.8),
            coords=[
                ("time", date_range(start="2020-1-1", periods=2)),
                ("y", [1, 2, 3]),
                ("x", [1, 2, 3]),
            ],
        )
        press_xr.loc[{"time": "2020-01-02"}] = 85.17
        # Apply the et.calc_press function to each element of the DataArray
        calculated_psy = et.calc_psy(press_xr)
        # Create expected results as DataArrays
        expected_psy = DataArray(
            full((2, 3, 3), 0.054),
            coords=[
                ("time", date_range(start="2020-1-1", periods=2)),
                ("y", [1, 2, 3]),
                ("x", [1, 2, 3]),
            ],
        )
        expected_psy.loc[{"time": "2020-01-02"}] = 0.0566
        # Check that the results are as expected
        testing.assert_allclose(calculated_psy.round(3), expected_psy.round(3))

    def test_e0_calc(self):
        # Based on Example 3 and 4, p. 36 FAO.
        tmean0 = Series(
            [24.5, 15, 19.75, 19.5], index=date_range(start="2020-1-1", periods=4)
        )
        calculated_e0 = et.calc_e0(tmean0)
        expected_e0 = Series(
            [3.075, 1.705, 2.3, 2.267], index=date_range(start="2020-1-1", periods=4)
        )
        testing.assert_allclose(expected_e0.round(1), calculated_e0.round(1))

        # Create a 3D DataArray with dimensions 'x', 'y', and 'time'
        tmean_xr = DataArray(
            full((4, 3, 3), 24.5),
            coords=[
                ("time", date_range(start="2020-1-1", periods=4)),
                ("y", [1, 2, 3]),
                ("x", [1, 2, 3]),
            ],
        )
        tmean_xr.loc[{"time": "2020-01-02"}] = 15
        tmean_xr.loc[{"time": "2020-01-03"}] = 19.75
        tmean_xr.loc[{"time": "2020-01-04"}] = 19.5
        # Apply the et.calc_press function to each element of the DataArray
        calculated_e0 = et.calc_e0(tmean_xr)
        # Create expected results as DataArrays
        expected_e0 = DataArray(
            full((4, 3, 3), 3.075),
            coords=[
                ("time", date_range(start="2020-1-1", periods=4)),
                ("y", [1, 2, 3]),
                ("x", [1, 2, 3]),
            ],
        )
        expected_e0.loc[{"time": "2020-01-02"}] = 1.705
        expected_e0.loc[{"time": "2020-01-03"}] = 2.3
        expected_e0.loc[{"time": "2020-01-04"}] = 2.267
        # Check that the results are as expected
        testing.assert_allclose(calculated_e0.round(1), expected_e0.round(1))

    def test_es_calc(self):
        # Based on Example 3 and 6, p. 36 FAO.
        es = et.calc_es(tmax=24.5, tmin=15.0)
        self.assertAlmostEqual(es, 2.39, 3)

    def test_ea_calc(self):
        # Based on Example 5, p. 39 FAO.
        ea1 = et.calc_ea(tmax=25.0, tmin=18.0, rhmax=82.0, rhmin=54.0)
        rhmean = (82 + 54) / 2
        ea2 = et.calc_ea(tmax=25.0, tmin=18.0, rh=rhmean)
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
        doy = et.day_of_year(DatetimeIndex(["2015-09-03"]))
        self.assertAlmostEqual(float(doy.iloc[0]), 246, 1)
        # Based on ASCD Table C-3
        dindex = date_range("2020-07-01", "2020-07-10")
        doy1 = et.day_of_year(dindex).tolist()
        doyr = [183, 184, 185, 186, 187, 188, 189, 190, 191, 192]
        self.assertEqual(doy1, doyr, 1)

    def test_extraterrestrial_r(self):
        # Based on Example 8, p. 47 FAO.
        extrar = et.extraterrestrial_r(DatetimeIndex(["2015-09-03"]), -0.35)
        self.assertAlmostEqual(float(extrar), 32.2, 1)

    def test_daylight_hours(self):
        # Based on Example 9, p. 47 FAO.
        dayhours = et.daylight_hours(DatetimeIndex(["2015-09-03"]), -0.35)
        self.assertAlmostEqual(float(dayhours), 11.7, 1)

    def test_calc_rad_long(self):
        # Based on Example 10, p. 52 FAO.
        rs = Series([14.5], index=DatetimeIndex(["2015-05-15"]))
        tmax = Series([25.1], index=DatetimeIndex(["2015-05-15"]))
        tmin = Series([19], index=DatetimeIndex(["2015-05-15"]))
        ea = Series([2.1], index=DatetimeIndex(["2015-05-15"]))
        rso = Series([18.8], index=DatetimeIndex(["2015-05-15"]))
        rnl = et.calc_rad_long(rs, tmax=tmax, tmin=tmin, ea=ea, rso=rso)
        self.assertAlmostEqual(float(rnl), 3.5, 1)

    def test_calc_rad_sol_in(self):
        # Based on example 10, p 50 TestFAO56
        lat = -22.9 * pi / 180
        n = Series(7.1, DatetimeIndex(["2021-5-15"]))
        rad_sol_in = et.calc_rad_sol_in(n, lat)
        self.assertAlmostEqual(float(rad_sol_in.iloc[0]), 14.5, 1)

    def test_et_fao56(self):
        # Based on Example 18, p. 72 FAO.
        wind = Series([2.078], index=DatetimeIndex(["2015-07-06"]))
        tmax = Series([21.5], index=DatetimeIndex(["2015-07-06"]))
        tmin = Series([12.3], index=DatetimeIndex(["2015-07-06"]))
        tmean = (tmax + tmin) / 2
        rhmax = Series([84], index=DatetimeIndex(["2015-07-06"]))
        rhmin = Series([63], index=DatetimeIndex(["2015-07-06"]))
        rs = Series([22.07], index=DatetimeIndex(["2015-07-06"]))
        n = 9.25
        nn = 16.1
        elevation = 100
        lat = 50.80 * pi / 180
        et56 = et.pm_fao56(
            tmean,
            wind,
            elevation=elevation,
            lat=lat,
            rs=rs,
            tmax=tmax,
            tmin=tmin,
            rhmax=rhmax,
            rhmin=rhmin,
            n=n,
            nn=nn,
        )
        self.assertAlmostEqual(float(et56.iloc[0]), 3.9, 1)
        # Create an xarray Dataset with DataArrays
        # Create an xarray Dataset with DataArrays
        data = {
            "wind": DataArray(wind, dims=("time",)),
            "tmax": DataArray(tmax, dims=("time",)),
            "tmin": DataArray(tmin, dims=("time",)),
            "tmean": DataArray(tmean, dims=("time",)),
            "rhmax": DataArray(rhmax, dims=("time",)),
            "rhmin": DataArray(rhmin, dims=("time",)),
            "rs": DataArray(rs, dims=("time",)),
            "et0": DataArray([3.9], dims=("time",)),
        }

        # Create the xarray Dataset
        dataset = Dataset(data)
        # Add additional variables (n, nn, elevation, lat)
        dataset["n"] = n
        dataset["nn"] = nn
        dataset["elevation"] = elevation
        dataset["lat"] = lat
        et56_xr = et.pm_fao56(
            dataset["tmean"],
            dataset["wind"],
            elevation=dataset["elevation"],
            lat=dataset["lat"],
            rs=dataset["rs"],
            tmax=dataset["tmax"],
            tmin=dataset["tmin"],
            rhmax=dataset["rhmax"],
            rhmin=dataset["rhmin"],
            n=dataset["n"],
            nn=dataset["nn"],
        )

        testing.assert_allclose(et56_xr.round(1), dataset["et0"])
