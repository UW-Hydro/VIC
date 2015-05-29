import pytest
import datetime
import calendar
import ctypes
import numpy as np
from vic.vic import (fractional_day_from_dmy, julian_day_from_dmy,
                     no_leap_day_from_dmy, all_leap_from_dmy,
                     all_30_day_from_dmy, make_lastday,
                     initialize_time, valid_date, dt_seconds_to_time_units,
                     dmy_struct, print_dmy,
                     global_param, MISSING,
                     CALENDAR_360_DAY, CALENDAR_ALL_LEAP, CALENDAR_NOLEAP,
                     CALENDAR_365_DAY, CALENDAR_GREGORIAN,
                     CALENDAR_PROLEPTIC_GREGORIAN, CALENDAR_366_DAY,
                     CALENDAR_JULIAN, CALENDAR_STANDARD, TIME_UNITS_DAYS,
                     TIME_UNITS_HOURS, TIME_UNITS_MINUTES, TIME_UNITS_SECONDS)

try:
    from netcdftime import netcdftime
    from netCDF4 import date2num as nc_date2num
    nctime_unavailable = False
except ImportError:
    nctime_unavailable = True


calendars = {'standard': CALENDAR_STANDARD,
             'gregorian': CALENDAR_GREGORIAN,
             'proleptic_gregorian': CALENDAR_PROLEPTIC_GREGORIAN,
             'noleap': CALENDAR_NOLEAP,
             'julian': CALENDAR_JULIAN,
             'all_leap': CALENDAR_ALL_LEAP,
             '365_day': CALENDAR_365_DAY,
             '366_day': CALENDAR_366_DAY,
             '360_day': CALENDAR_360_DAY}

units = {'days': TIME_UNITS_DAYS,
         'hours': TIME_UNITS_HOURS,
         'minutes': TIME_UNITS_MINUTES,
         'seconds': TIME_UNITS_SECONDS}


def datetime_to_dmy(dt):
    dmy = dmy_struct()
    dmy.year = dt.year
    dmy.month = dt.month
    dmy.day = dt.day
    dmy.day_seconds = (dt.hour * 3600 + dt.minute * 60 + dt.second +
                       dt.microsecond / 100000)
    dmy.day_in_year = dt.timetuple().tm_yday
    return dmy


# # def test_netcdf4():
#     # This test should be skipped
#     from netCDF4 import netcdftime
#     print(netcdftime)


@pytest.fixture()
def feb3_noon(scope='module'):
    return datetime.datetime(2015, 2, 3, 12)


@pytest.fixture()
def dmy_feb_3_noon(feb3_noon, scope='module'):
    return datetime_to_dmy(feb3_noon)


@pytest.fixture()
def dmy_june31(scope='module'):
    dmy = dmy_struct()
    dmy.year = 1984
    dmy.month = 6
    dmy.day = 31
    dmy.day_seconds = 0
    dmy.day_in_year = 180
    return dmy


@pytest.fixture()
def dmy_now(scope='function'):
    now = datetime.datetime.now()
    return datetime_to_dmy(now)


def test_fractional_day_from_dmy(dmy_feb_3_noon):
    print_dmy(ctypes.byref(dmy_feb_3_noon))
    frac_day = fractional_day_from_dmy(ctypes.byref(dmy_feb_3_noon))
    assert frac_day == 3.5


@pytest.mark.skipif(nctime_unavailable,
                    reason='netCDF4.netcdftime not available')
def test_julian_day_from_dmy(feb3_noon, dmy_feb_3_noon):
    for cal in ['standard', 'gregorian', 'proleptic_gregorian', 'julian']:
        d_vic = julian_day_from_dmy(ctypes.byref(dmy_feb_3_noon),
                                    calendars[cal])
        d_nc = netcdftime.JulianDayFromDate(feb3_noon, cal)
        np.testing.assert_approx_equal(d_vic, d_nc)


@pytest.mark.skipif(nctime_unavailable,
                    reason='netCDF4.netcdftime not available')
def test_no_leap_day_from_dmy(feb3_noon, dmy_feb_3_noon):
    d_vic = no_leap_day_from_dmy(ctypes.byref(dmy_feb_3_noon))
    d_nc = netcdftime._NoLeapDayFromDate(feb3_noon)
    np.testing.assert_approx_equal(d_vic, d_nc)


@pytest.mark.skipif(nctime_unavailable,
                    reason='netCDF4.netcdftime not available')
def test_all_leap_from_dmy(feb3_noon, dmy_feb_3_noon):
    d_vic = all_leap_from_dmy(ctypes.byref(dmy_feb_3_noon))
    d_nc = netcdftime._AllLeapFromDate(feb3_noon)
    np.testing.assert_approx_equal(d_vic, d_nc)


@pytest.mark.skipif(nctime_unavailable,
                    reason='netCDF4.netcdftime not available')
def test_all_30_day_from_dmy(feb3_noon, dmy_feb_3_noon):
    d_vic = all_30_day_from_dmy(ctypes.byref(dmy_feb_3_noon))
    d_nc = netcdftime._360DayFromDate(feb3_noon)
    np.testing.assert_approx_equal(d_vic, d_nc)


@pytest.mark.skipif(nctime_unavailable,
                    reason='netCDF4.netcdftime not available')
def test_date2num():
    pass


@pytest.mark.skipif(nctime_unavailable,
                    reason='netCDF4.netcdftime not available')
def test_dmy_julian_day():
    pass


@pytest.mark.skipif(nctime_unavailable,
                    reason='netCDF4.netcdftime not available')
def test_dmy_no_leap_day():
    pass


@pytest.mark.skipif(nctime_unavailable,
                    reason='netCDF4.netcdftime not available')
def test_dmy_all_leap():
    pass


@pytest.mark.skipif(nctime_unavailable,
                    reason='netCDF4.netcdftime not available')
def test_dmy_all_30_day():
    pass


@pytest.mark.skipif(nctime_unavailable,
                    reason='netCDF4.netcdftime not available')
def test_num2date():
    pass


@pytest.mark.skipif(nctime_unavailable,
                    reason='netCDF4.netcdftime not available')
def test_make_lastday():
    lastday = (ctypes.c_uint * 12)()
    for year in np.arange(1900, 2100):
        for cal in ['standard', 'gregorian', 'proleptic_gregorian', 'julian']:
            make_lastday(calendars[cal], year, ctypes.byref(lastday))
            cal_dpm = [calendar.monthrange(year, month)[1] for month in
                       np.arange(1, 13)]
            np.testing.assert_equal(cal_dpm, lastday)
        for cal in ['noleap', '365_day']:
            make_lastday(calendars[cal], year, ctypes.byref(lastday))
            cal_dpm = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
            np.testing.assert_equal(cal_dpm, lastday)
        for cal in ['all_leap', '366_day']:
            make_lastday(calendars[cal], year, ctypes.byref(lastday))
            cal_dpm = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
            np.testing.assert_equal(cal_dpm, lastday)
        for cal in ['all_leap', '366_day']:
            make_lastday(calendars[cal], year, ctypes.byref(lastday))
            cal_dpm = [30] * 12
            np.testing.assert_equal(cal_dpm, lastday)


@pytest.mark.skipif(nctime_unavailable,
                    reason='netCDF4.netcdftime not available')
def test_initialize_time():
    assert initialize_time() is None
    assert global_param.time_origin_num != MISSING
    np.testing.assert_almost_equal(
        global_param.time_origin_num,
        nc_date2num(datetime.datetime(1, 1, 1, 0, 0), 'days since 0001-01-01'))


def test_valid_date(dmy_feb_3_noon, dmy_june31):
    for cal in calendars:
        assert valid_date(calendars[cal], ctypes.byref(dmy_feb_3_noon)) == 0
        assert valid_date(calendars[cal], ctypes.byref(dmy_june31)) > 0


@pytest.mark.skipif(nctime_unavailable,
                    reason='netCDF4.netcdftime not available')
def test_dt_seconds_to_time_units():
    dt_time_units = ctypes.c_double()
    for tu in units:
        dt_time_units.value = 0.
        dt_seconds_to_time_units(units[tu], 3600, ctypes.byref(dt_time_units))
        assert dt_time_units > 0.
