import pytest
import datetime
import calendar
import numpy as np
from vic.vic import ffi
from vic import lib as vic_lib

try:
    from netcdftime import netcdftime
    from netCDF4 import date2num as nc_date2num
    nctime_unavailable = False
    nc_reason = 'installed'
except ImportError:
    nctime_unavailable = True
    nc_reason = 'netCDF4.netcdftime not available'


calendars = {'standard': vic_lib.CALENDAR_STANDARD,
             'gregorian': vic_lib.CALENDAR_GREGORIAN,
             'proleptic_gregorian': vic_lib.CALENDAR_PROLEPTIC_GREGORIAN,
             'noleap': vic_lib.CALENDAR_NOLEAP,
             'julian': vic_lib.CALENDAR_JULIAN,
             'all_leap': vic_lib.CALENDAR_ALL_LEAP,
             '365_day': vic_lib.CALENDAR_365_DAY,
             '366_day': vic_lib.CALENDAR_366_DAY,
             '360_day': vic_lib.CALENDAR_360_DAY}

units = {'days': vic_lib.TIME_UNITS_DAYS,
         'hours': vic_lib.TIME_UNITS_HOURS,
         'minutes': vic_lib.TIME_UNITS_MINUTES,
         'seconds': vic_lib.TIME_UNITS_SECONDS}


def datetime_to_dmy(dt):
    dmy = ffi.new('dmy_struct *')
    dmy[0].year = dt.year
    dmy[0].month = dt.month
    dmy[0].day = dt.day
    dmy[0].dayseconds = int(dt.hour * 3600 + dt.minute * 60 + dt.second +
                            dt.microsecond / 100000)
    dmy[0].day_in_year = dt.timetuple().tm_yday
    return dmy


@pytest.fixture()
def feb3_noon(scope='module'):
    return datetime.datetime(2015, 2, 3, 12)


@pytest.fixture()
def dmy_feb_3_noon(feb3_noon, scope='module'):
    return datetime_to_dmy(feb3_noon)


@pytest.fixture()
def dmy_june31(scope='module'):
    dmy = ffi.new('dmy_struct *')
    dmy[0].year = 1984
    dmy[0].month = 6
    dmy[0].day = 31
    dmy[0].dayseconds = 0
    dmy[0].day_in_year = 180
    return dmy


@pytest.fixture()
def dmy_now(scope='function'):
    now = datetime.datetime.now()
    return datetime_to_dmy(now)


def test_fractional_day_from_dmy(dmy_feb_3_noon):
    vic_lib.print_dmy(dmy_feb_3_noon)
    frac_day = vic_lib.fractional_day_from_dmy(dmy_feb_3_noon)
    assert frac_day == 3.5


@pytest.mark.skipif(nctime_unavailable, reason=nc_reason)
def test_julian_day_from_dmy(feb3_noon, dmy_feb_3_noon):
    for cal in ['standard', 'gregorian', 'proleptic_gregorian', 'julian']:
        d_vic = vic_lib.julian_day_from_dmy(dmy_feb_3_noon,
                                            calendars[cal])
        d_nc = netcdftime.JulianDayFromDate(feb3_noon, cal)
        np.testing.assert_approx_equal(d_vic, d_nc)


@pytest.mark.skipif(nctime_unavailable, reason=nc_reason)
def test_no_leap_day_from_dmy(feb3_noon, dmy_feb_3_noon):
    d_vic = vic_lib.no_leap_day_from_dmy(dmy_feb_3_noon)
    d_nc = netcdftime._NoLeapDayFromDate(feb3_noon)
    np.testing.assert_approx_equal(d_vic, d_nc)


@pytest.mark.skipif(nctime_unavailable, reason=nc_reason)
def test_all_leap_from_dmy(feb3_noon, dmy_feb_3_noon):
    d_vic = vic_lib.all_leap_from_dmy(dmy_feb_3_noon)
    d_nc = netcdftime._AllLeapFromDate(feb3_noon)
    np.testing.assert_approx_equal(d_vic, d_nc)


@pytest.mark.skipif(nctime_unavailable, reason=nc_reason)
def test_all_30_day_from_dmy(feb3_noon, dmy_feb_3_noon):
    d_vic = vic_lib.all_30_day_from_dmy(dmy_feb_3_noon)
    d_nc = netcdftime._360DayFromDate(feb3_noon)
    np.testing.assert_approx_equal(d_vic, d_nc)


@pytest.mark.skipif(nctime_unavailable, reason=nc_reason)
def test_date2num():
    pass


@pytest.mark.skipif(nctime_unavailable, reason=nc_reason)
def test_dmy_julian_day():
    pass


@pytest.mark.skipif(nctime_unavailable, reason=nc_reason)
def test_dmy_no_leap_day():
    pass


@pytest.mark.skipif(nctime_unavailable, reason=nc_reason)
def test_dmy_all_leap():
    pass


@pytest.mark.skipif(nctime_unavailable, reason=nc_reason)
def test_dmy_all_30_day():
    pass


@pytest.mark.skipif(nctime_unavailable, reason=nc_reason)
def test_num2date():
    pass


@pytest.mark.xfail
@pytest.mark.skipif(nctime_unavailable, reason=nc_reason)
def test_make_lastday():
    lastday = np.zeros(12)
    for year in np.arange(1900, 2100):
        for cal in ['standard', 'gregorian', 'proleptic_gregorian', 'julian']:
            vic_lib.make_lastday(calendars[cal], year,
                                 ffi.addressof(lastday.ctypes.data))
            cal_dpm = [calendar.monthrange(year, month)[1] for month in
                       np.arange(1, 13)]
            np.testing.assert_equal(cal_dpm, lastday)
        for cal in ['noleap', '365_day']:
            vic_lib.make_lastday(calendars[cal], year,
                                 ffi.addressof(lastday.ctypes.data))
            cal_dpm = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
            np.testing.assert_equal(cal_dpm, lastday)
        for cal in ['all_leap', '366_day']:
            vic_lib.make_lastday(calendars[cal], year,
                                 ffi.addressof(lastday.ctypes.data))
            cal_dpm = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
            np.testing.assert_equal(cal_dpm, lastday)
        for cal in ['all_leap', '366_day']:
            vic_lib.make_lastday(calendars[cal], year,
                                 ffi.addressof(lastday.ctypes.data))
            cal_dpm = [30] * 12
            np.testing.assert_equal(cal_dpm, lastday)


@pytest.mark.xfail
@pytest.mark.skipif(nctime_unavailable, reason=nc_reason)
def test_initialize_time():
    assert vic_lib.initialize_time() is None
    assert vic_lib.global_param.time_origin_num != -99999
    np.testing.assert_almost_equal(
        vic_lib.global_param.time_origin_num,
        nc_date2num(datetime.datetime(1, 1, 1, 0, 0), 'days since 0001-01-01'))


def test_valid_date(dmy_feb_3_noon, dmy_june31):
    for cal in calendars:
        assert vic_lib.valid_date(calendars[cal], dmy_feb_3_noon) == 0
        assert vic_lib.valid_date(calendars[cal], dmy_june31) > 0


@pytest.mark.skipif(nctime_unavailable, reason=nc_reason)
def test_dt_seconds_to_time_units():
    dt_time_units = ffi.new('double *')
    for tu in units:
        dt_time_units[0] = 0.
        vic_lib.dt_seconds_to_time_units(units[tu], 3600, dt_time_units)
        assert dt_time_units[0] > 0.
