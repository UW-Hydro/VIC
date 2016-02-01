import pytest
import datetime
import calendar
import numpy as np
import pandas as pd
from vic.vic import ffi
from vic import lib as vic_lib

try:
    from netcdftime import (netcdftime, utime, JulianDayFromDate)
    from netcdftime.netcdftime import (_NoLeapDayFromDate,
                                       _AllLeapFromDate, _360DayFromDate)
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

vic_default_units = 'days since 0001-01-01'


def datetime_to_dmy(dt):
    dmy = ffi.new('dmy_struct *')
    dmy[0].year = dt.year
    dmy[0].month = dt.month
    dmy[0].day = dt.day
    dmy[0].dayseconds = int(dt.hour * 3600 + dt.minute * 60 + dt.second +
                            dt.microsecond / 100000)
    dmy[0].day_in_year = dt.timetuple().tm_yday
    return dmy


def dmy_to_datetime(dmy):
    day = datetime.datetime(dmy[0].year, dmy[0].month, dmy[0].day)
    td = datetime.timedelta(seconds=dmy[0].dayseconds)
    return day + td


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
    frac_day = vic_lib.fractional_day_from_dmy(dmy_feb_3_noon)
    assert frac_day == 3.5


def test_leap_year():
    for year in np.arange(1900, 2100):
        actual = calendar.isleap(year)
        assert not bool(vic_lib.leap_year(year, calendars['noleap']))
        assert actual == bool(vic_lib.leap_year(year, calendars['standard']))


@pytest.mark.skipif(nctime_unavailable, reason=nc_reason)
def test_julian_day_from_dmy(feb3_noon, dmy_feb_3_noon):
    for cal in ['standard', 'gregorian', 'proleptic_gregorian', 'julian']:
        d_vic = vic_lib.julian_day_from_dmy(dmy_feb_3_noon,
                                            calendars[cal])
        d_nc = netcdftime.JulianDayFromDate(feb3_noon, cal)
        np.testing.assert_allclose(d_vic, d_nc)


@pytest.mark.skipif(nctime_unavailable, reason=nc_reason)
def test_no_leap_day_from_dmy(feb3_noon, dmy_feb_3_noon):
    d_vic = vic_lib.no_leap_day_from_dmy(dmy_feb_3_noon)
    d_nc = netcdftime._NoLeapDayFromDate(feb3_noon)
    np.testing.assert_allclose(d_vic, d_nc)


@pytest.mark.skipif(nctime_unavailable, reason=nc_reason)
def test_all_leap_from_dmy(feb3_noon, dmy_feb_3_noon):
    d_vic = vic_lib.all_leap_from_dmy(dmy_feb_3_noon)
    d_nc = netcdftime._AllLeapFromDate(feb3_noon)
    np.testing.assert_allclose(d_vic, d_nc)


@pytest.mark.skipif(nctime_unavailable, reason=nc_reason)
def test_all_30_day_from_dmy(feb3_noon, dmy_feb_3_noon):
    d_vic = vic_lib.all_30_day_from_dmy(dmy_feb_3_noon)
    d_nc = netcdftime._360DayFromDate(feb3_noon)
    np.testing.assert_allclose(d_vic, d_nc)


@pytest.mark.skipif(nctime_unavailable, reason=nc_reason)
def test_date2num(feb3_noon, dmy_feb_3_noon):
    for cal, cal_num in calendars.items():
        for unit, unit_num in units.items():
            unit_string = '{0} since 0001-01-01'.format(unit)
            ut = utime(unit_string, calendar=cal)
            expected = ut.date2num(feb3_noon)
            actual = vic_lib.date2num(ut._jd0, dmy_feb_3_noon, 0,
                                      cal_num, unit_num)
            np.testing.assert_allclose(actual, expected)


@pytest.mark.skipif(nctime_unavailable, reason=nc_reason)
def test_dmy_julian_day(feb3_noon):
    dmy_struct = ffi.new("dmy_struct *")
    expected = feb3_noon
    for cal in ['julian', 'standard', 'gregorian', 'proleptic_gregorian']:
        expected_jd = JulianDayFromDate(expected, calendar=cal)
        vic_lib.dmy_julian_day(expected_jd, calendars[cal], dmy_struct)
        actual = dmy_to_datetime(dmy_struct)
        # assert that the difference is less than one second
        assert abs(expected - actual) < datetime.timedelta(seconds=1)


@pytest.mark.skipif(nctime_unavailable, reason=nc_reason)
def test_dmy_julian_day_timeseries():
    # Regression test for GH298, updated to handle all standard calendars
    dmy_struct = ffi.new("dmy_struct *")

    for cal in ['standard', 'gregorian', 'proleptic_gregorian']:
        dates = pd.date_range('1970-01-01', '1975-12-31', freq='6H').to_pydatetime()
        for expected in dates:
            expected_jday = JulianDayFromDate(expected)
            vic_lib.dmy_julian_day(expected_jday, calendars[cal], dmy_struct)
            actual = dmy_to_datetime(dmy_struct)
            print(cal, expected_jday, expected, actual)
            # assert that the difference is less than one second
            assert abs(expected - actual) < datetime.timedelta(seconds=1)


@pytest.mark.skipif(nctime_unavailable, reason=nc_reason)
def test_dmy_no_leap_day(feb3_noon):
    dmy_struct = ffi.new("dmy_struct *")
    expected = feb3_noon
    expected_jd = _NoLeapDayFromDate(expected)
    vic_lib.dmy_no_leap_day(expected_jd, dmy_struct)
    actual = dmy_to_datetime(dmy_struct)
    # assert that the difference is less than one second
    assert abs(expected - actual) < datetime.timedelta(seconds=1)


@pytest.mark.skipif(nctime_unavailable, reason=nc_reason)
def test_dmy_all_leap(feb3_noon):
    dmy_struct = ffi.new("dmy_struct *")
    expected = feb3_noon
    expected_jd = _AllLeapFromDate(expected)
    vic_lib.dmy_all_leap(expected_jd, dmy_struct)
    actual = dmy_to_datetime(dmy_struct)
    # assert that the difference is less than one second
    assert abs(expected - actual) < datetime.timedelta(seconds=1)


@pytest.mark.skipif(nctime_unavailable, reason=nc_reason)
def test_dmy_all_30_day(feb3_noon):
    dmy_struct = ffi.new("dmy_struct *")
    expected = feb3_noon
    expected_jd = _360DayFromDate(expected)
    vic_lib.dmy_all_30_day(expected_jd, dmy_struct)
    actual = dmy_to_datetime(dmy_struct)
    # assert that the difference is less than one second
    assert abs(expected - actual) < datetime.timedelta(seconds=1)


@pytest.mark.skipif(nctime_unavailable, reason=nc_reason)
def test_num2date():

    dates = pd.date_range(start="1900-01-01", end="2100-12-31",
                          freq='MS').to_pydatetime()

    dmy_struct = ffi.new("dmy_struct *")
    for cal, cal_num in calendars.items():
        ut = utime(vic_default_units, calendar=cal)
        for date in dates:
            num = ut.date2num(date)
            vic_lib.num2date(ut._jd0, num, 0., cal_num, units['days'],
                             dmy_struct)
            actual = dmy_to_datetime(dmy_struct)
            assert abs(date - actual) < datetime.timedelta(seconds=1)


@pytest.mark.skipif(nctime_unavailable, reason=nc_reason)
def test_make_lastday():
    # wasn't able to map lastday to a numpy array, don't know why...
    lastday = ffi.new('unsigned short int [12]', [0] * 12)
    for year in np.arange(1900, 2100):
        for cal in ['standard', 'gregorian', 'proleptic_gregorian']:
            vic_lib.make_lastday(calendars[cal], year, lastday)
            cal_dpm = [calendar.monthrange(year, month)[1] for month in
                       np.arange(1, 13)]
            np.testing.assert_equal(cal_dpm, list(lastday))
        for cal in ['noleap', '365_day']:
            vic_lib.make_lastday(calendars[cal], year, lastday)
            cal_dpm = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
            np.testing.assert_equal(cal_dpm, list(lastday))
        for cal in ['all_leap', '366_day']:
            vic_lib.make_lastday(calendars[cal], year, lastday)
            cal_dpm = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
            np.testing.assert_equal(cal_dpm, list(lastday))
        for cal in ['360_day']:
            vic_lib.make_lastday(calendars[cal], year, lastday)
            cal_dpm = [30] * 12
            np.testing.assert_equal(cal_dpm, list(lastday))


@pytest.mark.skipif(nctime_unavailable, reason=nc_reason)
def test_initialize_time():
    for cal, cal_num in calendars.items():
        for unit, unit_num in units.items():
            unit_string = '{0} since 0001-01-01'.format(unit)
            ut = utime(unit_string, calendar=cal)
            vic_lib.global_param.calendar = cal_num
            vic_lib.global_param.time_units = unit_num
            assert vic_lib.initialize_time() is None
            assert vic_lib.global_param.time_origin_num != -99999
            print(unit_string, cal)
            np.testing.assert_allclose(vic_lib.global_param.time_origin_num,
                                       ut._jd0)


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


def test_calendar_from_chars():
    for cal, expected in calendars.items():
        actual = vic_lib.calendar_from_chars(ffi.new('char[]', cal.lower().encode()))
        assert actual == expected

        actual = vic_lib.calendar_from_chars(ffi.new('char[]', cal.upper().encode()))
        assert actual == expected


def timeunits_from_chars():
    for unit, expected in units.items():
        actual = vic_lib.calendar_from_chars(ffi.new('char[]', unit.lower().encode()))
        assert actual == expected
        actual = vic_lib.calendar_from_chars(ffi.new('char[]', unit.upper().encode()))
        assert actual == expected


def test_parse_nc_time_units():
    actual_unit = ffi.new('unsigned short int *')
    actual_dmy = ffi.new('dmy_struct *')
    unit_chars = ffi.new('char[]', b'days since 0001-01-01')
    vic_lib.parse_nc_time_units(ffi.addressof(unit_chars), actual_unit, actual_dmy)
    assert actual_unit[0] == units['days']
    assert actual_dmy[0].year == 1
    assert actual_dmy[0].month == 1
    assert actual_dmy[0].day == 1
    assert actual_dmy[0].dayseconds == 0
