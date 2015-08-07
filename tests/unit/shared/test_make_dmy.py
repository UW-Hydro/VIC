# import ctypes
from vic import lib as vic_lib
from vic import ffi


def test_make_dmy():
    vic_lib.global_param.model_steps_per_day = 1
    vic_lib.global_param.dt = 86400
    vic_lib.global_param.startsec = 0
    vic_lib.global_param.startday = 1
    vic_lib.global_param.startmonth = 1
    vic_lib.global_param.startyear = 2015
    vic_lib.global_param.nrecs = 30

    for calendar in ['CALENDAR_STANDARD', 'CALENDAR_GREGORIAN',
                     'CALENDAR_PROLEPTIC_GREGORIAN', 'CALENDAR_NOLEAP',
                     'CALENDAR_365_DAY', 'CALENDAR_360_DAY', 'CALENDAR_JULIAN',
                     'CALENDAR_ALL_LEAP', 'CALENDAR_366_DAY']:
        vic_lib.global_param.calendar = getattr(vic_lib, calendar)
        dmy = vic_lib.make_dmy(ffi.addressof(vic_lib.global_param))

        assert dmy is not None
        for i in range(vic_lib.global_param.nrecs):
            assert dmy[i].year == 2015
            assert dmy[i].month == 1
            assert dmy[i].day == i + 1
            assert dmy[i].dayseconds == 0
