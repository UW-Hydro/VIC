# import ctypes
from vic.vic import make_dmy, global_param


def test_make_dmy():
    global_param.startsec = 0
    global_param.startday = 1
    global_param.startmonth = 1
    global_param.startyear = 2015
    global_param.nrecs = 30

    # dmy = make_dmy(global_param)
