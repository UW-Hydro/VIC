from vic import VIC_DRIVER
from vic import lib as vic_lib
from vic.pycompat import pylong


def test_print_version():
    assert VIC_DRIVER == b'Python'
    assert vic_lib.print_version(VIC_DRIVER) is None


def test_set_global_scalars():
    vic_lib.flag = 1
    vic_lib.NR = 6
    vic_lib.NF = 3

    assert isinstance(vic_lib.flag, int)
    assert isinstance(vic_lib.NR, pylong)
    assert isinstance(vic_lib.NF, pylong)
    assert vic_lib.flag == 1
    assert vic_lib.NR == 6
    assert vic_lib.NF == 3


def test_globals_are_initialized():
    assert vic_lib.global_param.atmos_dt == -99999.
    assert vic_lib.options.AboveTreelineVeg == -1
    assert vic_lib.param.LAPSE_RATE == -0.0065
