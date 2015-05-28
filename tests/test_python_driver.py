import os
import ctypes


def test_vic_shared_object_present():
    from vic.vic import _libs
    from vic.vic import vic_core

    assert vic_core in _libs
    assert os.path.isfile(_libs[vic_core]._name)


def test_set_global_scalars():
    from vic.vic import flag, NR, NF

    flag.value = 1
    NR.value = 6
    NF.value = 3

    assert type(flag) == ctypes.c_int
    assert type(NR) == ctypes.c_ulong
    assert type(NF) == ctypes.c_ulong

    assert flag.value == 1
    assert NR.value == 6
    assert NF.value == 3


def test_print_version():
    from vic.vic import print_version, VIC_DRIVER

    assert VIC_DRIVER == "Python"
    assert print_version(VIC_DRIVER) is None


def test_globals_are_initialized():
    from vic.vic import MISSING, global_param, options, param

    assert global_param.atmos_dt == MISSING
    assert options.AboveTreelineVeg == -1
    assert param.LAPSE_RATE == -0.0065
