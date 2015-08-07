from vic.vic import ffi
from vic import lib as vic_lib


def test_snow_albedo_new_snow():
    assert vic_lib.snow_albedo(
        0.1, 1., 0.7, -77, 3600., 0,
        ffi.cast('char', False)) == vic_lib.param.SNOW_NEW_SNOW_ALB


def test_snow_albedo_old_snow():
    orig_albedo = 0.7
    assert vic_lib.snow_albedo(0., 1., orig_albedo, -77, 3600., 23,
                               ffi.cast('char', False)) < orig_albedo
    assert vic_lib.snow_albedo(0., 1., orig_albedo, -77, 3600., 23,
                               ffi.cast('char', True)) < orig_albedo
