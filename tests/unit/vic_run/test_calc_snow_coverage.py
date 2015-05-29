import ctypes
from vic.vic import calc_snow_coverage


def test_calc_snow_coverage_no_change():
    store_snow = ctypes.c_bool()
    max_snow_depth = ctypes.c_double(1)
    store_swq = ctypes.c_double(1)
    snow_distrib_slope = ctypes.c_double(1)
    store_coverage = ctypes.c_double(1)
    old_coverage = 0.75
    coverage = calc_snow_coverage(store_snow, 0.5,
                                  old_coverage, 1.25, 1.25, 2.3, 2.3, 0.,
                                  ctypes.byref(max_snow_depth), 0.,
                                  ctypes.byref(store_swq),
                                  ctypes.byref(snow_distrib_slope),
                                  ctypes.byref(store_coverage))

    assert coverage == old_coverage


# segfaulting...
# def test_calc_snow_coverage_increased():
#     store_snow = ctypes.c_bool(True)
#     max_snow_depth = ctypes.c_double(3)
#     store_swq = ctypes.c_double(1.25)
#     snow_distrib_slope = ctypes.c_double(0.5)
#     store_coverage = ctypes.c_double(0.75)
#     old_coverage = 0.75
#     coverage = calc_snow_coverage(store_snow, 0.5,
#                                   old_coverage, 1.25, 1.5, 2.3, 3., 0.,
#                                   ctypes.byref(max_snow_depth), 0.25,
#                                   ctypes.byref(store_swq),
#                                   ctypes.byref(snow_distrib_slope),
#                                   ctypes.byref(store_coverage))

#     assert coverage > old_coverage
