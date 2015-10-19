from vic.vic import ffi
from vic import lib as vic_lib


def test_calc_snow_coverage_no_change():
    store_snow = ffi.new('_Bool *')
    max_snow_depth = ffi.new('double *')
    store_swq = ffi.new('double *')
    snow_distrib_slope = ffi.new('double *')
    store_coverage = ffi.new('double *')
    old_coverage = 0.75
    coverage = vic_lib.calc_snow_coverage(
        store_snow, 0.5, old_coverage, 1.25, 1.25, 2.3, 2.3, 0.,
        max_snow_depth, 0., store_swq, snow_distrib_slope, store_coverage)
    assert coverage == old_coverage


def test_calc_snow_coverage_increased():
    store_snow = ffi.new('_Bool *')
    store_snow[0] = True
    max_snow_depth = ffi.new('double *')
    max_snow_depth[0] = 3.
    store_swq = ffi.new('double *')
    store_swq[0] = 0.5
    snow_distrib_slope = ffi.new('double *')
    snow_distrib_slope[0] = 0.5
    store_coverage = ffi.new('double *')
    store_coverage[0] = 0.75
    old_coverage = 0.75
    coverage = vic_lib.calc_snow_coverage(
        store_snow, 0.5, old_coverage, 1.25, 1.5, 2.3, 3., 0., max_snow_depth,
        0.25, store_swq, snow_distrib_slope, store_coverage)

    assert coverage > old_coverage
