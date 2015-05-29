import ctypes
from vic.vic import calc_Nscale_factors


def test_calc_nscale_factors():
    canopy_layer_bnd = ctypes.c_double()
    nscale_factor = ctypes.c_double()

    assert calc_Nscale_factors(True,
                               ctypes.byref(canopy_layer_bnd),
                               3., 45, 0., 0., 45,
                               ctypes.byref(nscale_factor)) is None
