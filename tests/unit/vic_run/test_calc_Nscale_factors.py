import numpy as np
from vic.vic import ffi
from vic import lib as vic_lib


def test_calc_nscale_factors():
    canopy_layer_bnd = np.ones(3, dtype=np.float64)
    nscale_factor = np.empty(3, dtype=np.float64)
    flag = ffi.cast('char', True)

    assert vic_lib.calc_Nscale_factors(
        flag, ffi.cast('double *', canopy_layer_bnd.ctypes.data), 3., 
        0.5, ffi.cast('double *', nscale_factor.ctypes.data)) is None
