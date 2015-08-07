import numpy as np
from vic.vic import ffi
from vic import lib as vic_lib


def test_calc_aerodynamic():

    r_a = np.zeros(3, dtype=np.float64)
    u = np.zeros(3, dtype=np.float64)
    displacement = np.zeros(3, dtype=np.float64)
    ref_height = np.zeros(3, dtype=np.float64)
    roughness = np.zeros(3, dtype=np.float64)

    assert vic_lib.CalcAerodynamic(
        ffi.cast('_Bool', True), 40., 0.2, 0.0005, 0.001, 0.5,
        ffi.cast('double *', r_a.ctypes.data),
        ffi.cast('double *', u.ctypes.data),
        ffi.cast('double *', displacement.ctypes.data),
        ffi.cast('double *', ref_height.ctypes.data),
        ffi.cast('double *', roughness.ctypes.data)) == 0.
