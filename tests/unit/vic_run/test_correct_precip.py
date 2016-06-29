import numpy as np
from vic.vic import ffi
from vic import lib as vic_lib

# TODO: include compile time constants to lib (requires enum in vic_def.h)
RAIN = 0
SNOW = 1


def test_correct_precip():
    guage_correction = np.zeros(2, dtype=np.float64)

    assert vic_lib.correct_precip(
        ffi.cast('double *', guage_correction.ctypes.data),
        4., 2., 0.4, 0.005) is None
    assert guage_correction[SNOW] > 1.
    assert guage_correction[RAIN] > 1.
    assert guage_correction[SNOW] > guage_correction[RAIN]
