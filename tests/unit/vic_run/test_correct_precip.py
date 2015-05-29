import ctypes
from vic.vic import correct_precip, RAIN, SNOW


def test_correct_precip():
    guage_correction = (2 * ctypes.c_double)(0)
    assert correct_precip(ctypes.byref(guage_correction),
                          4., 2., 0.4, 0.005) is None
    assert guage_correction[SNOW] > 1.
    assert guage_correction[RAIN] > 1.
    assert guage_correction[SNOW] > guage_correction[RAIN]
