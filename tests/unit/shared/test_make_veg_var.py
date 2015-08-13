from vic import lib as vic_lib


def test_make_veg_var():
    assert vic_lib.make_veg_var(4) is not None


def test_make_veg_var_1snowband():
    vic_lib.options.SNOW_BAND = 1
    assert vic_lib.make_veg_var(4) is not None


def test_make_veg_var_5snowband():
    vic_lib.options.SNOW_BAND = 5
    assert vic_lib.make_veg_var(4) is not None


def test_make_veg_var_5snowbands_carbon_is_true():
    vic_lib.options.SNOW_BAND = 5
    vic_lib.options.CARBON = True
    assert vic_lib.make_veg_var(4) is not None
