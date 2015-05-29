from vic.vic import make_veg_var, options


def test_make_veg_var():
    assert make_veg_var(4) is not None


def test_make_veg_var_1snowband():
    options.SNOW_BAND = 1
    assert make_veg_var(4) is not None


def test_make_veg_var_5snowband():
    options.SNOW_BAND = 5
    assert make_veg_var(4) is not None


def test_make_veg_var_5snowbands_carbon_is_true():
    options.SNOW_BAND = 5
    options.carbon = True
    assert make_veg_var(4) is not None
