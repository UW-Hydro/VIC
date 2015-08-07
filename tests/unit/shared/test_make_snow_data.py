from vic import lib as vic_lib


def test_make_snow_data():
    assert vic_lib.make_snow_data(4) is not None


def test_make_snow_data_1snowband():
    vic_lib.options.SNOW_BAND = 1
    assert vic_lib.make_snow_data(4) is not None


def test_make_snow_data_5snowband():
    vic_lib.options.SNOW_BAND = 5
    assert vic_lib.make_snow_data(4) is not None
