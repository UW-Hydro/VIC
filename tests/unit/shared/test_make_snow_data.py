from vic.vic import make_snow_data, options


def test_make_snow_data():
    assert make_snow_data(4) is not None


def test_make_snow_data_1snowband():
    options.SNOW_BAND = 1
    assert make_snow_data(4) is not None


def test_make_snow_data_5snowband():
    options.SNOW_BAND = 5
    assert make_snow_data(4) is not None
