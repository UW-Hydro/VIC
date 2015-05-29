from vic.vic import initialize_options, options


def test_initialize_options():
    assert initialize_options() is None
    assert options.AboveTreelineVeg == -1


def test_initialize_options_bools():
    assert initialize_options() is None
    assert not options.BLOWING
    assert not options.PRT_SNOW_BAND
