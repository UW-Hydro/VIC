from vic.vic import initialize_parameters, param


def test_initialize_parameters():
    assert initialize_parameters() is None
    assert param.LAPSE_RATE == -0.0065
