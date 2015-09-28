from vic import lib as vic_lib


def test_initialize_parameters():
    assert vic_lib.initialize_parameters() is None
    assert vic_lib.param.LAPSE_RATE == -0.0065
