from vic import lib as vic_lib


def test_initialize_options():
    assert vic_lib.initialize_options() is None
    assert vic_lib.options.AboveTreelineVeg == -1


def test_initialize_options_bools():
    assert vic_lib.initialize_options() is None
    assert not vic_lib.options.BLOWING
    assert vic_lib.options.TFALLBACK
