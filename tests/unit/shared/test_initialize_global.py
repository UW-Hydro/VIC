from vic import lib as vic_lib


def test_initialize_global():
    assert vic_lib.initialize_global() is None
    assert vic_lib.global_param.dt == -99999
