from vic.vic import initialize_global, global_param, MISSING, SEC_PER_DAY


def test_initialize_global():
    assert initialize_global() is None
    assert global_param.dt == MISSING
    assert global_param.out_dt == SEC_PER_DAY
