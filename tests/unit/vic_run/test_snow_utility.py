from vic.vic import snow_albedo, param


def test_snow_albedo_new_snow():
    assert snow_albedo(0.1, 1., 0.7, -77, 3600., 0, False) == \
        param.SNOW_NEW_SNOW_ALB


def test_snow_albedo_old_snow():
    orig_albedo = 0.7
    assert snow_albedo(0., 1., orig_albedo, -77, 3600., 23, False) < \
        orig_albedo
    assert snow_albedo(0., 1., orig_albedo, -77, 3600., 23, True) < \
        orig_albedo
