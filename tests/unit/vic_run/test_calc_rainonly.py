from vic.vic import calc_rainonly


def test_rainonly():
    assert calc_rainonly(0., 4., 1., -1.) == 2.  # 50/50 split
    assert calc_rainonly(3., 4., 1., -1.) == 4.  # all rain
    assert calc_rainonly(-2., 4., 1., -1.) == 0.  # no rain
