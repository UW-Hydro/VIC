from vic import lib as vic_lib


def test_rainonly():
    assert vic_lib.calc_rainonly(0., 4., 1., -1.) == 2.  # 50/50 split
    assert vic_lib.calc_rainonly(3., 4., 1., -1.) == 4.  # all rain
    assert vic_lib.calc_rainonly(-2., 4., 1., -1.) == 0.  # no rain
