from vic import lib as vic_lib


def test_get_dist():
    d0 = vic_lib.get_dist(0., 0., 0, 1.)
    assert d0 > 111000
    assert d0 < 112000

    d1 = vic_lib.get_dist(0., 0., 1., 1.)
    assert d1 > d0

    d2 = vic_lib.get_dist(75., 0., 75., 1.)
    assert d2 < d0

    assert vic_lib.get_dist(0., 0., 1., 0) == vic_lib.get_dist(74., 0., 75., 0.)

    assert vic_lib.get_dist(1, 0., -1, 0) > 0.
    assert vic_lib.get_dist(0, 0, 0, -1) > 0.
