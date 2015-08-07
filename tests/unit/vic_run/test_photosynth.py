from vic import lib as vic_lib


def test_hitinhib():
    temps = [-10, 0., 10., 20., 30., 55.]
    for t in temps:
        assert vic_lib.hiTinhib(t) >= 0.5


def test_darkinhib():
    irrs = [-10, 0., 10., 20., 30., 55.]
    for i in irrs:
        if i < 0.:
            assert vic_lib.darkinhib(i) == 0.
        else:
            assert vic_lib.darkinhib(i) > 0.
