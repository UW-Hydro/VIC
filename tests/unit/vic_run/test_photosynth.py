from vic.vic import hiTinhib, darkinhib


def test_hitinhib():
    temps = [-10, 0., 10., 20., 30., 55.]
    for t in temps:
        assert hiTinhib(t) >= 0.5


def test_darkinhib():
    irrs = [-10, 0., 10., 20., 30., 55.]
    for i in irrs:
        if i < 0.:
            assert darkinhib(i) == 0.
        else:
            assert darkinhib(i) > 0.
