from vic.vic import svp, svp_slope


def test_svp():
    assert svp(0.) > 0.
    assert svp(0.) > svp(-1.)


def test_svp_slope():
    assert svp_slope(0.) > 0.
    assert svp_slope(0.) > svp_slope(-1.)
