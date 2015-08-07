from vic import lib as vic_lib


def test_svp():
    assert vic_lib.svp(0.) > 0.
    assert vic_lib.svp(0.) > vic_lib.svp(-1.)


def test_svp_slope():
    assert vic_lib.svp_slope(0.) > 0.
    assert vic_lib.svp_slope(0.) > vic_lib.svp_slope(-1.)
