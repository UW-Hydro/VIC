from vic import lib as vic_lib


def test_stability_correction():
    assert vic_lib.StabilityCorrection(10, 2., 0., 2., 2., 0.5) > 0.
