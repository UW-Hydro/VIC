from vic.vic import StabilityCorrection


def test_stability_correction():
    assert StabilityCorrection(10, 2., 0., 2., 2., 0.5) > 0.
