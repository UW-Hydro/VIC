from vic.vic import advected_sensible_heat
import numpy as np


def test_advected_sensible_heat():
    # negative when soil temperature is less than air temperature
    assert advected_sensible_heat(0.5, 1.225, -3, 0.5, 1000) < 0
    assert advected_sensible_heat(0.5, 1.225, -3, -1., 1000) < 0
    # positive when soil temperature is g.t. air temperature
    assert advected_sensible_heat(0.5, 1.225, 3, 0., 1000) > 0
    # zero when soil and air temp is the same
    assert advected_sensible_heat(0.5, 1.225, 3, 3., 1000) == 0.
    # zero when snow coverage is 1
    assert advected_sensible_heat(1., 1.225, 3, 0., 1000) == 0.
    # .inf when snow coverage is 0
    assert advected_sensible_heat(0, 1.225, 3, 0., 1000) == np.inf
