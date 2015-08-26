from vic import lib as vic_lib
import numpy as np


# TODO: Include physical constants in library
JOULES_PER_CAL = 4.1868
GRAMS_PER_KG = 1000
CONST_LATVAP = 2.501e6


def test_calc_latent_heat_of_sublimation():
    ls0 = vic_lib.calc_latent_heat_of_sublimation(0)
    assert ls0 == (677. * JOULES_PER_CAL * GRAMS_PER_KG)
    ls20 = vic_lib.calc_latent_heat_of_sublimation(20)
    assert ls20 < ls0


def test_calc_latent_heat_of_vaporization():
    lv0 = vic_lib.calc_latent_heat_of_vaporization(0)
    assert lv0 == CONST_LATVAP
    lv20 = vic_lib.calc_latent_heat_of_vaporization(20)
    assert lv20 < lv0


def test_calc_outgoing_longwave():
    lwout273_1e = vic_lib.calc_outgoing_longwave(273., 1.)
    assert lwout273_1e > 0.
    lwout273_97e = vic_lib.calc_outgoing_longwave(273., 0.97)
    assert lwout273_97e < lwout273_1e
    np.testing.assert_almost_equal(lwout273_97e, lwout273_1e * 0.97)


def test_calc_scale_height():
    h00 = vic_lib.calc_scale_height(0., 0.)
    assert h00 > 0.
    h90 = vic_lib.calc_scale_height(9., 0.)
    assert h90 > h00
    h09 = vic_lib.calc_scale_height(0., 9.)
    assert h09 < h00


def test_calc_sensible_heat():
    assert vic_lib.calc_sensible_heat(2.75, 0., 1., 100) < 0.
    assert vic_lib.calc_sensible_heat(2.75, 0., 0., 100) == 0.
    assert vic_lib.calc_sensible_heat(2.75, 1., 0., 100) > 0.
