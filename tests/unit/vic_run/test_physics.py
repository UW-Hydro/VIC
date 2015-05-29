from vic.vic import (calc_latent_heat_of_sublimation,
                     calc_latent_heat_of_vaporization,
                     calc_outgoing_longwave,
                     calc_scale_height,
                     calc_sensible_heat,
                     JOULES_PER_CAL, GRAMS_PER_KG, CONST_LATVAP)
import numpy as np


def test_calc_latent_heat_of_sublimation():
    ls0 = calc_latent_heat_of_sublimation(0)
    assert ls0 == (677. * JOULES_PER_CAL * GRAMS_PER_KG)
    ls20 = calc_latent_heat_of_sublimation(20)
    assert ls20 < ls0


def test_calc_latent_heat_of_vaporization():
    lv0 = calc_latent_heat_of_vaporization(0)
    assert lv0 == CONST_LATVAP
    lv20 = calc_latent_heat_of_vaporization(20)
    assert lv20 < lv0


def test_calc_outgoing_longwave():
    lwout273_1e = calc_outgoing_longwave(273., 1.)
    assert lwout273_1e > 0.
    lwout273_97e = calc_outgoing_longwave(273., 0.97)
    assert lwout273_97e < lwout273_1e
    np.testing.assert_almost_equal(lwout273_97e, lwout273_1e * 0.97)


def test_calc_scale_height():
    h00 = calc_scale_height(0., 0.)
    assert h00 > 0.
    h90 = calc_scale_height(9., 0.)
    assert h90 > h00
    h09 = calc_scale_height(0., 9.)
    assert h09 < h00


def test_calc_sensible_heat():
    assert calc_sensible_heat(2.75, 0., 1., 100) < 0.
    assert calc_sensible_heat(2.75, 0., 0., 100) == 0.
    assert calc_sensible_heat(2.75, 1., 0., 100) > 0.
