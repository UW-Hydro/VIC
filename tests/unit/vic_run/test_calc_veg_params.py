from vic import lib as vic_lib


def test_calc_veg_displacement():
    assert vic_lib.calc_veg_displacement(0.) == 0.
    assert vic_lib.calc_veg_displacement(10.) < 6.8


def test_calc_veg_height():
    assert vic_lib.calc_veg_height(0.) == 0.
    assert vic_lib.calc_veg_height(10) > 10.


def test_calc_veg_roughness():
    assert vic_lib.calc_veg_roughness(0.) == 0.
    assert vic_lib.calc_veg_roughness(10.) < 10.
