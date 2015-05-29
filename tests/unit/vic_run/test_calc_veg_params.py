from vic.vic import calc_veg_displacement, calc_veg_height, calc_veg_roughness


def test_calc_veg_displacement():
    assert calc_veg_displacement(0.) == 0.
    assert calc_veg_displacement(10.) < 6.8


def test_calc_veg_height():
    assert calc_veg_height(0.) == 0.
    assert calc_veg_height(10) > 10.


def test_calc_veg_roughness():
    assert calc_veg_roughness(0.) == 0.
    assert calc_veg_roughness(10.) < 10.
