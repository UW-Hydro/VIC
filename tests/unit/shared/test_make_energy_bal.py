from vic import lib as vic_lib


def test_make_energy_bal():
    assert vic_lib.make_energy_bal(4) is not None


def test_make_energy_bal_1snowband():
    vic_lib.options.SNOW_BAND = 1
    assert vic_lib.make_energy_bal(4) is not None


def test_make_energy_bal_5snowband():
    vic_lib.options.SNOW_BAND = 5
    assert vic_lib.make_energy_bal(4) is not None
