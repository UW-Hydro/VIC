from vic.vic import make_energy_bal, options


def test_make_energy_bal():
    assert make_energy_bal(4) is not None


def test_make_energy_bal_1snowband():
    options.SNOW_BAND = 1
    assert make_energy_bal(4) is not None


def test_make_energy_bal_5snowband():
    options.SNOW_BAND = 5
    assert make_energy_bal(4) is not None
