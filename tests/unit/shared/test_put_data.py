import numpy as np
from vic import lib as vic_lib

np.random.seed(1234)


def test_calc_water_balance_error():
    steps = 100
    recs = np.arange(0, steps, dtype=np.int)
    inflows = np.random.random(steps)
    outflows = np.random.random(steps)
    storages = np.random.random(steps) - 0.5
    last_storages = np.roll(storages, 1)

    assert vic_lib.calc_water_balance_error(0., 0., 0., 0.) == 0.
    for rec in recs:
        assert vic_lib.calc_water_balance_error(
            inflows[rec], outflows[rec], storages[rec], last_storages[rec]) != 0.


def test_calc_energy_balance_error():
    steps = 100
    recs = np.arange(0, steps, dtype=np.int)
    net_rads = np.random.random(steps)
    latents = np.random.random(steps) - 0.25
    sensibles = np.random.random(steps) - 0.2
    grnd_fluxes = np.random.random(steps) - 0.3
    snow_fluxes = np.random.random(steps) - 0.15
    assert vic_lib.calc_energy_balance_error(0., 0., 0., 0., 0.) == 0.
    for rec in recs:
        assert vic_lib.calc_energy_balance_error(
            net_rads[rec], latents[rec], sensibles[rec], grnd_fluxes[rec],
            snow_fluxes[rec]) != 0.
