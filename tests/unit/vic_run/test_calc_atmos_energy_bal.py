import ctypes
from vic.vic import calc_atmos_energy_bal


def test_calc_atmos_energy_bal():
    # note that this test uses default options
    in_over_sensible = 15.  # (double)
    in_under_sensible = 5.  # (double)
    latent_heat_over = 15.  # (double)
    latent_heat_under = 5.  # (double)
    latent_heat_sub_over = 5.  # (double)
    latent_heat_sub_under = 1.  # (double)
    net_long_over = 50.  # (double)
    net_long_under = 30.  # (double)
    net_short_over = 0.  # (double)
    net_short_under = 0.  # (double)
    r_a = 500.  # (double)
    t_air = 2.  # (double)
    atmos_density = 1.225  # (double)
    error = ctypes.c_double()  # (* double)
    latent_heat = ctypes.c_double()  # (* double)
    latent_heat_sub = ctypes.c_double()  # (* double)
    net_long_atmos = ctypes.c_double()  # (* double)
    net_short_atmos = ctypes.c_double()  # (* double)
    sensible_heat = ctypes.c_double()  # (* double)
    tcanopy_fbflag = ctypes.c_bool()  # (* bool)
    tcanopy_fbcount = ctypes.c_uint()  # (* unsigned)
    calc_atmos_energy_bal(in_over_sensible,
                          in_under_sensible,
                          latent_heat_over,
                          latent_heat_under,
                          latent_heat_sub_over,
                          latent_heat_sub_under,
                          net_long_over,
                          net_long_under,
                          net_short_over,
                          net_short_under,
                          r_a,
                          t_air,
                          atmos_density,
                          ctypes.byref(error),
                          ctypes.byref(latent_heat),
                          ctypes.byref(latent_heat_sub),
                          ctypes.byref(net_long_atmos),
                          ctypes.byref(net_short_atmos),
                          ctypes.byref(sensible_heat),
                          ctypes.byref(tcanopy_fbflag),
                          ctypes.byref(tcanopy_fbcount))

    assert error.value != 0.
    assert latent_heat.value != 0.
    assert latent_heat_sub.value != 0.
    assert net_long_atmos.value != 0.
    assert net_short_atmos.value == 0.
    assert sensible_heat.value == 0.
    assert not tcanopy_fbflag.value
    assert tcanopy_fbcount.value == 0.
