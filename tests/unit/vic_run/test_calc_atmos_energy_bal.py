from vic.vic import ffi
from vic import lib as vic_lib


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
    error = ffi.new('double *')
    error[0] = 0  # (* double)
    latent_heat = ffi.new('double *')
    latent_heat[0] = 0  # (* double)
    latent_heat_sub = ffi.new('double *')
    latent_heat_sub[0] = 0  # (* double)
    net_long_atmos = ffi.new('double *')
    net_long_atmos[0] = 0  # (* double)
    net_short_atmos = ffi.new('double *')
    net_short_atmos[0] = 0  # (* double)
    sensible_heat = ffi.new('double *')
    sensible_heat[0] = 0  # (* double)
    tcanopy_fbflag = ffi.new('_Bool *')
    tcanopy_fbflag[0] = 0  # (* bool)
    tcanopy_fbcount = ffi.new('unsigned *')
    tcanopy_fbcount[0] = 0  # (* unsigned)
    vic_lib.calc_atmos_energy_bal(in_over_sensible,
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
                                  error,
                                  latent_heat,
                                  latent_heat_sub,
                                  net_long_atmos,
                                  net_short_atmos,
                                  sensible_heat,
                                  tcanopy_fbflag,
                                  tcanopy_fbcount)

    assert error[0] != 0.
    assert latent_heat[0] != 0.
    assert latent_heat_sub[0] != 0.
    assert net_long_atmos[0] != 0.
    assert net_short_atmos[0] == 0.
    assert sensible_heat[0] == 0.
    assert not tcanopy_fbflag[0]
    assert tcanopy_fbcount[0] == 0.
