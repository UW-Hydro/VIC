import os
import subprocess
import pytest


class CESMTestException(Exception):
    pass

if os.getenv('TESTID', 'notset') != 'cesm':
    wrong_driver = True
else:
    wrong_driver = False
    test_dir = os.path.dirname(__file__)
    shared_object = os.path.join(test_dir, os.pardir, os.pardir,
                                 'vic', 'drivers', 'cesm', 'lndlib.a')
    if os.path.isfile(shared_object):
        lib_built = True
    else:
        raise CESMTestException('lndlib.a has not been built!')


def nm(shared_object, *args):
    '''wrapper of the unix nm utility'''
    if not os.path.isfile(shared_object):
        raise FileNotFoundError('%s is not a file' % shared_object)

    cmd = ['nm', '-j']
    cmd.extend(args)
    cmd.append(shared_object)
    output = subprocess.check_output(cmd)

    return output.decode().split()


@pytest.mark.skipif(wrong_driver, reason='Only run for CESM driver builds')
def test_vic_lib_present():
    assert lib_built


@pytest.mark.skipif(wrong_driver, reason='Only run for CESM driver builds')
def test_vic_run_symbols_present():
    lib_symbols = nm(shared_object)

    vic_run_symbols = ['_advect_carbon_storage',
                       '_advect_snow_storage',
                       '_advect_soil_veg_storage',
                       '_advected_sensible_heat',
                       '_alblake',
                       '_arno_evap',
                       '_calc_atmos_energy_bal',
                       '_calc_density',
                       '_calc_energy_balance_error',
                       '_calc_latent_heat_of_sublimation',
                       '_calc_latent_heat_of_vaporization',
                       '_calc_layer_average_thermal_props',
                       '_calc_outgoing_longwave',
                       '_calc_scale_height',
                       '_calc_sensible_heat',
                       '_calc_Nscale_factors',
                       '_calc_rainonly',
                       '_calc_rc',
                       '_calc_rc_ps',
                       '_calc_snow_coverage',
                       '_calc_soil_thermal_fluxes',
                       '_calc_surf_energy_bal',
                       '_calc_veg_displacement',
                       '_calc_veg_height',
                       '_calc_veg_roughness',
                       '_calc_water_balance_error',
                       '_CalcAerodynamic',
                       '_CalcBlowingSnow',
                       '_CalcIcePackEnergyBalance',
                       '_CalcSnowPackEnergyBalance',
                       '_CalcSubFlux',
                       '_canopy_assimilation',
                       '_canopy_evap',
                       '_colavg',
                       '_collect_eb_terms',
                       '_collect_wb_terms',
                       '_compute_coszen',
                       '_compute_pot_evap',
                       '_compute_runoff_and_asat',
                       '_compute_soil_resp',
                       '_compute_soil_layer_thermal_properties',
                       '_compute_zwt',
                       '_correct_precip',
                       '_darkinhib',
                       '_distribute_node_moisture_properties',
                       '_eddy',
                       '_energycalc',
                       '_error_calc_atmos_energy_bal',
                       '_error_calc_atmos_moist_bal',
                       '_error_calc_canopy_energy_bal',
                       '_error_calc_surf_energy_bal',
                       '_error_print_atmos_energy_bal',
                       '_error_print_atmos_moist_bal',
                       '_error_print_canopy_energy_bal',
                       '_error_print_solve_T_profile',
                       '_error_print_surf_energy_bal',
                       '_error_solve_T_profile',
                       '_ErrorIcePackEnergyBalance',
                       '_ErrorPrintIcePackEnergyBalance',
                       '_ErrorPrintSnowPackEnergyBalance',
                       '_ErrorSnowPackEnergyBalance',
                       '_estimate_layer_ice_content',
                       '_estimate_layer_ice_content_quick_flux',
                       '_estimate_T1',
                       '_faparl',
                       '_fda_heat_eqn',
                       '_fdjac3',
                       '_find_0_degree_fronts',
                       '_func_atmos_energy_bal',
                       '_func_atmos_moist_bal',
                       '_func_canopy_energy_bal',
                       '_func_surf_energy_bal',
                       '_funcd',
                       '_get_depth',
                       '_get_prob',
                       '_get_sarea',
                       '_get_shear',
                       '_get_thresh',
                       '_get_volume',
                       '_hiTinhib',
                       '_initialize_lake',
                       '_ice_melt',
                       '_IceEnergyBalance',
                       '_iceform',
                       '_icerad',
                       '_lakeice',
                       '_latent_heat_from_snow',
                       '_latsens',
                       '_linear_interp',
                       '_lkdrag',
                       '_MassRelease',
                       '_maximum_unfrozen_water',
                       '_new_snow_density',
                       '_newt_raph',
                       '_penman',
                       '_photosynth',
                       '_polint',
                       '_prepare_full_energy',
                       '_put_data',
                       '_qromb',
                       '_sub_with_height',
                       '_rescale_snow_energy_fluxes',
                       '_rescale_snow_storage',
                       '_rescale_soil_veg_fluxes',
                       '_rhoinit',
                       '_root_brent',
                       '_rtnewt',
                       '_runoff',
                       '_set_node_parameters',
                       '_shear_stress',
                       '_snow_albedo',
                       '_snow_density',
                       '_snow_intercept',
                       '_snow_melt',
                       '_SnowPackEnergyBalance',
                       '_soil_carbon_balance',
                       '_soil_conductivity',
                       '_soil_thermal_eqn',
                       '_solve_atmos_energy_bal',
                       '_solve_atmos_moist_bal',
                       '_solve_canopy_energy_bal',
                       '_solve_lake',
                       '_solve_snow',
                       '_solve_surf_energy_bal',
                       '_solve_T_profile',
                       '_solve_T_profile_implicit',
                       '_specheat',
                       '_StabilityCorrection',
                       '_sub_with_height',
                       '_surface_fluxes',
                       '_svp',
                       '_svp_slope',
                       '_temp_area',
                       '_tracer_mixer',
                       '_transpiration',
                       '_transport_with_height',
                       '_trapzd',
                       '_funcd',
                       '_tridia',
                       '_tridiag',
                       '_vic_run',
                       '_volumetric_heat_capacity',
                       '_water_balance',
                       '_water_energy_balance',
                       '_water_under_ice',
                       '_wrap_compute_zwt',
                       '_write_layer',
                       '_write_vegvar',
                       '_zero_output_list']

    for symbol in vic_run_symbols:
        assert symbol in lib_symbols


@pytest.mark.skipif(wrong_driver, reason='Only run for CESM driver builds')
def test_shared_symbols_present():
    lib_symbols = nm(shared_object)

    shared_symbols = ['_all_30_day_from_dmy',
                      '_all_leap_from_dmy',
                      '_calc_root_fractions',
                      '_calendar_from_char',
                      '_compute_treeline',
                      '_cmd_proc',
                      '_compress_files',
                      '_get_current_datetime',
                      '_date2num',
                      '_dmy_all_30_day',
                      '_dmy_all_leap',
                      '_dmy_julian_day',
                      '_dmy_no_leap_day',
                      '_dt_seconds_to_time_units',
                      '_display_current_settings',
                      '_fractional_day_from_dmy',
                      '_free_all_vars',
                      '_free_dmy',
                      '_free_vegcon',
                      '_get_dist',
                      '_get_parameters',
                      '_initialize_filenames',
                      '_initialize_fileps',
                      '_initialize_global',
                      '_initialize_options',
                      '_initialize_parameters',
                      '_initialize_snow',
                      '_initialize_soil',
                      '_initialize_time',
                      '_initialize_veg',
                      '_julian_day_from_dmy',
                      '_leap_year',
                      '_make_all_vars',
                      '_make_cell_data',
                      '_make_dmy',
                      '_make_energy_bal',
                      '_make_lastday',
                      '_make_snow_data',
                      '_make_veg_var',
                      '_no_leap_day_from_dmy',
                      '_num2date',
                      '_open_file',
                      '_print_cell_data',
                      '_print_dmy',
                      '_print_energy_bal',
                      '_print_filenames',
                      '_print_filep',
                      '_print_force_type',
                      '_print_global_param',
                      '_print_lake_con',
                      '_print_lake_var',
                      '_print_layer_data',
                      '_print_license',
                      '_print_option',
                      '_print_out_data',
                      '_print_out_data_file',
                      '_print_param_set',
                      '_print_parameters',
                      '_print_save_data',
                      '_print_snow_data',
                      '_print_soil_con',
                      '_print_veg_con',
                      '_print_veg_lib',
                      '_print_veg_var',
                      '_print_version',
                      '_print_usage',
                      '_soil_moisture_from_water_table',
                      '_valid_date',
                      '_validate_parameters']

    for symbol in shared_symbols:
        assert symbol in lib_symbols


@pytest.mark.skipif(wrong_driver, reason='Only run for CESM driver builds')
def test_cesm_symbols_present():
    lib_symbols = nm(shared_object)

    cesm_symbols = ['_add_nveg_to_global_domain',
                    '_alloc_atmos',
                    '_alloc_veg_hist',
                    '_air_density',
                    '_average',
                    '_create_output_list',
                    '_free_atmos',
                    '_free_out_data',
                    '_free_veg_hist',
                    '_get_global_domain',
                    '_get_global_param',
                    '_get_nc_dimension',
                    '_get_nc_field_double',
                    '_get_nc_field_float',
                    '_get_nc_field_int',
                    '_init_output_list',
                    '_initialize_domain',
                    '_initialize_energy',
                    '_initialize_history_file',
                    '_initialize_location',
                    '_initialize_l2x_data',
                    '_initialize_model_state',
                    '_initialize_soil_con',
                    '_initialize_state_file',
                    '_initialize_veg_con',
                    '_initialize_x2l_data',
                    '_open_file',
                    '_parse_output_info',
                    '_print_atmos_data',
                    '_print_domain',
                    '_print_l2x_data',
                    '_print_location',
                    '_print_nc_file',
                    '_print_nc_var',
                    '_print_vic_clock',
                    '_print_veg_con_map',
                    '_print_x2l_data',
                    '_put_nc_field_double',
                    '_put_nc_field_int',
                    '_q_to_vp',
                    '_sprint_location',
                    '_vic_alloc',
                    '_vic_nc_info',
                    '_vic_finalize',
                    '_vic_force',
                    '_vic_cesm_init',
                    '_vic_cesm_init_mpi',
                    '_vic_cesm_run',
                    '_vic_cesm_run_model',
                    '_vic_cesm_final',
                    '_vic_init',
                    '_vic_init_output',
                    '_vic_restore',
                    '_vic_start',
                    '_vic_store',
                    '_vic_write',
                    '_will_it_snow']

    for symbol in cesm_symbols:
        assert symbol in lib_symbols


@pytest.mark.skipif(wrong_driver, reason='Only run for CESM driver builds')
def test_vic_mpi_symbols_present():
    lib_symbols = nm(shared_object)

    vic_mpi_symbols = ['_create_MPI_global_struct_type',
                       '_create_MPI_location_struct_type',
                       '_create_MPI_nc_file_struct_type',
                       '_create_MPI_option_struct_type',
                       '_create_MPI_param_struct_type',
                       '_gather_put_nc_field_double',
                       '_gather_put_nc_field_int',
                       '_get_scatter_nc_field_double',
                       '_get_scatter_nc_field_float',
                       '_get_scatter_nc_field_int',
                       '_initialize_mpi',
                       '_map',
                       '_mpi_map_decomp_domain']

    for symbol in vic_mpi_symbols:
        assert symbol in lib_symbols


@pytest.mark.skipif(wrong_driver, reason='Only run for CESM driver builds')
def test_fortran_symbols_present():
    lib_symbols = nm(shared_object)

    vic_fortran_symbols = ['_vic_cesm_init',
                           '_vic_cesm_run',
                           '_vic_cesm_final']

    for symbol in vic_fortran_symbols:
        assert symbol in lib_symbols
