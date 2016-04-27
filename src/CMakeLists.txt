cmake_minimum_required(VERSION 2.4)
project(vicNl)

set(CMAKE_C_COMPILER gcc)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}  -I. -g -Wall -Wno-unused")
include_directories( "." )

set(SOURCE_FILES
    advected_sensible_heat.c
    alloc_atmos.c
    alloc_veg_hist.c
    arno_evap.c
    calc_air_temperature.c
    calc_atmos_energy_bal.c
    calc_longwave.c
    calc_Nscale_factors.c
    calc_rainonly.c
    calc_root_fraction.c
    calc_snow_coverage.c
    calc_surf_energy_bal.c
    calc_veg_params.c
    calc_water_energy_balance_errors.c
    CalcAerodynamic.c
    CalcBlowingSnow.c
    canopy_assimilation.c
    canopy_evap.c
    check_files.c
    check_state_file.c
    close_files.c
    cmd_proc.c
    compress_files.c
    compute_coszen.c
    compute_pot_evap.c
    compute_soil_resp.c
    compute_treeline.c
    compute_zwt.c
    correct_precip.c
    display_current_settings.c
    estimate_T1.c
    faparl.c
    free_all_vars.c
    free_vegcon.c
    frozen_soil.c
    full_energy.c
    func_atmos_energy_bal.c
    func_atmos_moist_bal.c
    func_canopy_energy_bal.c
    func_surf_energy_bal.c
    get_dist.c
    get_force_type.c
    get_global_param.c
    global.h
    ice_melt.c
    IceEnergyBalance.c
    initialize_atmos.c
    initialize_global.c
    initialize_lake.c
    initialize_model_state.c
    initialize_snow.c
    initialize_soil.c
    initialize_veg.c
    LAKE.h
    lakes.eb.c
    latent_heat_from_snow.c
    make_all_vars.c
    make_cell_data.c
    make_dmy.c
    make_energy_bal.c
    make_in_and_outfiles.c
    make_snow_data.c
    make_veg_var.c
    massrelease.c
    modify_Ksat.c
    mtclim_constants_vic.h
    mtclim_parameters_vic.h
    mtclim_vic.c
    mtclim_wrapper.c
    newt_raph_func_fast.c
    nrerror.c
    open_file.c
    open_state_file.c
    output_list_utils.c
    parse_output_info.c
    penman.c
    photosynth.c
    prepare_full_energy.c
    print_library.c
    put_data.c
    read_atmos_data.c
    read_forcing_data.c
    read_initial_model_state.c
    read_lakeparam.c
    read_snowband.c
    read_soilparam.c
    read_veglib.c
    read_vegparam.c
    root_brent.c
    runoff.c
    set_output_defaults.c
    snow.h
    snow_intercept.c
    snow_melt.c
    snow_utility.c
    SnowPackEnergyBalance.c
    soil_carbon_balance.c
    soil_conduction.c
    soil_thermal_eqn.c
    solve_snow.c
    StabilityCorrection.c
    surface_fluxes.c
    svp.c
    vicerror.c
    vicNl.h
    vicNl_def.h
    water_energy_balance.c
    water_under_ice.c
    write_data.c
    write_forcing_file.c
    write_header.c
    write_layer.c
    write_model_state.c
    write_vegvar.c estimate_sw_by_sd.c)
set(MAIN_FILE vicNl.c)

add_library(vic SHARED ${SOURCE_FILES})
set_target_properties(vic PROPERTIES VERSION 4.2 SOVERSION 0.1)

add_executable(vicNl ${MAIN_FILE})
TARGET_LINK_LIBRARIES(vic m)
TARGET_LINK_LIBRARIES(vicNl m)
TARGET_LINK_LIBRARIES(vicNl vic)
