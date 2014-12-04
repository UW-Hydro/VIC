# VIC Makefile
# Modifications:
# 27-May-2003 Replaced read_vegparam by read_vegparam_LAI			KAC
# 12-Nov-2003 Added "make depend" to the "all" and "default" options.
#             This way, if a user always types "make", the user is
#             guaranteed to have a .depend file and therefore *.o will
#             always be recompiled whenever a .h file is updated.  The
#             user can override this behavior by typing "make model",
#             which doesn't invoke "make depend".
# 24-Mar-2005 Added 2 new files: conv_force_vic2alma.c and
#	      conv_results_vic2alma.c.						TJB
# 17-Apr-2005 Added vicInterp target.						TJB
# 2006-Sep-23 (Port from 4.0.6) Implemented flexible output configuration.	TJB
#	      Removed 2 files:
#		conv_force_vic2alma.c
#		conv_results_vic2alma.c
#	      Added the following files:
#		calc_water_energy_balance_errors.c
#		output_list_utils.c
#		parse_output_info.c
#		set_output_defaults.c
# 2006-Nov-30 Changed ".c" to ".o" for:
#               output_list_utils.o
#               parse_output_info.o
#               set_output_defaults.o
# 2007-Jan-15 Added PRT_HEADER option; added write_header.c.			TJB
# 2007-Apr-24 Added newt_raph_func_fast.c for IMPLICIT option.			JCA
# 2007-Nov-06 Added get_dist.c.							TJB
# 2008-Feb-14 Removed -g from normal compiling option.  Changed "vicInterp"
#	      to "vicDisagg".							TJB
# 2009-Jun-09 Added compute_pot_evap.c.						TJB
# 2009-Jul-31 Removed wetland_energy.c.						TJB
# 2010-Dec-01 Added compute_zwt.c.						TJB
# 2011-Nov-04 Renamed mtclim* files to remove version number from filenames.	TJB
# 2012-Jan-16 Removed files associated with LINK_DEBUG code:
#             open_debug.c 
#             store_moisture_for_debug.c
#             write_atmosdata.c
#             write_debug.c
# 2013-Jul-25 Added compute_coszen.c.						TJB
# 2013-Jul-25 Added calc_Nscale_factors.c, canopy_assimilation.c, faparl.c,
# 	      photosynth.c.							TJB
# 2013-Jul-25 Added compute_soil_resp.c and soil_carbon_balance.c.		TJB
# 2013-Dec-26 Removed calc_forcing_stats.c.					TJB
# 2014-Mar-25 Removed calc_cloud_cover_fraction.c       BN
# 2014-Mar-26 Removed compute_dz.c
# 2014-Mar-26 Removed files with unused functions:
#             write_snowparam.c
#             write_soilparam.c
#             write_vegparam.c
# 2014-Sep-03 Added CN (ability to compile CN F90 code).              MAB
#
# -----------------------------------------------------------------------

# -----------------------------------------------------------------------
# SET ENVIRONMENT-SPECIFIC OPTIONS HERE
# -----------------------------------------------------------------------

# VIC DRIVER PATH (include and src are subdirs of this)
DRIVERPATH = /home/brunke/VICCN5/drivers/classic

# VIC RUN PATH
VICPATH = /home/brunke/VICCN5/vic_run

# CN PATH
CNPATH = /home/brunke/VICCN5/cn

# Set SHELL = your shell here
SHELL = /bin/bash

# Set CC = your compiler here
CC = gcc

# Set F90 = your Fortran compiler here
F90 = ifort

# set includes
INCLUDES = -I ${DRIVERPATH}/include -I ${VICPATH}/include

# Uncomment for normal optimized code flags (fastest run option)
#CFLAGS  = -O3 -Wall -Wno-unused
LIBRARY = -lm

# Uncomment to include debugging information
CFLAGS  =  ${INCLUDES} -g -Wall -Wno-unused 
FFLAGS  =  -fp-model precise -convert big_endian -assume byterecl -ftz
#LIBRARY = -lm

# Uncomment to include execution profiling information
#CFLAGS  = ${INCLUDES} -O3 -pg -Wall -Wno-unused
#LIBRARY = -lm

# Uncomment to debug memory problems using electric fence (man efence)
#CFLAGS  = ${INCLUDES} -g -Wall -Wno-unused
#LIBRARY = -lm -lefence -L/usr/local/lib

# -----------------------------------------------------------------------
# MOST USERS DO NOT NEED TO MODIFY BELOW THIS LINE
# -----------------------------------------------------------------------

.SUFFIXES: .c .F90 .o

.c.o:
	$(CC) -c $(CFLAGS) $<

.F90.o:
	$(F90) -c $(FFLAGS) $<

HDRS = ${VICPATH}/include/vic_def.h \
	${VICPATH}/include/vic_run.h \
	${DRIVERPATH}/include/vic_driver_classic.h \
	${DRIVERPATH}/include/mtclim_constants_vic.h \
	${DRIVERPATH}/include/mtclim_parameters_vic.h

OBJS = \
	${DRIVERPATH}/src/alloc_atmos.o \
	${DRIVERPATH}/src/alloc_cn.o \
	${DRIVERPATH}/src/alloc_veg_hist.o \
	${DRIVERPATH}/src/calc_air_temperature.o \
	${DRIVERPATH}/src/calc_longwave.o \
	${DRIVERPATH}/src/calc_root_fraction.o \
	${DRIVERPATH}/src/check_files.o \
	${DRIVERPATH}/src/check_state_file.o \
	${DRIVERPATH}/src/close_files.o \
	${DRIVERPATH}/src/cmd_proc.o \
	${DRIVERPATH}/src/compress_files.o \
	${DRIVERPATH}/src/compute_treeline.o \
	${DRIVERPATH}/src/display_current_settings.o \
	${DRIVERPATH}/src/free_all_vars.o \
	${DRIVERPATH}/src/free_vegcon.o \
	${DRIVERPATH}/src/get_dist.o \
	${DRIVERPATH}/src/get_force_type.o \
	${DRIVERPATH}/src/get_global_param.o \
	${DRIVERPATH}/src/initialize_atmos.o \
	${DRIVERPATH}/src/initialize_global.o \
	${DRIVERPATH}/src/initialize_model_state.o \
	${DRIVERPATH}/src/initialize_snow.o \
	${DRIVERPATH}/src/initialize_soil.o \
	${DRIVERPATH}/src/initialize_veg.o \
	${DRIVERPATH}/src/make_all_vars.o \
	${DRIVERPATH}/src/make_cell_data.o \
	${DRIVERPATH}/src/make_dmy.o \
	${DRIVERPATH}/src/make_energy_bal.o \
	${DRIVERPATH}/src/make_in_and_outfiles.o \
	${DRIVERPATH}/src/make_snow_data.o \
	${DRIVERPATH}/src/make_veg_var.o \
	${DRIVERPATH}/src/mtclim_vic.o \
	${DRIVERPATH}/src/mtclim_wrapper.o \
	${DRIVERPATH}/src/open_file.o \
	${DRIVERPATH}/src/open_state_file.o \
	${DRIVERPATH}/src/output_list_utils.o \
	${DRIVERPATH}/src/parse_output_info.o \
	${DRIVERPATH}/src/print_library.o \
	${DRIVERPATH}/src/read_atmos_data.o \
	${DRIVERPATH}/src/read_forcing_data.o \
	${DRIVERPATH}/src/read_initial_model_state.o \
	${DRIVERPATH}/src/read_lakeparam.o \
	${DRIVERPATH}/src/read_snowband.o \
	${DRIVERPATH}/src/read_soilparam.o \
	${DRIVERPATH}/src/read_veglib.o \
	${DRIVERPATH}/src/read_vegparam.o \
	${DRIVERPATH}/src/set_output_defaults.o \
	${DRIVERPATH}/src/vic_classic.o \
	${DRIVERPATH}/src/write_data.o \
	${DRIVERPATH}/src/write_forcing_file.o \
	${DRIVERPATH}/src/write_header.o \
	${DRIVERPATH}/src/write_model_state.o \
	${DRIVERPATH}/src/write_output.o \
	${VICPATH}/src/CalcAerodynamic.o \
	${VICPATH}/src/CalcBlowingSnow.o \
	${VICPATH}/src/IceEnergyBalance.o \
	${VICPATH}/src/SnowPackEnergyBalance.o \
	${VICPATH}/src/StabilityCorrection.o \
	${VICPATH}/src/advected_sensible_heat.o \
	${VICPATH}/src/arno_evap.o \
	${VICPATH}/src/calc_Nscale_factors.o \
	${VICPATH}/src/calc_atmos_energy_bal.o \
	${VICPATH}/src/calc_rainonly.o \
	${VICPATH}/src/calc_snow_coverage.o \
	${VICPATH}/src/calc_surf_energy_bal.o \
	${VICPATH}/src/calc_veg_params.o \
	${VICPATH}/src/calc_water_energy_balance_errors.o \
	${VICPATH}/src/canopy_assimilation.o \
	${VICPATH}/src/canopy_evap.o \
	${VICPATH}/src/compute_coszen.o \
	${VICPATH}/src/compute_pot_evap.o \
	${VICPATH}/src/compute_soil_resp.o \
	${VICPATH}/src/compute_zwt.o \
	${VICPATH}/src/correct_precip.o \
	${VICPATH}/src/estimate_T1.o \
	${VICPATH}/src/faparl.o \
	${VICPATH}/src/frozen_soil.o \
	${VICPATH}/src/func_atmos_energy_bal.o \
	${VICPATH}/src/func_atmos_moist_bal.o \
	${VICPATH}/src/func_canopy_energy_bal.o \
	${VICPATH}/src/func_surf_energy_bal.o \
	${VICPATH}/src/ice_melt.o \
	${VICPATH}/src/initialize_lake.o \
	${VICPATH}/src/lakes.eb.o \
	${VICPATH}/src/latent_heat_from_snow.o \
	${VICPATH}/src/massrelease.o \
	${VICPATH}/src/modify_Ksat.o \
	${VICPATH}/src/newt_raph_func_fast.o \
	${VICPATH}/src/nrerror.o \
	${VICPATH}/src/penman.o \
	${VICPATH}/src/photosynth.o \
	${VICPATH}/src/prepare_full_energy.o \
	${VICPATH}/src/put_data.o \
	${VICPATH}/src/root_brent.o \
	${VICPATH}/src/runoff.o \
	${VICPATH}/src/snow_intercept.o \
	${VICPATH}/src/snow_melt.o \
	${VICPATH}/src/snow_utility.o \
	${VICPATH}/src/soil_carbon_balance.o \
	${VICPATH}/src/soil_conduction.o \
	${VICPATH}/src/soil_thermal_eqn.o \
	${VICPATH}/src/solve_snow.o \
	${VICPATH}/src/surface_fluxes.o \
	${VICPATH}/src/svp.o \
	${VICPATH}/src/vic_run.o \
	${VICPATH}/src/water_energy_balance.o \
	${VICPATH}/src/water_under_ice.o \
	${VICPATH}/src/write_layer.o \
	${VICPATH}/src/write_vegvar.o \
	${VICPATH}/src/zero_output_list.o \
	${CNPATH}/src/VICCNInterface.o \
	${CNPATH}/src/pft2cov.o \
	${CNPATH}/src/shr_kind_mod.o \
	${CNPATH}/src/shr_const_mod.o \
	${CNPATH}/src/shr_orb_mod.o \
	${CNPATH}/src/clm_varpar.o \
	${CNPATH}/src/clm_varcon.o \
	${CNPATH}/src/nanMod.o \
	${CNPATH}/src/clmtype.o \
	${CNPATH}/src/clm_time_manager.o \
	${CNPATH}/src/pftvarcon.o \
	${CNPATH}/src/ndepStreamMod.o \
	${CNPATH}/src/SurfaceAlbedoMod.o \
	${CNPATH}/src/SurfaceRadiationMod.o \
	${CNPATH}/src/QSatMod.o \
	${CNPATH}/src/CanopyFluxesMod.o \
	${CNPATH}/src/subgridAveMod.o \
	${CNPATH}/src/CNAllocationMod.o \
	${CNPATH}/src/CNAnnualUpdateMod.o \
	${CNPATH}/src/CNBalanceCheckMod.o \
	${CNPATH}/src/CNCStateUpdate1Mod.o \
	${CNPATH}/src/CNCStateUpdate2Mod.o \
	${CNPATH}/src/CNCStateUpdate3Mod.o \
	${CNPATH}/src/CNDecompMod.o \
	${CNPATH}/src/CNFireMod.o \
	${CNPATH}/src/CNGRespMod.o \
	${CNPATH}/src/CNGapMortalityMod.o \
	${CNPATH}/src/CNMRespMod.o \
	${CNPATH}/src/CNNDynamicsMod.o \
	${CNPATH}/src/CNNStateUpdate1Mod.o \
	${CNPATH}/src/CNNStateUpdate2Mod.o \
	${CNPATH}/src/CNNStateUpdate3Mod.o \
	${CNPATH}/src/CNPhenologyMod.o \
	${CNPATH}/src/CNPrecisionControlMod.o \
	${CNPATH}/src/CNVegStructUpdateMod.o \
	${CNPATH}/src/CNSetValueMod.o \
	${CNPATH}/src/CNSummaryMod.o \
	${CNPATH}/src/CNWoodProductsMod.o \
	${CNPATH}/src/CNiniSpecial.o \
	${CNPATH}/src/clmtypeInitMod.o \
	${CNPATH}/src/iniTimeConst.o \
	${CNPATH}/src/CNEcosystemDynMod.o \
	${CNPATH}/src/CNiniTimeVar.o \
	${CNPATH}/src/clm_initializeMod.o \
	${CNPATH}/src/vic2clmtypeMod.o

OBJS2 = \
	alloc_atmos.o \
	alloc_cn.o \
	alloc_veg_hist.o \
	calc_air_temperature.o \
	calc_longwave.o \
	calc_root_fraction.o \
	check_files.o \
	check_state_file.o \
	close_files.o \
	cmd_proc.o \
	compress_files.o \
	compute_treeline.o \
	display_current_settings.o \
	free_all_vars.o \
	free_vegcon.o \
	get_dist.o \
	get_force_type.o \
	get_global_param.o \
	initialize_atmos.o \
	initialize_global.o \
	initialize_model_state.o \
	initialize_snow.o \
	initialize_soil.o \
	initialize_veg.o \
	make_all_vars.o \
	make_cell_data.o \
	make_dmy.o \
	make_energy_bal.o \
	make_in_and_outfiles.o \
	make_snow_data.o \
	make_veg_var.o \
	mtclim_vic.o \
	mtclim_wrapper.o \
	open_file.o \
	open_state_file.o \
	output_list_utils.o \
	parse_output_info.o \
	print_library.o \
	read_atmos_data.o \
	read_forcing_data.o \
	read_initial_model_state.o \
	read_lakeparam.o \
	read_snowband.o \
	read_soilparam.o \
	read_veglib.o \
	read_vegparam.o \
	set_output_defaults.o \
	vic_classic.o \
	write_data.o \
	write_forcing_file.o \
	write_header.o \
	write_model_state.o \
	write_output.o \
	CalcAerodynamic.o \
	CalcBlowingSnow.o \
	IceEnergyBalance.o \
	SnowPackEnergyBalance.o \
	StabilityCorrection.o \
	advected_sensible_heat.o \
	arno_evap.o \
	calc_Nscale_factors.o \
	calc_atmos_energy_bal.o \
	calc_rainonly.o \
	calc_snow_coverage.o \
	calc_surf_energy_bal.o \
	calc_veg_params.o \
	calc_water_energy_balance_errors.o \
	canopy_assimilation.o \
	canopy_evap.o \
	compute_coszen.o \
	compute_pot_evap.o \
	compute_soil_resp.o \
	compute_zwt.o \
	correct_precip.o \
	estimate_T1.o \
	faparl.o \
	frozen_soil.o \
	func_atmos_energy_bal.o \
	func_atmos_moist_bal.o \
	func_canopy_energy_bal.o \
	func_surf_energy_bal.o \
	ice_melt.o \
	initialize_lake.o \
	lakes.eb.o \
	latent_heat_from_snow.o \
	massrelease.o \
	modify_Ksat.o \
	newt_raph_func_fast.o \
	nrerror.o \
	penman.o \
	photosynth.o \
	prepare_full_energy.o \
	put_data.o \
	root_brent.o \
	runoff.o \
	snow_intercept.o \
	snow_melt.o \
	snow_utility.o \
	soil_carbon_balance.o \
	soil_conduction.o \
	soil_thermal_eqn.o \
	solve_snow.o \
	surface_fluxes.o \
	svp.o \
	vic_run.o \
	water_energy_balance.o \
	water_under_ice.o \
	write_layer.o \
	write_vegvar.o \
	zero_output_list.o \
	VICCNInterface.o \
	pft2cov.o \
	shr_kind_mod.o \
	shr_const_mod.o \
	shr_orb_mod.o \
	clm_varpar.o \
	clm_varcon.o \
	nanMod.o \
	clmtype.o \
	clm_time_manager.o \
	pftvarcon.o \
	ndepStreamMod.o \
	SurfaceAlbedoMod.o \
	SurfaceRadiationMod.o \
	QSatMod.o \
	CanopyFluxesMod.o \
	subgridAveMod.o \
	CNAllocationMod.o \
	CNAnnualUpdateMod.o \
	CNBalanceCheckMod.o \
	CNCStateUpdate1Mod.o \
	CNCStateUpdate2Mod.o \
	CNCStateUpdate3Mod.o \
	CNDecompMod.o \
	CNFireMod.o \
	CNGRespMod.o \
	CNGapMortalityMod.o \
	CNMRespMod.o \
	CNNDynamicsMod.o \
	CNNStateUpdate1Mod.o \
	CNNStateUpdate2Mod.o \
	CNNStateUpdate3Mod.o \
	CNPhenologyMod.o \
	CNPrecisionControlMod.o \
	CNVegStructUpdateMod.o \
	CNSetValueMod.o \
	CNSummaryMod.o \
	CNWoodProductsMod.o \
	CNiniSpecial.o \
	clmtypeInitMod.o \
	iniTimeConst.o \
	CNEcosystemDynMod.o \
	CNiniTimeVar.o \
	clm_initializeMod.o \
	vic2clmtypeMod.o

SRCS = $(OBJS:%.o=%.c) 

#$(SRCS):
#	co $@

#all:
#	make depend
#	make model

#default:
#	make depend
#	make model

#full:
#	make clean
#	make depend
#	make tags
#	make model

#clean::
#	/bin/rm -f ${DRIVERPATH}/src/*.o ${VICPATH}/src/*.o core log ${DRIVERPATH}/src/*~ ${VICPATH}/src/*~ vic_classic

#model: $(OBJS)
#	$(CC) -o vic_classic$(EXT) $(OBJS) $(CFLAGS) $(LIBRARY)
model:  $(OBJS)
	$(F90) -o vic_classic -nofor_main vic_classic$(EXT) $(OBJS2)

# -------------------------------------------------------------
# tags
# so we can find our way around
# -------------------------------------------------------------
#tags:	TAGS
#TAGS:	$(SRCS) $(HDRS)
#	etags $(SRCS) $(HDRS)
#clean::
#	\rm -f TAGS	       


# -------------------------------------------------------------
# depend
# -------------------------------------------------------------
#depend: .depend
#.depend: $(SRCS) $(HDRS)
#	$(CC) $(CFLAGS) -M $(SRCS) > $@

#clean::
#	\rm -f .depend	     
