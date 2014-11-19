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
# 2014-Mar-25 Removed calc_cloud_cover_fraction.c                                BN
# 2014-Mar-26 Removed compute_dz.c                                               BN
# 2014-Mar-26 Removed files with unused functions:
#             write_snowparam.c
#             write_soilparam.c
#             write_vegparam.c                                                   BN
# 2014-Mar-28 Removed DIST_PRCP option, and the files:
#             dist_prec.c
#             free_dist_prcp.c
#             make_dist_prcp.c
#             initialize_new_storm.c
#             redistribute_during_storm.c					TJB
# 2014-Apr-25 Added alloc_veg_hist.c.						TJB
#
# $Id$
#
# -----------------------------------------------------------------------

# -----------------------------------------------------------------------
# SET ENVIRONMENT-SPECIFIC OPTIONS HERE
# -----------------------------------------------------------------------

# Set SHELL = your shell here
SHELL = /bin/bash

# Set CC = your compiler here
CC = gcc

# Uncomment for normal optimized code flags (fastest run option)
#CFLAGS  = -I. -O3 -Wall -Wno-unused
LIBRARY = -lm

# Uncomment to include debugging information
CFLAGS  = -I. -g -Wall -Wno-unused
#LIBRARY = -lm

# Uncomment to include execution profiling information
#CFLAGS  = -I. -O3 -pg -Wall -Wno-unused
#LIBRARY = -lm

# Uncomment to debug memory problems using electric fence (man efence)
#CFLAGS  = -I. -g -Wall -Wno-unused
#LIBRARY = -lm -lefence -L/usr/local/lib

# -----------------------------------------------------------------------
# MOST USERS DO NOT NEED TO MODIFY BELOW THIS LINE
# -----------------------------------------------------------------------

HDRS = vicNl.h vicNl_def.h global.h snow.h mtclim_constants_vic.h mtclim_parameters_vic.h LAKE.h

OBJS =  CalcAerodynamic.o CalcBlowingSnow.o SnowPackEnergyBalance.o \
        StabilityCorrection.o advected_sensible_heat.o alloc_atmos.o \
        alloc_veg_hist.o arno_evap.o calc_air_temperature.o \
	calc_atmos_energy_bal.o calc_longwave.o calc_Nscale_factors.o \
	calc_rainonly.o calc_root_fraction.o calc_snow_coverage.o \
	calc_surf_energy_bal.o calc_veg_params.o \
	calc_water_energy_balance_errors.o canopy_assimilation.o canopy_evap.o \
	check_files.o check_state_file.o close_files.o cmd_proc.o \
	compress_files.o compute_coszen.o compute_pot_evap.o \
	compute_soil_resp.o compute_treeline.o compute_zwt.o correct_precip.o \
	display_current_settings.o estimate_T1.o faparl.o free_all_vars.o \
	free_vegcon.o frozen_soil.o full_energy.o func_atmos_energy_bal.o \
	func_atmos_moist_bal.o func_canopy_energy_bal.o \
	func_surf_energy_bal.o get_dist.o get_force_type.o get_global_param.o \
	initialize_atmos.o initialize_model_state.o \
	initialize_global.o initialize_snow.o \
	initialize_soil.o initialize_veg.o latent_heat_from_snow.o \
	make_cell_data.o make_all_vars.o make_dmy.o make_energy_bal.o \
	make_in_and_outfiles.o make_snow_data.o make_veg_var.o massrelease.o \
	modify_Ksat.o mtclim_vic.o mtclim_wrapper.o newt_raph_func_fast.o \
	nrerror.o open_file.o open_state_file.o \
	output_list_utils.o parse_output_info.o penman.o photosynth.o \
	prepare_full_energy.o print_library.o put_data.o \
	read_atmos_data.o read_forcing_data.o read_initial_model_state.o \
	read_snowband.o read_soilparam.o read_veglib.o \
	read_vegparam.o root_brent.o runoff.o \
	set_output_defaults.o snow_intercept.o snow_melt.o \
	snow_utility.o soil_carbon_balance.o soil_conduction.o \
	soil_thermal_eqn.o solve_snow.o \
	surface_fluxes.o svp.o vicNl.o vicerror.o \
	write_data.o write_forcing_file.o write_header.o write_layer.o \
	write_model_state.o write_vegvar.o lakes.eb.o initialize_lake.o \
	read_lakeparam.o ice_melt.o IceEnergyBalance.o water_energy_balance.o \
	water_under_ice.o

SRCS = $(OBJS:%.o=%.c) 

#$(SRCS):
#	co $@

all:
	make depend
	make model

disagg:
	sed -i.bak 's/OUTPUT_FORCE FALSE/OUTPUT_FORCE TRUE/' user_def.h
	make clean
	make depend
	make vicDisagg
	sed -i.bak 's/OUTPUT_FORCE TRUE/OUTPUT_FORCE FALSE/' user_def.h
	make clean
	make depend

default:
	make depend
	make model

full:
	make clean
	make depend
	make tags
	make model
	make disagg

clean::
	/bin/rm -f *.o core log *~

model: $(OBJS)
	$(CC) -o vicNl$(EXT) $(OBJS) $(CFLAGS) $(LIBRARY)

vicDisagg: $(OBJS)
	$(CC) -o vicDisagg $(OBJS) $(CFLAGS) $(LIBRARY)

# -------------------------------------------------------------
# tags
# so we can find our way around
# -------------------------------------------------------------
tags:	TAGS
TAGS:	$(SRCS) $(HDRS)
	etags $(SRCS) $(HDRS)
clean::
	\rm -f TAGS	       


# -------------------------------------------------------------
# depend
# -------------------------------------------------------------
depend: .depend
.depend:	$(SRCS) $(HDRS)
	$(CC) $(CFLAGS) -M $(SRCS) > $@

clean::
	\rm -f .depend	     
