# VIC Makefile
# Modifications:
# 27-May-2003 Replaced read_vegparam by read_vegparam_LAI       KAC
# 12-Nov-2003 Added "make depend" to the "all" and "default" options.
#             This way, if a user always types "make", the user is
#             guaranteed to have a .depend file and therefore *.o will
#             always be recompiled whenever a .h file is updated.  The
#             user can override this behavior by typing "make model",
#             which doesn't invoke "make depend".		TJB
# 2006-Sep-11 Changes for flexible output configuration. TJB
#             Added the following files:
#               calc_water_energy_balance_errors.c
#               output_list_utils.c
#               parse_output_info.c
#               set_output_defaults.c
#             
# -----------------------------------------------------------------------

# -----------------------------------------------------------------------
# SET ENVIRONMENT-SPECIFIC OPTIONS HERE
# -----------------------------------------------------------------------

# Set SHELL = your shell here
SHELL = /bin/csh

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

HDRS = vicNl.h vicNl_def.h global.h snow.h user_def.h mtclim42_vic.h

OBJS =  CalcAerodynamic.o SnowPackEnergyBalance.o StabilityCorrection.o \
	alloc_atmos.o arno_evap.o calc_air_temperature.o \
	calc_cloud_cover_fraction.o calc_longwave.o calc_rainonly.o \
	calc_root_fraction.o calc_surf_energy_bal.o calc_veg_params.o \
	calc_water_energy_balance_errors.o canopy_evap.o check_files.o \
	check_state_file.o close_files.o \
	cmd_proc.o compress_files.o compute_dz.o compute_treeline.o \
	correct_precip.o display_current_settings.o dist_prec.o \
	estimate_T1.o free_dist_prcp.o free_vegcon.o frozen_soil.o \
	full_energy.o func_surf_energy_bal.o \
	get_force_type.o get_global_param.o initialize_atmos.o \
	initialize_model_state.o initialize_global.o \
	initialize_new_storm.o initialize_snow.o initialize_soil.o \
	initialize_veg.o make_cell_data.o make_dist_prcp.o make_dmy.o \
	make_energy_bal.o make_in_and_outfiles.o make_snow_data.o \
	make_veg_var.o massrelease.o modify_Ksat.o mtclim42_vic.o \
	mtclim42_wrapper.o nrerror.o open_debug.o open_file.o \
	open_state_file.o output_list_utils.o parse_output_info.o \
	penman.o prepare_full_energy.o put_data.o \
	read_arcinfo_ascii.o read_atmos_data.o read_forcing_data.o \
	read_initial_model_state.o read_snowband.o \
	read_soilparam.o read_soilparam_arc.o read_veglib.o read_vegparam.o \
	redistribute_during_storm.o root_brent.o runoff.o \
	set_output_defaults.o snow_intercept.o snow_melt.o \
	snow_utility.o soil_conduction.o soil_thermal_eqn.o \
	solve_snow.o store_moisture_for_debug.o surface_fluxes.o svp.o \
	vicNl.o vicerror.o write_atmosdata.o write_data.o write_debug.o \
	write_forcing_file.o write_layer.o write_model_state.o \
	write_soilparam.o write_vegparam.o write_vegvar.o 

SRCS = $(OBJS:%.o=%.c) 

all:
	make depend
	make model

default:
	make depend
	make model

full:
	make clean
	make depend
	make tags
	make model

clean::
	/bin/rm -f *.o core log *~

model: $(OBJS)
	$(CC) -o vicNl $(OBJS) $(CFLAGS) $(LIBRARY)

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
