# VIC image Makefile
# -----------------------------------------------------------------------

# -----------------------------------------------------------------------
# SET ENVIRONMENT-SPECIFIC OPTIONS HERE
# -----------------------------------------------------------------------

# VIC DRIVER PATH (include and src are subdirs of this)
DRIVERPATH = ..

# VIC RUN PATH
VICPATH = ../../../vic_run

# NETCDF include .. change the path depending on travis or other location
ifeq (true, ${TRAVIS})
NETCDFPATH = /usr
else ifeq (hydro-c1, $(shell hostname))
NETCDFPATH = /usr/local/netcdf-4.3.0+gcc-4.7.2
else
NETCDFPATH = /opt/local
endif

# Set SHELL = your shell here
SHELL = /bin/bash

# Set CC = your compiler here
CC = gcc

# set includes
INCLUDES = -I ${DRIVERPATH}/include -I ${VICPATH}/include -I ${DRIVERPATH}/../classic/include -I ${NETCDFPATH}/include

# Uncomment for normal optimized code flags (fastest run option)
#CFLAGS  = -O3 -Wall -Wno-unused
LIBRARY = -lm -L${NETCDFPATH}/lib -lnetcdf

# Uncomment to include debugging information
CFLAGS  =  ${INCLUDES} -ggdb -O0 -Wall 
# -Wno-unused
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

HDRS = ${VICPATH}/include/vic_def.h \
	${VICPATH}/include/vic_run.h \
	${DRIVERPATH}/include/vic_driver_image.h \
	${DRIVERPATH}/../classic/include/vic_driver_classic.h \
	${DRIVERPATH}/../classic/include/mtclim_constants_vic.h

OBJS = \
	${DRIVERPATH}/src/alloc_atmos.o \
	${DRIVERPATH}/src/alloc_veg_hist.o \
	${DRIVERPATH}/src/calc_root_fraction.o \
	${DRIVERPATH}/src/cmd_proc.o \
	${DRIVERPATH}/src/display_current_settings.o \
	${DRIVERPATH}/src/free_all_vars.o \
	${DRIVERPATH}/src/get_global_domain.o \
	${DRIVERPATH}/src/get_global_param.o \
	${DRIVERPATH}/src/get_nc_dimension.o \
	${DRIVERPATH}/src/get_nc_field.o \
	${DRIVERPATH}/src/init_library.o \
	${DRIVERPATH}/src/initialize_energy.o \
	${DRIVERPATH}/src/initialize_global.o \
	${DRIVERPATH}/src/initialize_model_state.o \
	${DRIVERPATH}/src/initialize_snow.o \
	${DRIVERPATH}/src/initialize_soil.o \
	${DRIVERPATH}/src/initialize_veg.o \
	${DRIVERPATH}/src/make_all_vars.o \
	${DRIVERPATH}/src/make_cell_data.o \
	${DRIVERPATH}/src/make_dmy.o \
	${DRIVERPATH}/src/make_energy_bal.o \
	${DRIVERPATH}/src/make_snow_data.o \
	${DRIVERPATH}/src/make_veg_var.o \
	${DRIVERPATH}/src/open_file.o \
	${DRIVERPATH}/src/output_list_utils.o \
	${DRIVERPATH}/src/parse_output_info.o \
	${DRIVERPATH}/src/print_library.o \
	${DRIVERPATH}/src/put_nc_field.o \
	${DRIVERPATH}/src/soil_moisture_from_water_table.o \
	${DRIVERPATH}/src/vic_alloc.o \
	${DRIVERPATH}/src/vic_finalize.o \
	${DRIVERPATH}/src/vic_force.o \
	${DRIVERPATH}/src/vic_image.o \
	${DRIVERPATH}/src/vic_image_run.o \
	${DRIVERPATH}/src/vic_init.o \
	${DRIVERPATH}/src/vic_init_output.o \
	${DRIVERPATH}/src/vic_nc_info.o \
	${DRIVERPATH}/src/vic_restore.o \
	${DRIVERPATH}/src/vic_start.o \
	${DRIVERPATH}/src/vic_store.o \
	${DRIVERPATH}/src/vic_write.o \
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
	${VICPATH}/src/zero_output_list.o

SRCS = $(OBJS:%.o=%.c)

#$(SRCS):
#	co $@

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
	/bin/rm -f ${DRIVERPATH}/src/*.o ${VICPATH}/src/*.o core log ${DRIVERPATH}/src/*~ ${VICPATH}/src/*~ vic_image

model: $(OBJS)
	$(CC) -o vic_image$(EXT) $(OBJS) $(CFLAGS) $(LIBRARY)

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
.depend: $(SRCS) $(HDRS)
	$(CC) $(CFLAGS) -M $(SRCS) > $@

clean::
	\rm -f .depend
