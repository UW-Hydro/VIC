 ##############################################################################
 # @section DESCRIPTION
 #
 # VIC Image Driver Makefile
 #
 # @section LICENSE
 #
 # The Variable Infiltration Capacity (VIC) macroscale hydrological model
 # Copyright (C) 2014 The Land Surface Hydrology Group, Department of Civil
 # and Environmental Engineering, University of Washington.
 #
 # The VIC model is free software; you can redistribute it and/or
 # modify it under the terms of the GNU General Public License
 # as published by the Free Software Foundation; either version 2
 # of the License, or (at your option) any later version.
 #
 # This program is distributed in the hope that it will be useful,
 # but WITHOUT ANY WARRANTY; without even the implied warranty of
 # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 # GNU General Public License for more details.
 #
 # You should have received a copy of the GNU General Public License along with
 # this program; if not, write to the Free Software Foundation, Inc.,
 # 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 ##############################################################################

# -----------------------------------------------------------------------
# SET ENVIRONMENT-SPECIFIC OPTIONS HERE
# -----------------------------------------------------------------------

# VIC DRIVER PATH (include and src are subdirs of this)
DRIVERPATH = .

# VIC RUN PATH
SHAREDPATH = ../shared

# VIC RUN PATH
VICPATH = ../../vic_run

# NETCDF include .. change the path depending on travis or other location
ifeq (true, ${TRAVIS})
NETCDFPATH = /usr
else ifeq (hydro-c1, $(shell hostname))
NETCDFPATH = /usr/local/netcdf-4.3.0+gcc-4.7.2
else
NETCDFPATH = /opt/local
endif

# MPI include .. change the path depending on travis or other location
ifeq (true, ${TRAVIS})
MPIPATH = /usr
else ifeq (hydro-c1, $(shell hostname))
MPIPATH = /usr/local/netcdf-4.3.0+gcc-4.7.2
else
MPIPATH = /opt/local
endif

# Set SHELL = your shell here
SHELL = /bin/bash

# Set CC = your compiler here
ifeq (true, ${TRAVIS})
CC = mpicc
else
CC = mpicc-mpich-mp
endif

# set includes
INCLUDES = -I ${DRIVERPATH}/include -I ${VICPATH}/include -I ${SHAREDPATH}/include -I ${NETCDFPATH}/include -I ${MPIPATH}/include

# Uncomment for normal optimized code flags (fastest run option)
#CFLAGS  = -O3 -Wall -Wno-unused
LIBRARY = -lm -L${NETCDFPATH}/lib -lnetcdf -L${MPIPATH}/lib -lmpi

# Uncomment to include debugging information
CFLAGS  =  ${INCLUDES} -ggdb -O0 -Wall -Wextra
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
	${SHAREDPATH}/include/vic_driver_shared.h \
	${VICPATH}/include/vic_physical_constants.h \
	${VICPATH}/include/vic_log.h

OBJS = \
	${DRIVERPATH}/src/alloc_atmos.o \
	${DRIVERPATH}/src/alloc_veg_hist.o \
	${DRIVERPATH}/src/display_current_settings.o \
	${DRIVERPATH}/src/get_global_domain.o \
	${DRIVERPATH}/src/get_global_param.o \
	${DRIVERPATH}/src/get_nc_dimension.o \
	${DRIVERPATH}/src/get_nc_field.o \
	${DRIVERPATH}/src/init_library.o \
	${DRIVERPATH}/src/initialize_energy.o \
	${DRIVERPATH}/src/initialize_model_state.o \
	${DRIVERPATH}/src/output_list_utils.o \
	${DRIVERPATH}/src/parse_output_info.o \
	${DRIVERPATH}/src/print_library_image.o \
	${DRIVERPATH}/src/put_nc_field.o \
	${DRIVERPATH}/src/vic_alloc.o \
	${DRIVERPATH}/src/vic_finalize.o \
	${DRIVERPATH}/src/vic_force.o \
	${DRIVERPATH}/src/vic_image.o \
	${DRIVERPATH}/src/vic_image_run.o \
	${DRIVERPATH}/src/vic_init.o \
	${DRIVERPATH}/src/vic_init_output.o \
	${DRIVERPATH}/src/vic_mpi_support.o \
	${DRIVERPATH}/src/vic_nc_info.o \
	${DRIVERPATH}/src/vic_restore.o \
	${DRIVERPATH}/src/vic_start.o \
	${DRIVERPATH}/src/vic_store.o \
	${DRIVERPATH}/src/vic_write.o \
	${SHAREDPATH}/src/calc_root_fraction.o \
	${SHAREDPATH}/src/cmd_proc.o \
	${SHAREDPATH}/src/compress_files.o \
	${SHAREDPATH}/src/compute_treeline.o \
	${SHAREDPATH}/src/free_all_vars.o \
	${SHAREDPATH}/src/free_vegcon.o \
	${SHAREDPATH}/src/get_dist.o \
	${SHAREDPATH}/src/get_parameters.o \
	${SHAREDPATH}/src/initialize_files.o \
	${SHAREDPATH}/src/initialize_global.o \
	${SHAREDPATH}/src/initialize_options.o \
	${SHAREDPATH}/src/initialize_parameters.o \
	${SHAREDPATH}/src/initialize_snow.o \
	${SHAREDPATH}/src/initialize_soil.o \
	${SHAREDPATH}/src/initialize_veg.o \
	${SHAREDPATH}/src/make_all_vars.o \
	${SHAREDPATH}/src/make_cell_data.o \
	${SHAREDPATH}/src/make_dmy.o \
	${SHAREDPATH}/src/make_energy_bal.o \
	${SHAREDPATH}/src/make_snow_data.o \
	${SHAREDPATH}/src/make_veg_var.o \
	${SHAREDPATH}/src/open_file.o \
	${SHAREDPATH}/src/print_library_shared.o \
	${SHAREDPATH}/src/soil_moisture_from_water_table.o \
	${SHAREDPATH}/src/vic_log.o \
	${SHAREDPATH}/src/vic_time.o \
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
	${VICPATH}/src/penman.o \
	${VICPATH}/src/photosynth.o \
	${VICPATH}/src/physics.o \
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
	/bin/rm -f ${DRIVERPATH}/src/*.o ${SHAREDPATH}/src/*.o  ${VICPATH}/src/*.o core log ${DRIVERPATH}/src/*~ ${SHAREDPATH}/src/*.~ ${VICPATH}/src/*~ vic_image

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
