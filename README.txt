# README.txt - Release Notes
#
# $Id$
#------------------------------------------------------------------------

August 15, 2000: VIC release 4.0.3

This release fixes a problem with the implementation of va_arg that 
causes run time errors on some systems.  Previous releases of the code 
worked correctly on the LINUX and freeBSD systems where it was tested. 
However on some systems (including Sun Ultra-2s) character variables 
passed with va_arg are changed into integers so reading a character from 
the argument list does not produce the value sent to the routine.  The
character flags used by VIC to indicate if there is snow present and if
the frozen soil algorithm has been activated have now been converted to 
integers, which should make the va_arg call work on all systems.

Also fixed in this release was a check in dist_prec.c to see if it is 
still raining which actually used the memory address of the 
precipitation variable rather than the daily value in the check.

MODIFIED FILES:
	read_atmos_data.c	-	Fixed input file time step check
	write_forcing_files.c	-	Added free statements for pointers
	calc_surf_energy_bal.c	-	Converted char flags to int
	dist_prec.c		-	Fixed logical statement error
	frozen_soil.c		-	Converted char flags to int
	func_surf_energy_bal.c	-	Converted char flags to int
	initialize_atmos.c	-	Added flag for output forcing
	vicNl.h			-	Converted char flags to int
	vicNl_def.h		-	Converted char flags to int

July 19, 2000: VIC release 4.0.2

Two new pre-processor options have been added to VIC as well as minor 
modifications to two subroutines.  

If set to TRUE the NO_REWIND pre-processor option stops the VIC model from
rewinding the soil and vegetation parameter input files for each new grid
cell.  This reduces run times but requires that all input files are in the
same order as the soil parameter file.  

If set TRUE the OUTPUT_FORCE pre-processor option blocks off the main
model and only reads the provided forcing files. Once VIC has estimated
the missing forcing parameters the full forcing data set for the defined
simulation period is output to a series of gridded forcing files.  The
gridded forcing files are written to the RESULTS directory defined in the
global control file with the prefix "full_data_".  The new files are in
Binary or ASCII depending on the setting of BINARY_OUTPUT.

The error messages in get_global_param.c have been modified so that the
correct file is referenced when telling the user to change values found in
the model source code.

In read_soilparam.c, the soil parameters are defined only if the current
grid cell is run, otherwise the line in the file is skipped and soil_con
is retruned without new data values.

May 30, 2000: VIC release 4.0.1

Increased use of the released VIC model code has lead to the
discovery of a couple of minor bugs.  This release fixes those bugs as
well as introducing a improved precipitation correction algorithm based on 
Yang et al. 1998.  Unless you have encountered these problems or are
trying to correct precipitation undercatch due to wind, in the VIC
model, your results will not be impacted by these fixes.

MODIFIED FILES:
    correct_precip.c       - changed to WMO correction equation for
                             NWS 8" standard gauge.
    full_energy.c          - modified to handle WMO correction.
    initialize_atmos.c     - modified to handle WMO correction.
                             fixed error in estimating minimum daily 
                             temperature from sub-daily temperatures.
    make_in_and_outfiles.c - removed line that opened the state file
                             again for each new grid cell.
    open_state_file.c      - modified comments.
    put_data.c             - modified to handle WMO correction.
    snow_utility.c         - cleaned source code.
    solve_snow.c           - modified to handle WMO correction.
    surface_fluxes.c       - modified to handle WMO correction.
    vicNl.h                - modified to handle WMO correction.
    vicNl_def.h            - modified to handle WMO correction.

REFERENCE:

    Yang, D., B. E. Goodison, J. R. Metcalfe, V. S. Golubev, R.
    Bates, T. Pangburn, and C. L. Hanson, Accuracy of NWS 8" Standard 
    Nonrecording Precipitation Gauge: Results and Application of WMO 
    Intercomparison, J. Atmos. Oceanic Tech., 15, 54-68, 1998.
             

Date: May 16, 2000

From: Keith Cherkauer

Topic: Release of VIC 4.0.0

The code for VIC release 4.0.0 has undergone several months of tests (as
VIC release 3.3.0 Beta) and has now been deemed ready for release to the
general public.  This document is designed to provide information
concerning changes in the model between the last release version (3.2.1)
and the current version.

There is no formal users manual, information about how to use the current
version can be found at
http://www.hydro.washington.edu/Lettenmaier/Models/VIC/VIChome.html.  
Information about the basic model design can be found in Liang, et al.
(1994), while the rewrite of the source code as well as the addition of
cold season processes is described in Cherkauer and Lettenmaier (1999).

The VIC macroscale hydrologic model was developed as a research tool.  As
such it is the users responsibility to verify that it is working correctly
when applied to new situations.  When using the VIC model reference should
be made to Liang, et al. (1994) and Cherkauer and Lettenmaier (1999) as
well as an acknowledgment that the code was received from the University
of Washington.  Other important papers relating to the development of the
VIC model are included on the home page and in the source code.

Possible bugs in the code should be reported to
vicadmin@hydro.washington.edu.  ALWAYS CHECK YOUR INPUT FILES!  Most
"bugs" are actually caused by trying to run the model with bad parameters
or forcing data.  The VIC model will run limited checks to find common
major errors but in most cases it will attempt to run with the bad values.  
If after checking all of your input data you still believe you have found
a bug in the model, send an e-mail including the complete output from the
model as well as a description of the problem and the files necessary to
run the model to recreate the code (if files are large please put a
compressed tar file in
ftp://ftp.ce.washington.edu/pub/HYDRO/vicadmin/TEMP).  Outdated and
modified versions of the code are the responsibility of the user to debug.  
Modifications made to the code, which may improve the general model
performance, may be submitted to vicadmin@hydro.washington.edu for
possible inclusion in future versions of the code.

VIC release 4.0.0 represents a major change to the source code from
version 3.2.1.  It is strongly recommended that if you were using version
3.2.1 or earlier versions that you update with a complete copy of the new
code.

Major changes from release version 3.2.1 to 4.0.0:

- Radiation Estimation Update: The routines to estimate shortwave and
longwave radiation as well as vapor pressure from daily minimum and
maximum temperature and precipitation have been updated to correspond to
the algorithm described by Thornton and Running (1999).  These routines
provide significantly improved radiation estimates especially in regions
outside the continental United States.

- Model Core Update: The core of the VIC model was rewritten so that all
modes of the model (water balance, energy balance, etc.) make use of the
same model code.  This makes it easier to modify the model and have
modifications apply to sll modes, it also allows the model to run with new
combinations of algorithms (i.e. full energy balance mode with the finite
difference ground heat flux solution).

- Soil Moisture Transport Update: The frozen and thawed sub-layers added
to the model for the original frozen soil algorithm have been removed.  
This makes the soil drainage routine cleaner and faster.  Frozen soils now
estimate full layer ice contents from the ice content at each soil thermal
node.  Without being confined by sub-layers, the frozen soil algorithm can
now be applied to regions of permafrost.

- Forcing File Control Added: Version 4.0.0 moves controls of the forcing
file format and data types into the global control file.  The model can
now handle most ASCII column and short int Binary files without writing
new subroutines and recompiling the source code.

- Pre-processor Options Added: There are now more option flags in the
source code headers to control which parts of the model are in fact
compiled.  This allows the model functionality to be adjusted without the
addition of computationally intensive conditional switching statements.

- Model State File: With the release of version 4.0.0 separate snow and
soil initialization files have been combined into a single model state
file.  The state file can be created outside the model for starting
simulations with prescribed initial conditions, or the model state can be
saved by VIC at a specified date.  Note that currently there will be small
differences between a full and a warm started simulation because radiation
and vapor pressure are estimated using forcing data from the simulation
period, not from the full dataset included in the forcing file.  It also
does not store wet and dry fraction information, when running with
distributed precipitation the model is restarted using average soil
conditions.

References:

Liang, X., D. P. Lettenmaier, E. F. Wood, and S. J. Burges, A simple
hydrologically based model of land surface water and energy fluxes for
GSMs, J. Geophys. Res., 99(D7), 14,415-14,428, 1994.

Cherkauer, K. A., and D. P. Lettenmaier, Hydrologic effects of frozen
soils in the upper Mississippi River basin, J. Geophys. Res., 104(D16),
19,599-19,610, 1999.

Thornton, P.E. and S.W. Running, An improved algorithm for estimating
incident daily solar radiation from measurements of temperature, humidity,
and precipitation, Ag. For. Met., 93, 211-228, 1999.

File List:

CalcAerodynamic.c               modify_Ksat.c
Makefile                        mtclim42_vic.c
SnowPackEnergyBalance.c         mtclim42_vic.h
StabilityCorrection.c           mtclim42_wrapper.c
alloc_atmos.c                   nrerror.c
arno_evap.c                     open_debug.c
calc_air_temperature.c          open_file.c
calc_cloud_cover_fraction.c     open_state_file.c
calc_longwave.c                 penman.c
calc_rainonly.c                 prepare_full_energy.c
calc_root_fraction.c            put_data.c
calc_surf_energy_bal.c          read_arcinfo_ascii.c
calc_veg_params.c               read_atmos_data.c
canopy_evap.c                   read_forcing_data.c
check_files.c                   read_initial_model_state.c
check_state_file.c              read_snowband.c
close_files.c                   read_soilparam.c
cmd_proc.c                      read_soilparam_arc.c
compress_files.c                read_veglib.c
compute_dz.c                    read_vegparam.c
correct_precip.c                redistribute_during_storm.c
dist_prec.c                     root_brent.c
estimate_T1.c                   runoff.c
free_dist_prcp.c                snow.h
free_vegcon.c                   snow_intercept.c
frozen_soil.c                   snow_melt.c
full_energy.c                   snow_utility.c
func_surf_energy_bal.c          soil_conduction.c
get_force_type.c                soil_thermal_eqn.c
get_global_param.c              solve_snow.c
global.h                        store_moisture_for_debug.c
initialize_atmos.c              surface_fluxes.c
initialize_global.c             svp.c
initialize_model_state.c        user_def.h
initialize_new_storm.c          vicNl.c
initialize_snow.c               vicNl.h
initialize_soil.c               vicNl_def.h
initialize_veg.c                vicerror.c
make_cell_data.c                write_atmosdata.c
make_dist_prcp.c                write_data.c
make_dmy.c                      write_debug.c
make_energy_bal.c               write_layer.c
make_in_and_outfiles.c          write_model_state.c
make_snow_data.c                write_soilparam.c
make_veg_var.c                  write_vegparam.c
massrelease.c                   write_vegvar.c
