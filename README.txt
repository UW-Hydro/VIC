# README.txt - Release notes
# $Id$
#-----------------------------------------------------------------------

Modifications:

	cmd_proc.c, display_current_settings.c, get_global_param.c,
	global.h, vicNl.c, vicNl.h, vicNl_def.h:
	    Added -o option, to display current run-time and
	    compile-time options.  Updated version string.	TJB
	Makefile:
	    Now automatically re-compiles all dependent files when
	    any .h files are updated.				TJB
	open_file.c:
	    Now announces when it opens a binary file for reading.
	write_data.c:
	    Corrected typo in output format for fluxes file.	TJB
	calc_rainonly.c:
	    (found by Justin Sheffield at Princeton)
	    Changed test
		else if(air_temp > MAX_SNOW_TEMP)
	    to
		else if(air_temp >= MAX_SNOW_TEMP)
	    to fix situation in which, if air_temp = MAX_SNOW_TEMP,
	    rainfall (rainonly) was set to 0 and snowfall was set
	    to 100% of precip, causing function to fail.	TJB
	compute_dz.c, initialize_atmos.c, read_soilparam.c,
	read_soilparam_arc.c:
	    (found by Justin Sheffield at Princeton)
	    Replaced rint(something) with (float)(int)(something + 0.5)
	    to handle rounding without resorting to rint(), which
	    isn't supported on all platforms.			TJB
	    Added print statement for current cell number.	TJB
	func_surf_energy_bal.c:
	    (found by Justin Sheffield at Princeton)
	    Added check that both FS_ACTIVE and FROZEN_SOIL are true
	    before adjusting *deltaH.  This fixes a bug caused when
	    FS_ACTIVE was false and FROZEN_SOIL was true, in which
	    non-zero ice content was calculated at the beginning of
	    the time step but always assumed zero at the point of
	    *deltaH adjustment, leading to *deltaH always being
	    calculated as if ice had melted, leading to lower soil
	    temperatures and snow packs that never melted.	TJB
	initialize_global.c:
	    (found by Justin Sheffield at Princeton)
	    Initialize ARC_SOIL, COMPRESS, and ARNO_PARAMS to FALSE.
	    Also changed limit on loop over forcing types from
	    hard-coded 17 to variable N_FORCING_TYPES.		TJB
	initialize_model_state.c:
	    (found by Justin Sheffield at Princeton)
	    Initialize soil_con->dz_node[Nnodes] to 0.0, since it is
	    accessed in set_node_parameters().			TJB
	open_debug.c:
	    (found by Justin Sheffield at Princeton)
	    Initialize debug_store_moist array when debug.PRT_MOIST
	    is true (as well as under the other previously-defined
	    conditions).					TJB
	open_state_file.c:
	    (Port from 4.1.0)
	    Modified the statefile name to contain year, month, day
	    rather than day, month, year.  This makes it consistent
	    with the planned release of 4.1.0.			TJB
	penman.c:
	    (found by Justin Sheffield at Princeton)
	    Changed
		if (vpd > 0.0 && evap < 0.0)
	    to
		if (vpd >= 0.0 && evap < 0.0)
	    to correctly handle evap when vpd == 0.0.		TJB
	read_atmos_data.c:
	    (Port from 4.1.0)
	    Replaced NF with global_param.dt in condition checking
	    whether forcing file contains enough records to cover
	    the time range of the simulation.			TJB
	read_initial_model_state.c:
	    (Port from 4.1.0)
	    Added check to verify that the sum of the defined nodes
	    equals the damping depth.				TJB
	read_soilparam.c, read_soilparam_arc.c:
	    (fix by Chunmei Zhu and Alan Hamlet)
	    Added check to make sure that wilting point
	    (porosity*Wpwp_FRACT) is greater than residual
	    moisture.						TJB
	read_soilparam_arc.c:
	    Modified to handle ARNO parameters.			TJB
	runoff.c:
	    (fix by Chunmei Zhu and Alan Hamlet)
	    Added check to make sure baseflow doesn't exceed
	    difference between liquid moisture and field
	    capacity.						TJB
	    Changed block that handles baseflow when soil moisture
	    drops below residual moisture.  Now, the block is only
	    entered if baseflow > 0 and soil moisture < residual,
	    and the amount of water taken out of baseflow and given
	    to the soil cannot exceed baseflow.  In addition, error
	    messages are no longer printed, since it isn't an error
	    to be in that block.				TJB
	soil_conduction.c:
	    (Port from 4.1.0)
	    set_node_parameters(): Modified to correct differences
	    between calculations to determine maximum node moisture
	    and node moisture, so that nodes on the boundary between
	    soil layers are computed the same way for both.	TJB
	    distribute_node_moisture_properties(): Modified to check
	    that node soil moisture is less than or equal to maximum
	    node soil moisture, otherwise an error is printed to the
	    screen and the model exits.				TJB
	surface_fluxes.c:
	    (found by Justin Sheffield at Princeton)
	    Fixed initialization of canopyevap to initialize for every
	    value of dist, rather than just dist 0.		TJB
	write_atmosdata.c:
	    (found by Justin Sheffield at Princeton)
	    No longer close the debug file, since the next cell
	    must write to it.					TJB

December 1, 2003: VIC release 4.0.4

This release covers fixes for a small number of additional bugs discovered
since the release of 4.0.4 beta r3, plus 2 new features.

New Features:

      - Typing 'vicNl -v' on the command line now displays the version
	of VIC you are running.  This provides a fast, simple way of
	determining the version of VIC without running a simulation.

      - The parameters {Ds, Dsmax, Ws and c} in the soil parameter file
	may now be replaced by the set of ARNO baseflow parameters
	{d1, d2, d3, d4}.  To enable VIC to read the ARNO baseflow
	parameters, you must add the following line to the global
	parameter file:
		ARNO_PARAMS TRUE
	and replace {Ds, Dsmax, Ws, and c} with {d1, d2, d3, and d4}
	in the soil parameter file.

	To disable this behavior, your soil file should contain
	{Ds, Dsmax, Ws, and c} and your global parameter file should
	contain:
		ARNO_PARAMS FALSE
	or you may omit the ARNO_PARAMS line entirely.

Modifications:

	read_vegparam.c:
	    Applied Alan Hamlet's fix for COMPUTE_TREELINE option,
	    which fixed a segmentation fault when COMPUTE_TREELINE=TRUE.
	    This consisted of removing the call to realloc and
	    instead allocating an extra veg class to begin with,
	    as well as assigning this extra veg class a very small
	    fraction of the grid cell's area to avoid changing the
	    results for areas below the treeline.		TJB
	runoff.c:
	    This fixed strange behavior (spikes) in bottom-layer
	    soil moisture and baseflow for low liquid moisture
	    conditions (e.g.  frozen soils).  Changed calculation
	    of dt_baseflow to go to zero when soil liquid moisture
	    <= residual moisture.  Changed block that handles case
	    of total soil moisture < residual moisture to not allow
	    dt_baseflow to go negative.				TJB
	write_model_state.c, read_initial_model_state.c:
	    (bugs found by Lifeng Luo at Princeton)
	    Read_initial_model_state was not parsing ASCII state
	    files according to the format written by
	    write_model_state.  Modified write_model_state.c to
	    add a "\n" after writing mu.  Modified
	    read_initial_model_state to loop over tmp_Nveg and
	    tmp_Nband when searching for the desired cellnum, and
	    to read the correct number of lines per veg class
	    once the desired cellnum has been found.		TJB
	vicNl_def.h, get_global_param.c, read_soilparam.c:
	    Added the ARNO_PARAMS option.  If TRUE, read_soilparam
	    reads the ARNO baseflow parameters d1, d2, d3, and d4
	    and converts them to Ds, Dsmax, Ws, and c.		TJB
	snow_utility.c:
	    In the snow_density function, modified the checks on
	    delta_depth so that the condition is
	    delta_depth > MAX_CHANGE*depth
	    Modified compression due to aging to only be calculated
	    if depth > 0.					TJB
	global.h, vicNl.c, cmd_proc.c:
	    Modified the version information banner to display the
	    version string defined in global.h.  Added -v option
	    to display version information as well.		TJB

September 5, 2003: VIC release 4.0.4beta r3

	snow_utility.c: modification by KAC
	    Added check to keep compression from aging from exceeding
	    the actual depth of the snowpack.
	get_global_param.c, initialize_atmos.c, read_vegparam.c,
	vicNl_def.h: modification by KAC
	    Moved COMPUTE_TREELINE flag from user_def.h to the
	    options structure.  Now when not set to FALSE, the
	    value indicates the default above treeline vegetation
	    if no usable vegetation types are in the grid cell
	    (i.e. everything has a canopy).  A negative value
	    will cause the model to use bare soil.  Make sure that
	    positive index value refer to a non-canopied vegetation
	    type in the vegetation library.

	    To activate treeline calculation, add the following to
	    your global control file:
		COMPUTE_TREELINE n
	    If n < 0, the default above treeline vegetation is
	    set to bare soil.  If n >= 0, the default vegetation type
	    is set to type n in the vegetation library.
	    To deactivate treeline calculation, add the following or remove
	    the line from the global control file:
		COMPUTE_TREELINE FALSE

	    WARNING #1: Since the model does not store default root zone
	    distributions, the default vegetation type will use the
	    values from the last defined vegetation type for the current
	    grid cell (i.e. veg type N-1, prior to the addition of the
	    new vegetation type).

	    WARNING #2: If you are using GLOBAL_LAI, than the LAI and
	    dew storage terms for the default vegetation type will be set
	    to the values used by the last occurrence of the default
	    vegetation type.

	write_model_state.c: bug found by AWW and KAC
	    Modified to print space before dz_node for ASCII state
	    file, this corrects a problem with state files created
	    for models using the Cherkauer and Lettenmaier (1999) heat
	    flux formulation.

June 4, 2003: VIC release 4.0.4beta r2

This covers bugs found during tests with the snow algorithm.

	solve_snow.c: Bug found by KAC
	    Counter for number of days since last snowfall was
	    incremented twice in the MELTING update.  This has been
	    fixed.
	solve_snow.c: modification by KAC
	    Added check so that MELTING flag is only TRUE if melt
	    occurs in the melt season - currently this is defined
	    between March 1 and October 1.  Otherwise the MELTING
	    flag can trigger rapid very early season melt
	write_model_state.c, read_initial_model_state.c, open_state_file.c,
	check_state_file.c: Modification by KAC
	    Modified to handle both ASCII and BINARY state files.
	    NOTE: ASCII mode is designed to make it easier to create
	    a state file for initializing a point model.  It includes
	    all features of the Binary state file, but values are
	    truncated so the resulting restart will not be identical.
	    If you need an exact restart, use the Binary files.  Also
	    removed ice content from the state file as it is computed
	    at the begining of each time step, so storing its value
	    is unnecessary.

April 23, 2003: VIC release 4.0.4beta r1

This covers bug fixes found by beta testers and fixed in the version of
the code bundled with this file.

	surface_fluxes.c: (found by Ingjerd Haddeland)
	    Indexing fix sent SNOW_STEP to calc_surf_energy_bal rather
	    than the model time step, meaning that without snow the
	    evaporation was computed for SNOW_STEP hours rather than a
	    full day.  This was fixed by introducing step_inc to
	    index the arrays, while step_dt keeps track of the correct
	    time step.

March 27, 2003: VIC release 4.0.4beta

This release covers patches for several bugs found in 4.0.3, which were 
never formally released (i.e. the downloadable source code was modified, 
but no major announcement was made).  It also includes other fixes and 
modifications that have been identified as being needed prior to releasing 
version 4.1.0, which will involve several significant changes (including 
lakes & wetlands, spatial snow & frost, and a closed canopy energy balance).

Modifications:

	Snow albedo update: (found by Keith)
	     In previous releases, the snow albedo function has been hyper 
	     sensitive to trace amounts of snowfall during the melt period.  
	     Whenever new snow falls or the cold content falls below 0, the 
	     albedo is switched from the ablation curve to the accumulation 
	     curve.  This curve is then followed until the cold content 
	     exceeds 0, indicating it is in the spring melt season.  This is 
	     fine when accounting for thin early season snowpacks or mid-
	     season melt events, however, a cold snap or light dusting of 
	     snow should not reset the snowpack albedo to much higher winter 
	     values for days or weeks at a time.  This release of the model 
	     monitors the state of pack with the variable MELTING.  This 
	     flag keeps the snowpack on the ablation curve unless new snow-
	     fall exceeds a threshold (TraceSnow in user_def.h) indicating 
	     that the top of the snowpack should be represented by the albedo 
	     of the new snow.
	Model initialization: (found by Ulysses and Keith)
	     In previous releases, the initialization of soil layer moist[] 
	     and ice[] was within a second set of loops for band and 
	     vegetation, using the same counters.  Because of the dual use
	     of counters, initialization was not completed correctly.  The
	     primary effect of this was that thermal node values beyond the
	     first vegetation type were not correctly initialized, which 
	     caused the model to crash during some simulations with full
             energy or frozen soils active.  Without frozen soil, most 
	     simulations would compensate for the problem within their
	     spin-up time, so it is unlikely that this bug impacts any 
	     simulations not employing frozen soil.
	Snow time step: (found by Andy, et al.)
	     The snow algorithm needs to run sub-daily for the energy balance
	     components to function properly.  This means that for daily 
	     simulations, the snow model must be solved at a finer (sub-daily)
	     time step.  In the previous release, initialize_atmos.c stored
	     sub-daily met data in each days variable using positions (e.g. 
             0,1,..8 for 3 hour data).  In surface_fluxes.c the model indexed 
	     the sub-daily time steps used by the snow algorithm with hours 
	     (e.g. 0,3,6,...21 for three hour data).  This means the arrays
	     were incorrectly indexed and the resulting model simulations 
	     would be wrong.  The fix implemented here has been tested under 
	     several model configurations and is deemed the official version.
	     WARNING: There are several versions of this fix circulating,
	     please update your code to this version - the previous fixes
	     may not work in all circumstances!
	FROZEN_SOIL active flag: (found by Ed and Justin)
	     The cause of the problem is a bug in the code that occurs when 
	     the global frozen soils flag (FROZEN_SOILS) is set to true but 
	     the individual cell frozen soil flag (FS_ACTIVE) is set to 
	     false. This causes the change in soil heat storage to be 
	     calculated incorrectly.  This was fixed by adding additional
	     conditions within frozen_soil.c and initialize_model_state.c,
	     which verify the FS_ACTIVE is TRUE before running 
	     estimate_layer_ice_content.  This avoids the problem of the 
	     soil layer ice content being set to a positive value, ignored
	     by the rest of the model.
	Vapor pressure: (Keith)
	     All internal vapor pressure calculations are now done in 
	     Pascals.  Previous release versions, switched between Pa
	     and kPa, so this simply removes the extra step.  The input
	     file format is unchanged, so there should be no change to
	     model output (model might run slightly faster, but it is
	     also unlikely that this will be witnessed by a normal user).
	Constant dew despite changing LAI: (Dave Peterson)
	     The modification of read_vegparam.c to update LAI based on 
	     a grid cell specific value did not change the values of Wdmax.  
	     Wdmax values were computed in read_veglib.c based on the 
	     default LAI values, so they did not necessarily reflect the 
	     actual LAI values used for the grid cell.  Values for Wdmax 
	     are now computed in read_vegparam.c whenever GLOBAL_LAI are 
	     provided.  The effects of this change will change in magnitude
	     based on how different the cell LAI values are from those in
	     the default library file.
	DRY_TIME error: (Reinur Schnur)
	     DRY_TIME in dist_prcp.c was incremented by the time step in 
	     hours.  Then to see if the current rain was part of the same 
	     storm or the start of a new one, DRY_TIME was checked to see 
	     if it was greater than or equal to 24/dt.  This compares 
	     DRY_TIME in hours to the model time step.  The "/dt" has been 
	     removed, so now DRY_TIME is checked versus the hours since the 
	     last storm.
	State file: (KAC)
	     *** WARNING: This may require modifications to your global file ***
	     The state file has been modified to account for model updates.
	     It has also been converted to write binary files - this makes 
	     them less convenient to edit, but means that model starts using
	     the same forcing, soil and vegetation files will produce the
	     exact same results.  There has also been a slight change in how
	     the global file is set up to restart the model.  The global
	     file should now have the same year, month, day and hour as the 
	     original global file - the VIC model will compute the number of
	     records to skip at the beginning to reach the point where the 
	     model state was saved.  This means that caluclations to yield
	     sub-daily meteological forcings from daily forcings will produce 
	     the exact same forcing values -> this also means that restarted
	     simulations will be exactly the same as the original run.  Slight
	     variations in the model results were also introduced because the
	     method for storing soil node depths led to the possibility of 
	     very small differences in dz_node for the restarted model.  
	     Previous versions also did not store the snowpack cold content,
	     this meant that for restarts during winter, the snowpack albedo
	     might start on the ablation curve (as cold content was initialized
	     to 0 and not to a value less than 0) rather than the accumulation
	     curve.  This could lead to slight differences in the snowpack
	     if no new snow fell and the snowpack was not melting - but
	     after 10 years of simulations the differences were minor.  As 
	     noted above the new version of the state file should allow the
	     model to be restarted and to produce exactly the same results as
	     the original complete result.  If there are cases where this is
	     not true, please report.  If you edit the read/write model state
	     functions - BE VERY CAREFUL to edit Nbytes to reflect any changes.
	COMPUTE_TREELINE: (KAC and LCB)
	     This is an added feature which computes the treeline elevation
	     in the current grid cell and does not include vegetation 
	     fractions with overstory in the grid cell averages for snow bands
	     that exceed the treeline elevation.  This feature was added
	     to the model to reduce the appearance of "glaciers" in high
	     elevation snow bands.  It computes average annual July air
	     temperatures using the temperature data from the atmospheric
	     forcing files (WARNING - elevation of treeline may change if
	     the period of forcing data used changes).  It then lapses the
	     average annual July air temperature to locate the elevation 
	     at which it equals 10C.  This is assumed to be the treeline,
	     so vegetation types with overstory in snow bands with average 
	     elevations higher than this, are excluded from the grid cell
	     average variables in put_data.c.  For the time being full 
	     energy and water balances are still computed for these 
	     vegetation fractions, and no attempt is made to verify that
	     a snow band has a non-overstory vegetation type that can be 
	     expanded to represent the coverage area lost due to the 
	     exclusion of the overstory fraction.
	     

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
is returned without new data values.

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
