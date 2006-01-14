# README.txt - Release Notes
#
# $Id$
#------------------------------------------------------------------------

--------------------------------------------------------------------------------
                     ***** VIC MODEL RELEASE NOTES *****
--------------------------------------------------------------------------------

Usage:
------

	vicNl [-v | -o | -g<global_parameter_file>]

	  v: display version information
	  o: display compile-time options settings (set in user_def.h)
	  g: read model parameters from <global_parameter_file>.
	     <global_parameter_file> is a file that contains all needed model
	     parameters as well as model option flags, and the names and
	     locations of all other files.


--------------------------------------------------------------------------------
***** Description of changes from VIC 4.1.0 beta r2 to VIC 4.1.0 beta r3 *****
--------------------------------------------------------------------------------


New Features:
-------------

Improved lake model

	Files Affected:
	CalcBlowingSnow.c, LAKE.h, IceEnergyBalance.c, ice_melt.c,
	initialize_atmos.c, initialize_lake.c, lakes.eb.c, read_lakeparam.c,
	solve_snow.c, surface_fluxes.c, vicNl.h, vicNl_def.h,
	water_energy_balance.c, water_under_ice.c

	Description:
	This fixes the following bugs:

		Lake model crashes when lakes fill up
		Lake model crashes when lakes dry out

	And introduces the following new features:

		Blowing snow sublimation computed over lakes
		Tuning of default lake profile (parabolic to square)

	The new parameter BETA in LAKE.h can be used to control the shape of the
	automatic lake profile (when the global option LAKE_PROFILE is FALSE).
	BETA=0.5 is strongly parabolic, while BETA=0.001 is almost square
	(vertical walls).

	If LAKE_PROFILE=FALSE, the lake parameter file from VIC 4.1.0 beta 2 can
	be used without any changes.  However, if LAKE_PROFILE=TRUE, the lake
	parameter file must now contain a (depth,area) pair for each lake node.
	See read_lakeparam.c for details.
	

VIC now supports ALMA input and output variables.

	Files affected:
	Makefile, conv_force_vic2alma.c, conv_results_vic2alma.c, close_files.c,
	display_current_settings.c, dist_prec.c, get_force_type.c,
	get_global_param.c, initialize_atmos.c, initialize_global.c,
	make_in_and_outfiles.c, put_data.c, read_atmos_data.c, vicNl.h,
	vicNl_def.h, write_data.c, write_forcing_file.c

	Description:
	VIC now supports ALMA input and output variables.  To have VIC read ALMA
	input (forcing) variables, all that is required is to identify them by
	the appropriate names in the forcing section of the global parameter file.
	The appropriate names are listed in vicNl_def.h.  To have VIC write ALMA
	output variables, it is necessary to include the option ALMA_OUTPUT in the
	global parameter file and set it to TRUE.  This option does not affect
	model physics.


State file is now written at the END of the final timestep of the date indicated
in the global parameter file.

	Files affected:
	dist_prec.c

	Description:
	Previous versions of VIC wrote the state file at the beginnning of the
	date indicated (via STATEYEAR, STATEMONTH, and STATEDAY) in the global
	parameter file.  For sub-daily model time steps, this meant that the
	simulation's end could not be aligned exactly with the writing of the
	state file, since the model would continue simluating for the remainder
	of the final day after writing the state file.  This created difficulties
	when running VIC in real-time simulations, in which the forcings are
	broken up into short time periods and the previous state file must be
	used as the initial condition for the next time period.

	Now that the state file is written at the end of the final timestep of the
	date indicated in the global parameter file, it is possible for the end of
	a forcing file, the end of the simulation, and the timing of the state file
	all to align on the same time.  This simplifies re-starting the model for
	a forcing file that begins immediately after the previous forcing file
	ended, since the state file now is equivalent to the initial condition at
	the beginning of the new forcing file.


EQUAL_AREA global parameter option

	Files affected:
	display_current_settings.c, get_global_param.c, initialize_global.c,
	read_lakeparam.c, vicNl_def.h

	Description:
	New global parameter option.  When EQUAL_AREA is TRUE, all grid cells are
	assumed to have the same area, and the global parameter RESOLUTION is
	interpreted to store the grid cell area in km^2.  When EQUAL_AREA is FALSE,
	grid cells are assumed to have the same side length in degrees, and
	RESOLUTION is interpreted to store the grid cell side length in degrees,
	as before.


Improved validation of global options

	Files affected:
	get_global_param.c

	Description:
	Now does not allow bad combinations of FULL_ENERGY, FROZEN_SOIL, and
	GRND_FLUX values.


Improved validation of soil param file

	Files affected:
	read_soilparam.c

	Description:
	Now validates depth_full_snow_cover and frost_slope.


New vicInterp executable

	Files affected:
	Makefile, check_files.c, close_files.c, get_global_param.c,
	read_soilparam.c, vicNl.c

	Description:
	Now users can create a stand-alone executable called "vicInterp" by
	typing "make interp".  This executable is a stripped-down version 
	of VIC that reads a set of forcing files and outputs diurnally-varying
	sub-daily forcings.  The interpolation scheme, based on the Thornton and
	Running (mtclim) algorithm, can convert from from daily Prec, Tmin, Tmax,
	and Wind to sub-daily Prec, AirTemp, Shortwave and Longwave radiation,
	Wind, Pressure, and air Density.  The interpolated forms of Prec, Wind,
	Pressure, and Density do not vary diurnally but are assumed constant
	during each day.

	vicInterp is vicNl, compiled with the OUTPUT_FORCE option set to TRUE.
	As such, vicInterp can read a standard VIC global parameter file.  While
	it will ignore many of the options, it understands (and requires) the
	following:
		TIMESTEP
		STARTYEAR
		STARTMONTH
		STARTDAY
		STARTHOUR
		ENDYEAR
		ENDMONTH
		ENDDAY
		GRID_DECIMAL
		FORCING1
		N_TYPES
		FORCE_TYPE
		FORCE_FORMAT
		FORCE_ENDIAN
		FORCE_DT
		FORCEYEAR
		FORCEMONTH
		FORCEDAY
		FORCEHOUR
		SOIL
		RESULT_DIR
	In addition, the following options are optional:
		BINARY_OUTPUT
		ALMA_OUTPUT

	For example, if FORCE_DT is 24 and TIMESTEP is 3, vicInterp will
	interpolate the input daily forcings to a 3-hour time step.  If
	ALMA_OUTPUT is TRUE, these output forcings will be standard ALMA
	variables.

	This new feature also involves an update to the OUTPUT_FORCE option
	that fixes a bug in which output files were not properly closed
	before exiting.  This prevented all of the output data from being
	completely written.

Added sample global param file to distribution

        Files affected:
        global.param.sample

Bug Fixes:
----------

Large water balance errors in daily water balance mode when snow is present

	Files affected:
	surface_fluxes.c

	Description:
	Fixed broken snow step in surface_fluxes.c.  Per-snow-step snow
	quantities were being reset at the beginning of each canopy-surface
	energy balance iteration, preventing snow characteristics from being
	accumulated over all snow steps, and resulting in large water balance
	errors.  Now per-iteration snow quantities are stored separately
	from per-snow-step quantities.


Aerodynamic resistance incorrect in output fluxes file

        Files affected:
        IceEnergyBalance.c, LAKE.h, SnowPackEnergyBalance.c,
        calc_surf_energy_bal.c, full_energy.c, func_canopy_energy_bal.c,
        func_surf_energy_bal.c, ice_melt.c, lakes.eb.c, put_data.c,
        snow_intercept.c, snow_melt.c, solve_snow.c, surface_fluxes.c,
        vicNl.h, vicNl_def.h, wetland_energy.c

        Description:
        In 4.1.0 beta r2 and earlier, the aerodynamic resistance written to the
        output fluxes file rarely reflected the value actually used in flux
        computations.  This has been fixed.

        VIC actually computes an array of 3 different aerodynamic resistances,
        as follows:
          aero_resist[0] : over vegetation or bare soil
          aero_resist[1] : over snow-filled overstory
          aero_resist[2] : over snow pack
        VIC determines which element of the array to use depending on the
	current vegetation type and whether snow is present.  In addition, in
	most cases, VIC applies a stability correction to this aerodynamic
	resistance before using it in flux computations.  Furthermore, when the
	current vegetation tile contains overstory and snow is present on the
	ground, aero_resist[2] is used for snow pack flux computations and
	either aero_resist[1] or aero_resist[0] is used for canopy flux
	computations, meaning that two different aerodynamic resistances are in
	use in the same time step.

        However, VIC 4.1.0 beta r2 always wrote the uncorrected value of
        aero_resist[0] to the fluxes file.

        In 4.1.0 beta r3 and later, the value written to the fluxes file is the
        actual value used in flux computations, including any corrections that
        were applied.  In the case mentioned above in which two different
        aerodynamic resistances are in use at the same time, the one used for
        the snow pack is written.  In addition, the aerodynamic resistance of
        the canopy air/surrounding atmosphere interface is not tracked.


Aerodynamic resistance not correctly aggregated for output

	Files affected:
	conv_results_vic2alma.c, put_data.c, vicNl_def.h

	Description:
	In previous releases, aerodynamic resistance (out_data->aero_resist)
	was aggregated by a simple area-weighted average over veg tiles.  This
	led to an aggregate value that was not the true effective resistance
	of the entire grid cell.  Since evaporation is proportional to
	1/aero_resist, it is (1/aero_resist), or the aerodynamic conductivity,
	that should be averaged over the grid cell.  Therefore, a new variable,
	out_data.aero_cond, was created for the purposes of aggregation.  After
	aggregation, out_data.aero_resist is computed as 1/out_data.aero_cond.

	The effect of the change is most pronounced when there is a large range
	of values of aerodynamic resistance among the veg tiles in a grid cell.
	The effective aerodynamic resistance, computed the new way, will tend
	to be smaller than the old way when the values cover a wide range.
	However, the effective aerodynamic resistance will never be smaller than
	the smallest value among the various veg tiles in the cell.
	

Attempts to skip deactivated grid cells fail when using a binary initial state
file.

	Files affected:
	write_model_state.c

	Description:
	write_model_state() now computes the correct number of bytes per record
	for a binary state file.  Before this fix, if a user de-activated some
	grid cells (by setting the first field in the soil parameter file to 0)
	and attempted to read from a binary state file that was written when the
	cells were active, VIC would not read the state file correctly and would
	crash.


State variables for SPATIAL_FROST and LAKE_MODEL options not stored in state
file.

	Files affected:
	dist_prec.c, initialize_model_state.c, read_initial_model_state.c,
	vicNl.h, write_model_state.c

	Description:
	State variables for SPATIAL_FROST and LAKE_MODEL options were not stored
	in state files, causing these variables to be reset to default values
	when reading from an initial state file.  This has been fixed.

	NOTE: this fix involves a change in the format of state files which
	makes state files used with 4.1.0 beta 3 incompatible with earlier
	releases (and vice versa).


State files not written/read correctly when QUICK_FLUX set to TRUE.

	Files affected:
	get_global_param.c, make_dist_prcp.c, make_energy_bal.c, vicNl.c,
	vicNl.h

	Description:
	In previous releases, if QUICK_FLUX=TRUE, the number of soil thermal
	nodes is changed after being recorded in the state file.  The resulting
	mismatch between Nnodes in the state file header and the actual number
	of node temperatures recorded in the state file prevents VIC from being
	able to read the state file.  This has been fixed.


Automatic re-adjustment of lake fraction when sum of veg and lake fractions
exceeds 1 can set lake fraction to 0.

	Files affected:
	LAKE.h, read_lakeparam.c, vicNl.c

	Description:
	In earlier releases, the automatic re-adjustment of the lake fraction
	when the sum of the veg and lake fractions exceeds 1 can reset the lake
	fraction to 0.  This has been changed so that the veg and lake
	fractions all share proportionally in the re-adjustment.


Fixed incorrect check on soil moisture in distribute_node_moisture_properties().

	Files affected:
	soil_conduction.c

	Description:
	distribute_node_moisture_properties() contained the following check:
		if(abs(moist_node[nidx]-max_moist_node[nidx]) > SMALL)
	Checking the absolute value here was incorrect.  The abs() has been
	removed.


Fixed incorrect check on soil node depths in read_initial_model_state().

	Files affected:
	read_initial_model_state.c

	Description:
	read_initial_model_state() contained the following check:
		if( abs( sum - soil_con->dp ) > SMALL )
	Checking the absolute value here was incorrect.  The abs() has been
	removed.


ARNO_PARAMS global parameter option changed to NIJSSEN2001_BASEFLOW

	Files affected:
	display_current_settings.c, get_global_param.c, initialize_global.c,
	read_soilparam.c, read_soilparam_arc.c, vicNl_def.h

	Description:
	Changed the name of the ARNO_PARAMS global parameter option to
	NIJSSEN2001_BASEFLOW.  The meaning of the ARNO_PARAMS option was
	actually opposite to its name: when ARNO_PARAMS was FALSE, VIC would
	interpret the first four parameters in the soil parameter file to be the
	standard ARNO soil parameters Ds, Dsmax, Ws, and c, while when ARNO_PARAMS
	was TRUE, VIC would interpret the first four parameters to be d1, d2, d3,
	and d4, the soil parameters used in Nijssen et al. (2001).  The new name
	for this option more accurately reflects its meaning: when
	NIJSSEN2001_BASEFLOW is TRUE, VIC assumes the soil parameter file contains
	d1, d2, d3, and d4.  When NIJSSEN2001_BASEFLOW is FALSE, VIC assumes the
	soil parameter file contains Ds, Dsmax, Ws, and c.

	As of the current release of VIC 4.1.0, VIC accepts both ARNO_PARAMS and
	NIJSSEN2001_BASEFLOW in the global parameter file.  But eventually
	ARNO_PARAMS will be phased out, and users are encouraged to replace
	ARNO_PARAMS with NIJSSEN2001_BASEFLOW in their global parameter files.


Replaced %i with %d in fscanf statements

	Files affected:
	check_state_file.c, get_global_param.c, read_arcinfo_ascii.c,
	read_initial_model_state.c, read_snowband.c, read_soilparam.c,
	read_veglib.c

	Description:
	Having %i in fscanf statements was causing input values of "08" to be
	interpreted as octal rather than decimal.  These instances of %i have
	been replaced with %d.

Removed some snowband output

	Files affected:
	write_data.c

	Description:
        Removed  net sw radiation, net lw, albedo, latent heat flux, 
	sensible heat flux, ground heat flux

STATE file option is now specified in global file, not in user_def.h at compile time

        Files affected:

        display_current_settings.c
        dist_prec.c
        get_global_param.c
        initialize_global.c
        open_state_file.c
        user_def.h
        vicNl.c
        vicNl.h
        vicNl_def.h
        write_model_state.c

        Previously, to be able to read/write state files, VIC required
        users to define SAVE_STATE to be TRUE in user_def.h and
        recompile VIC, in addition to having the correct settings of
        INIT_STATE, STATENAME, STATEYEAR, STATEMONTH, and STATEDAY in
        the global parameter file.

        This was both unnecessary and confusing.  Setting SAVE_STATE
        to FALSE did not improve VIC performance noticeably, and
        having both a SAVE_STATE compile option and state file
        information in the global parameter file only led to mistakes
        in which the user made a change in one place without realizing
        a change needed to be made elsewhere.

        Now VIC behaves as follows:
        1. To read an initial state file, the following line
        must be in the global parameter file:
          INIT_STATE init_state_filename

        where init_state_filename is the name of the initial state
        file.  If this line is absent or commented out (preceded by a
        "#"), VIC will start from default initial conditions.

        2. To write a state file during the simulation, the
        following lines must be in the global parameter file:
          STATENAME  state_file_name
          STATEYEAR  state_year
          STATEMONTH state_month
          STATEDAY   state_day

        where state_file_name is the path and prefix of the file to
        save the state information in (the simulation date on which
        the state file is saved will be appended to the state file
        name), and state_year, state_month, and state_day are the
        year, month, and day of the simulation on which to save state.
        The state will be saved AFTER the FINAL time step of that
        date.  If all of these lines are absent or commented out, VIC
        will not save a state file.  If some (but not all) of these
        lines are present, VIC will give an error.

Fixed bug in error trapping when INIT_STATE filename matches SAVE_STATE filename

        Files affected:
        get_global_param.c, open_state_file.c

        Previously, the checks for matches would occur in get_global_param even
        though the SAVE_STATE filename wasn't created until open_state_file.
        The output name setting was moved from open_state_file to get_global_param. GCT

Added checks for range/valid month days

        Files affected:
        get_global_param.c

        In previous versions the user could set a non-valid date for STATE files
        and the model would run without writing to STATE file. Code now checks
        for valid date.

Allow user to use NO_FLUX in addition to NOFLUX for NOFLUX option in global param file

        Files affected:
        get_global_param.c

Now reads extra_veg from state file

        Files affected:
        read_initial_model_state.c
               
Removed Trad from vicNl_def.h

        Files affected:
        vicNl_def.h

Modified to read lake nodal variables for each of the active nodes. (JCA)

        Files affected:
        read_initial_model_state.c
--------------------------------------------------------------------------------
***** Description of changes from VIC 4.1.0 beta r1 to VIC 4.1.0 beta r2 *****
--------------------------------------------------------------------------------


New Features:
-------------

"-v" and "-o" command-line options and display of run/compile-time options

        Files affected:
        cmd_proc.c, display_current_settings.c, get_global_param.c, global.h,
        vicNl.c, vicNl.h, vicNl_def.h

        Description:
        In VIC 4.1.0 beta r2, if VERBOSE is TRUE, all compile-time options (from
        user_def.h) and run-time options (from your global parameter file) are
        displayed for you at the beginning of a model run.  If you are saving
        your model output to a log file, this information will be stored with
        your log, so you will have a record of the option settings that produced
        the results of that particular model run.

	The new "-v" option will display the release number of your vicNl
	executable.  For example, typing:
		vicNl -v
	gives
		***** VIC Version 4.1.0 Beta Release 2 *****

        Meanwhile, the new "-o" option displays the current compile-time
        options.  One benefit of this option is that you can see what the
        options in user_def.h were set to when vicNl was compiled, without
        having to start a model run.  Since your version of user_def.h may have
        changed since you compiled vicNl, this is the most reliable way to see
        what these options are set to.  For example:

	gen.hydro.washington.edu 175: vicNl -o

	***** VIC Version 4.1.0 Beta Release 2 - Current Model Settings *****

	COMPILE-TIME OPTIONS (set in user_def.h)
	----------------------------------------

	Output to Screen:
	OUTPUT_FORCE_STATS      FALSE
	VERBOSE                 TRUE

	Input Files:
	NO_REWIND               TRUE

	Output Files:
	LDAS_OUTPUT             FALSE
	LINK_DEBUG              FALSE
	OPTIMIZE                FALSE
	OUTPUT_FORCE            FALSE
	SAVE_STATE              FALSE

	Simulation Parameters:
	CLOSE_ENERGY            FALSE
	COMPUTE_TREELINE        FALSE
	LAKE_MODEL              FALSE
	LOW_RES_MOIST           FALSE
	QUICK_FS                FALSE
	SPATIAL_FROST           FALSE
	SPATIAL_SNOW            FALSE

	Maximum Array Sizes:
	MAX_BANDS               10
	MAX_FRONTS               3
	MAX_LAKE_NODES          20
	MAX_LAYERS               3
	MAX_NODES               18
	MAX_VEG                 12

	Snow Constants:
	NEW_SNOW_ALB            0.850000
	SNOW_ALB_ACCUM_A        0.940000
	SNOW_ALB_ACCUM_B        0.580000
	SNOW_ALB_THAW_A         0.820000
	SNOW_ALB_THAW_B         0.460000
	TraceSnow               0.030000

	Other Constants:
	LAI_WATER_FACTOR        0.200000
	LWAVE_COR               1.000000
	MAXIT_FE                25


Automatic recompilation on updates to *.h

        Files affected:
        Makefile

        Description:
        In VIC 4.1.0 beta r1, updating a .h file and recompiling VIC would not
	result in recompilation of files that depend on the .h file, unless a
	"make clean" command was issued first.  Now, if any .h files are
        updated, all dependent .c files are recompiled on the next "make".


NEW_ARNO_TYPE global option is now ARNO_PARAMS

	Files affected:
	get_global_param.c, read_soilparam.c, vicNl_def.h

	Description:
	Changed the name of the global option NEW_ARNO_TYPE to be ARNO_PARAMS.
	"NEW_ARNO_TYPE" is confusing, since this option will not be new
	forever, and doesn't refer to a "type" but rather a set of parameters.

	NOTE: This change requires the user to replace all occurrences of
	"NEW_ARNO_TYPE" in their global control files with "ARNO_PARAMS".


Bug Fixes:
----------

Spurious condensation at low temperatures

        Files affected:
        arno_evap.c

        Description:
        Changed logic of evap limit check to avoid creating spurious
        condensation.  In VIC 4.1.0 beta r1, whenever evaporation > (liquid
        moisture - residual moisture), evaporation would be set to (liquid
        moisture - residual moisture).  However, at low temperatures, when most
        or all soil moisture is frozen, liquid moisture < residual moisture,
        causing (liquid moisture - residual moisture) to be negative.  Any non-
        negative evap would be greater than this, resulting in evap getting set
        to (liquid moisture - residual moisture), which would be negative (i.e.
        condensation).  This artificially created condensation in whatever
        amount necessary to bring liquid moisture up to residual, causing 1)
        large latent heat flux, 2) incorrect surface temperatures, 3) occasional
        inability for calc_surf_energy_bal to converge in root_brent, and 4)
        spuriously high runoff and baseflow.  Now there is an added condition
        that liquid moisture > residual moisture for evap to be capped at
        (liquid moisture - residual moisture).

        NOTE: This fix results in lower runoff and baseflow in unvegetated areas
        with frozen soils, and may require recalibration of soil parameters.


Validation of initial soil moisture

        Files affected:
        read_soilparam.c, read_soilparam_arc.c

        Description:
        Added check to read_soilparam.c and read_soilparam_arc.c to make sure
        that wilting point is greater than residual moisture.  Changed lower
        limit on initial soil moisture to be residual moisture rather than
        wilting point.  Made validation statements clearer.  Validation does
        not occur if soil moisture will be read from a model state file.
        (found and fixed by Chunmei Zhu, Alan Hamlet, and Ted Bohn)

        NOTE: Soil parameter files containing Wpwp_FRACT and resid_moist
        such that Wpwp_FRACT < resid_moist / (1.0 - bulk_density/soil_density)
        will now cause VIC to exit with an error message.


Incorrect baseflow limits

        Files affected:
        runoff.c

        Description:
        In 4.1.0 beta r1, runoff.c checked for the wrong bounds on baseflow,
        allowing baseflow to become negative when liquid soil moisture < residual
        moisture.  These bounds have been fixed in 4.1.0 beta r2, as follows:
        baseflow is not allowed to exceed (liquid soil moisture - residual moisture);
        when baseflow < 0, baseflow is set to 0; when baseflow > 0 and the resulting
        soil moisture < residual moisture, water is taken out of baseflow and
        given to the soil as follows:
          if baseflow > (residual moisture - soil moisture), then
            baseflow -= (residual moisture - soil moisture);
            soil moisture += (residual moisture - soil moisture);
          else
            soil moisture += baseflow;
            baseflow = 0;

        NOTE: This fix may result in small changes in baseflow and evaporation.


Runs using initial state files starting at state file date rather than global
start date

        Files affected:
        check_state_file.c

        Description:
        In 4.1.0 beta r1, check_state_file.c would increment the index of the
        forcing data array until it reached the record corresponding to the date
        stored in the state file.  This caused the simulation to start at the
        date at which the state file was saved rather than the start date
        specified in the global parameter file.  If the state file's date was
        earlier than the start date in the global parameter file, the index
        would be incremented until a segmentation fault occurred.  This has
        been fixed in 4.1.0 beta r2 so that the start date in the global parameter
        file is always the start date of the simulation.  The date stored in
        the initial state file is ignored.

        NOTE: If you have been relying on the state file to dictate when your
        simulations start, this fix may require you to change your global
        parameter file so that STARTYEAR, STARTMONTH, etc. reflect the start
        date/time you want.


Incorrect sublimation values for BLOWING option

	Files affected:
	CalcBlowingSnow.c, IceEnergyBalance.c, SnowPackEnergyBalance.c,
	calc_surf_energy_bal.c, func_surf_energy_bal.c, ice_melt.c, lakes.eb.c,
	latent_heat_from_snow.c, put_data.c, snow_melt.c, solve_snow.c,
	surface_fluxes.c, vicNl.h, vicNl_def.h

	Description:
	Fixed 3 bugs in the sublimation terms: 1) sub_surface was wrong when
	snow step was not 1 hour, 2) sub_blowing was wrong under certain
	conditions, and 3) sub_blowing and sub_surface did not contain the
	contribution from lakes, even though sub_total did.  The fix
	establishes the convention that the internal variables VaporMassFlux,
	BlowingMassFlux, and SurfaceMassFlux always have units of kg/m2s; and
	that internal variables vapor_flux, blowing_flux, and surface_flux
	always have units of m/timestep.  Unnecessary terms were removed from
	the parameter lists of several functions.

	NOTE: The effects of this fix on major water balance terms such as
	runoff, baseflow, evaporation, etc. should be very small.

	NOTE 2: The lake contribution to sub_blowing is currently set to 0.


Negative incoming shortwave radiation at high latitudes

        Files affected:
        mtclim42_vic.c

        Description:
        In 4.1.0 beta r1, when sub-daily shortwave radiation is estimated from
        daily min/max temperatures, negative values occasionally are calculated
        in coastal areas above the Arctic Circle in winter.  Now, if estimated
        sub-daily incident shortwave is negative, it is set to 0.0.


Undefined daily precipitation for deserts

        Files affected:
        mtclim42_vic.c

        Description:
        In 4.1.0 beta r1, if a grid cell's annual precipitation (specified in
        the soil parameter file) is 0, then the adjusted daily precipitation
        calculated in mtclim42_vic.c ends up being undefined.  In 4.1.0 beta r2
        this has been fixed.

        More specifically, the MTCLIM 4.2 algorithm (which VIC uses for
        estimating sub-daily forcing values) was originally set up to expect a
        base precipitation as an input, and to translate this base precip into
        a site-specific precip via an adjustment function.  The adjustment
        function multiplies the base precip by the ratio of site_isohyet to
        base_isohyet to get site precipitation.  In calling the MTCLIM 4.2
        functions, VIC sets site_isohyet and base_isohyet to the grid cell's
        annual precipitation.  So this ratio should always = 1.  However,
        when annual precipitation is 0.0, this ratio is undefined.  So, the
        fix is to set the ratio to 1 when both site_isohyet and base_isohyet
        are 0 (or very small). (found by Liz Clark)


Snow_flux incorrectly set to Tcanopy in fluxes output file

	Files affected:
	put_data.c

	Description:
	In 4.1.0 beta r1, the snow_flux output variable was incorrectly set to
	Tcanopy.  This has been corrected.


Incorrect value for sub_snow in fluxes output file

	Files affected:
	write_data.c

	Description:
	Replaced output of sub_snow[0] in fluxes file with sub_total.


Mixmax uninitialized

	Files affected:
	water_under_ice.c

	Description:
	The variable mixmax was used in a max-finding loop without first being
	initialized.  It is now initialized to 0 before being used.


Typo in SPATIAL_FROST code in snow_intercept.c

	Files affected:
	snow_intercept.c

	Description:
	Fixed typo.  Changed SPATIAL_FRoST to SPATIAL_FROST.


Incorrect skipping over masked-out cells when reading initial state file

	Files affected:
	read_initial_model_state.c

	Description:
	In 4.1.0 beta r1, the algorithm parsing the initial state file was
	incorrect when attempting to skip over records for cells that were
	excluded from the mask (by setting the "active" flag in the first field
	of the soil file to 0).  This has been fixed.

	More specifically, the algorithm has been modified to loop over
	tmp_Nveg and tmp_Nband when searching for desired cellnum in ASCII file,
	rather than over Nveg and Nbands.  As we skip over other records in the
	state file while searching for the desired record, the loop must parse
	each undesired record differently, according to how many veg classes and
	snow bands exist in the record (tmp_Nveg and tmp_Nband, respectively),
	rather than the number of veg classes and snow bands in the desired
	record (Nveg and Nbands, respectively).


Incorrect bounds for snow depth

	Files affected:
	snow_utility.c

	Description:
	Modified the checks on delta_depth so that the condition is
		delta_depth > MAX_CHANGE*depth
	Modified compression due to aging to only be calculated if depth > 0.


100% snow when air_temp = MAX_SNOW_TEMP

        Files affected:
        calc_rainonly.c

        Description:
        In 4.1.0 beta r1, when air_temp = MAX_SNOW_TEMP, the portion of
	precipitation that is snow was set to 100%.  This has been fixed.
	(found by Justin Sheffield at Princeton)


Special case in Penman equation

        Files affected:
        penman.c

        Description:
        Changed
                if (vpd > 0.0 && evap < 0.0)
        to
                if (vpd >= 0.0 && evap < 0.0)
        to correctly handle evap when vpd == 0.0. (found by Justin Sheffield
        at Princeton)


Rint() function not supported on all platforms

        Files affected:
        compute_dz.c, initialize_atmos.c, read_soilparam.c,
        read_soilparam_arc.c

        Description:
        Replaced rint(something) with (float)(int)(something + 0.5) to
        handle rounding without resorting to rint(), which isn't supported
        on all platforms.  (found by Justin Sheffield at Princeton)


Global parameter initialization

        Files affected:
        initialize_global.c

        Description:
        Initialize ARC_SOIL, COMPRESS, and ARNO_PARAMS to FALSE.  Also changed
        limit on loop over forcing types from hard-coded 17 to variable
        N_FORCING_TYPES.  (found by Justin Sheffield at Princeton)


Bottom soil node thickness initialization

        Files affected:
        initialize_model_state.c

        Description:
        Initialize soil_con->dz_node[Nnodes] to 0.0, since it is accessed in
        set_node_parameters().  (found by Justin Sheffield at Princeton)


Initialization of debug parameters

        Files affected:
        open_debug.c

        Description:
        Initialize debug_store_moist array when debug.PRT_MOIST is true (as
        well as under the other previously-defined conditions). (found by
        Justin Sheffield at Princeton)


Canopy evaporation and distributed precipitation

        Files affected:
        surface_fluxes.c

        Description:
        Fixed initialization of canopyevap to initialize for every value of
        dist, rather than just dist 0. (found by Justin Sheffield at Princeton)


Output debug file error

        Files affected:
        write_atmosdata.c

        Description:
        No longer close the debug file, since the next cell must write to it.
        (found by Justin Sheffield at Princeton)


Calculation of deltaH when FS_ACTIVE is FALSE

        Files affected:
        func_surf_energy_bal.c

        Description:
        Added check that both FS_ACTIVE and FROZEN_SOIL are true before
        adjusting *deltaH.  This is just a safety measure; ice and ice0 should
        both be 0 when FS_ACTIVE is FALSE.  (found by Justin Sheffield at
        Princeton)


Root_brent error message clarification

        Files affected:
        calc_atmos_energy_bal.c, calc_surf_energy_bal.c, frozen_soil.c,
	ice_melt.c, root_brent.c, snow_intercept.c, snow_melt.c, vicNl.h

        Description:
        Instead of printing warning or error messages, root_brent.c now
	passes descriptions of the errors to the functions that called it,
        leaving it to them to describe the specific error and its consequences.
	In 4.0.4, root_brent's messages sometimes were wrong.


Display current grid cell number for arc/info soil files

        Files affected:
        read_soilparam_arc.c

        Description:
        In 4.1.0 beta r1, for arc/info-format soil files, the current grid
        cell is not displayed (for regular-format soil files, it is displayed).
        Now the current grid cell number is always displayed.


Binary file opening message

        Files affected:
        open_file.c

        Description:
        In 4.1.0 beta r1, the message announcing the opening of a binary input
        file for reading was truncated.  This has been fixed.


Extra options in user_def.h

	Files affected:
	user_def.h

	Description:
	Removed SPATIAL_FROST_SLOPE and MAX_FULL_COVERAGE_DEPTH since they were
	not used anywhere.


Unnecessary freeing of LINK_DEBUG

	Files affected:
	alloc_atmos.c

	Description:
	Added check for LINK_DEBUG global option.  If LINK_DEBUG is TRUE
	atmospheric data is not dynamically allocated, so it should not be
	freed.


Lakeparam element of filenames structure overloaded

	Files affected:
	close_files.c, make_in_and_outfiles.c, vicNl_def.h

	Description:
	The lakeparam element of the filenames structure was used to contain
	both the name of the input lakeparam file and the output lake file.
	An extra element has been added to the filenames structure so that the
	input lakeparam file and the output lake file can be referred to
	separately.


Inconsistent format in state file

	Files affected:
	write_model_state.c

	Description:
	Removed initial space on veg/band info line in ASCII file.


Pending Issues:
---------------

Lake model crashes when lake level gets too high

Lake model crashes when lake level gets too low

COMPUTE_TREELINE option non-functional

	Description:
	This issue is waiting for completion of a treeline-computation scheme
	that allows landcover fractions to be stated explicitly as a function
	of elevation band.

Lakes smeared across multiple elevation bands

	Description:
	Currently, VIC assumes that the percentage of grid cell area covered by
	lakes is constant throughout the grid cell, and therefore when multiple
	elevation bands are specified for a grid cell, the lake percentage is
	distributed evenly across all elevations.  However, when modeling large
	lakes, this assumption does not hold; the bulk of the lake coverage may
	be due to a single large lake, which by definition cannot exist in
	multiple elevation bands.  Therefore, we need a way of explicitly
	specifying lake coverage as a function of elevation.


--------------------------------------------------------------------------------

Updates 7-30-2003: VIC 4.1.0 r1

    (1) Added support for ASCII as well as Binary state files (NOTE: ASCII
	state files were added for easier editing when trying to start the 
	model from observed state - Only Binary state files will produce
	identical runs if used to store and restart from the model state).
	Added BINARY_STATE_FILE to the global file, if TRUE state file I/O
	is in binary, if FALSE it is ASCII.
    (2) Made certain that local NOFLUX flag is set in calc_surf_energy_bal.c
	every time the routine is executed.
    (3) write_data.c: Corrected output of sub_snow variable to item [0] 
	rather than a point - will need to decide what parts of this array 
	are important to output.
    (4) Modified runoff.c so that only the top two layers are used when
	computing infiltration.  Previously, all but the bottom layer was 
	used.  Therefore this only affects you if you tried running with
	more than three layers.  This change was made as it makes more sense 
	for a situation where you are trying to improve the representation
        of the soil column by increasing the resolution (and perhaps depth)
	of the soil column.
    (5)	snow_utility.c (snow_density): Added check to keep compression from 
	aging from exceeding the actual depth of the snowpack.
    (6) solve_snow.c: Added check so that MELTING flag is only TRUE if melt
        occurs in the melt season - currently this is defined between March 1 
	and October 1.  Otherwise the MELTING flag can trigger rapid very 
	early season melt.

Updates 4-21-2003: Mering codes to form VIC 4.1.0 r0

-----------------------------------------------------------------------------

    From comparison of LatestAdminSource and LatestTestSource

    Dynamic allocation of as many arrays and structures as are currently 
    used by the VIC model cause electric fence (a memory debugging tool)
    to bog down/fail in the allocation process.  As the allocation of
    memory is unlikely to be the cause of new memory errors, hard wired 
    memory can be used instead, letting eleectric fence penetrate further
    into the code.
    alloc_atmos.c: (KAC) MEMORY_DEBUG update
	Modified so that dynamic allocation of the atmospheric data arrays
	does not occur when the debugging code is linked to the model.  
	This makes it possible to use electric fence to debug memory errors
	without changing the model code.  Previously, dynamic allocation
	of the atmospheric arrays was only skipped if VIC was configured 
	for optimization work.
    vicNl_def.h: (KAC) MEMORY_DEBUG update
        See alloc_atmos.c.

    Newest version of the model had water balance errors when using 
    distributed precipitation.  Model crashes were also witnessed.  Fixed
    errors introduced by other model updates.
    arno_evap.c: (KAC) DIST_PRCP fix
	Moved unit conversion of moist_resid outside of the distributed
	precipitation loop.  This prevents it from being multiplied by 
	D1 * 1000 twice for the DRY fraction of the grid cell.
    initialize_model_state.c: (KAC) DIST_PRCP fix
        Modified to initialize soil and vegetation parameters for the dry 
	grid cell fraction, if distributed precipitation is activated.
    snow_intercept.c: (KAC) DIST_PRCP fix
        Added check of intercepted before mass balance error calculation.  
	Since interception quantity can be greater than Wd_max when snow 
	is present, the routine could return dew values higher than maximum 
	in a time step when intercepted snow melted.  If there is no other 
	snow, and distributed precipitation is active this could cause the 
	model to crash in redistribute_during_storm as it will be unable to 
	conserve water when too much water is in the canopy.

    ** Make sure that all updates from surface_fluxes.c are accounted for
    after merge with 4.0.4beta **

-----------------------------------------------------------------------------

    From comparison of LatestAdminSource (updated above) and LaurasNewestSource

    IceEnergyBalance.c: (LCB) LAKE_WETLAND updates
        SnowFlux variable renamed to qf to match the lake-ice model 
	documentation.
    LAKE.h: (LCB) LAKE_WETLAND updates
        Modified some model parameters (need to check with Laura as to their
	applicability outside dissertation basins).  Also updated profiles
	for several subroutines.
    SnowPackEnergyBalance.c: (KAC) OTHER
        Certain constant definitions have been moved to vicNl_def.h to 
	standardize the location of model constants.  This is especially 
	useful for eliminating redundant and unused constant definitions.
    SnowPackEnergyBalance.c: (LCB) BLOWING_SNOW update
        New variables are passed to the routine and along to 
	latent_heat_from_snow to add the effects of blowing snow to the
	accumulation and ablation of the snowpack.
    calc_surf_energy_bal.c: (LCB) BLOWING_SNOW update
        Modified to include the effects of blowing snow in the surface 
	energy balance calulations.
    close_files.c: (LCB) LAKE_WETLAND updates
        Now closes the lake and wetland debugging file.
    full_energy.c: (LCB) BLOWING_SNOW updates
        Modified to handle blowing snow.
    full_energy.c: (LCB) LAKE_WETLAND updates
	Modified to output lake model variables during debugging.
    func_surf_energy_bal.c: (LCB) BLOWING_SNOW updates
        Modified to handle blowing snow.
    func_surf_energy_bal.c: (KAC) OTHER
        Modified so that the calculation of sensible heat flux so that 
	occurs in all model versions.  This eliminates a problem in
	WB mode where sensible heat flux was not set to 0, instead it
	showed the cumulative sensible heat flux from the snowpack.
    get_global_param.c: (LCB) BLOWING_SNOW update
        Added BLOWING_SNOW parameter
    get_global_param.c: (Jenny?) NEW_ARNO_TYPE
        Added parameter for reading Bart's new Arno parameters
    get_global_param.c: (LCB) LAKE_WETLAND updates
        Added PRT_LAKE parameter to output lake variables during debugging.
    ice_melt.c: (LCB) LAKE_WETLAND updates
        Modified method by which lake coverage fraction and ice height 
        are updated? ** Check with Laura **
    initialize_global.c: (LCB) LAKE_WETLAND updates
        Added the initialization of PRT_LAKE to the list of debugging flags.
    initialize_global.c: (LCB) BLOWING_SNOW update
        Added initialization of BLOWING_SNOW to the list of global parameters.
    initialize_lake.c: (LCB) LAKE_WETLAND updates
        Made improvements to the initialization process for lakes.  
        ** Check with Laura **
    initialize_model_state.c: (LCB) LAKE_WETLAND updates
        Modified to initialize lake variables.
    initialize_snow.c: (LCB) BLOWING_SNOW update
        Modified to initalize blowing_snow variable.
    initialize_soil.c: (LCB) LAKE_WETLAND updates
        Modified to initialize wetland soil moisture.
    initialize_veg.c: (LCB) LAKE_WETLAND updates
        Modified to get the maximum number of vegetation types passed to 
	it.  This allows the maximum number of vegetation types to include 
	the wetland vegetation fraction when the lake model is active.
    lakes.eb.c: (LCB) LAKE_WETLAND updates
        Modifications were made to improve handling of snow and ice and to 
	make the lake algorithm interact with the wetland algorithm.
    latent_heat_from_snow.c: (LCB) BLOWING_SNOW update
        Modified to handle the effects of blowing snow.
    make_dist_prcp.c: (LCB) LAKE_WETLAND updates
        Modified to allocate vegetation variables for the wetland 
	vegetation class.
    make_in_and_outfiles.c: (LCB) OTHER
        Modified to print notification that the output fluxes file will be 
	in a binary format.
    open_debug.c: (LCB) LAKE_WETLAND updates
	Modified to open lake model debugging file.
    put_data.c: (LCB) LAKE_WETLAND updates
        Updated output of lake variables to reflect algorithm changes.  
    put_data.c: (LCB) BLOWING_SNOW update
        Added output variables for blowing snow algorithm.
    read_lakeparam.c: (LCB) LAKE_WETLAND updates
        Modified to reflect updates to the lake and wetland algorithms.
    read_soilparam.c: (JA) OTHER
        Modified to convert from Bart's new Arno parameters into the 
        standard parameters (Dsmax, Ds, Ws, and c).
    read_vegparam.c: (LCB) BLOWING_SNOW update
        Added code to read in blowing snow parameters.
    snow_melt.c: (LCB) BLOWING_SNOW update
        Modified to handle blowing snow.
    soil_conduction.c (set_node_parameters): (KAC) OTHER
        Modified to correct differences between calculations to determine 
	maximum node moisture and node moisture, so that nodes on the 
	boundary between soil layers are computed the same way for both.
    soil_conduction.c (distribute_node_moisture_properties): (KAC) OTHER
        Modified to check that node soil moisture is less than or equal 
	to maximum node soil moisture, otherwise an error is printed to 
	the screen and the model exits.
    solve_snow.c: (LCB) BLOWING_SNOW update
        Modified to handle the effects of blowing snow.
    surface_fluxes.c: (LCB) BLOWING_SNOW update
        Modified to add the effects of blowing snow.
    user_def.h: (KAC) SPATIAL_SNOW
        Added TraceSnow to indicate the minimum depth of new snow required 
	to reset the snow surface albedo from the ablation to the 
	accumulation curve.
    vicNl.c: (LCB) LAKE_WETLAND updates
        Updated storage of lake water for water balance calculations.
    vicNl.h: (LCB,KAC,JA)
        Updated to reflect model changes.
    vicNl_def.h: (LCB,KAC,JA)
        Updated to reflect model changes.
    water_energy_balance.c: (LCB) LAKE_WETLAND updates
        Updated to reflect changes in algorithm structure.
    water_under_ice.c: (LCB) LAKE_WETLAND updates
        Updated to reflect changes in algorithm structure.
    write_data.c: (LCB) LAKE_WETLAND updates
        Updated output of model for lakes and wetlands algorithm.
        ** No lake variables are output when using LDAS format. **
    write_data.c: (KAC) BLOWING_SNOW update
        Added additional sublimation terms to LDAS and standard snow
        output files.

    Added:
    CalcBlowingSnow.c: (LCB) BLOWING_SNOW update
        Subroutine to compute the effects of blowing snow on snowpack 
	sublimation.
    wetland_energy.c: (LCB) LAKE_WETLAND updates
        Subroutine computes the surface energy balance for exposed 
	wetland vegetation.

-----------------------------------------------------------------------------

    From comparison of LatestAdminSource (updated above) and SOURCE_4.0.4beta

    CalcAerodynamic.c: (KAC) CANOPY_ENERGY_BALANCE update
        This routine was modified to store wind speed, aerodynamics
	resistance, reference height, roughness lengh and displacement
	height for three conditions: (1) snow-free, (2) snow-covered and
	(3) canopy wind speed.
    SnowPackEnergyBalance.c: No Changes Needed
    StabilityCorrection.c: (KAC) OTHER
        Moved definition of G (gravity) to vicNl_def.h to provide a 
	consistant location for all model constants.
    alloc_atmos.c: No Changes Needed
    arno_evap.c: No Changes Needed
    calc_longwave.c: No Changes Needed
    calc_rainonly.c: 
        Why does version 4.0.4 check to see if MAX_SNOW_TEMP <= MIN_RAIN_TEMP,
	rather than justt less than as the new version does?  Was this a bug
	or personal preference?  ** Check This **
    calc_surf_energy_bal.c: (KAC) CANOPY_ENERGY_BALANCE update
        Modified to work with the canopy energy balance updates.
    calc_surf_energy_bal.c: (KAC) SPATIAL_FROST update
        Modified to work with the spatial frost updates.
    calc_surf_energy_bal.c: (KAC) QUICK_SOLVE update
        Modified determine the minimum number of soil thermal nodes to 
	capture the maximum amount of the surface energy flux exchange.
	This reduces solution time with frozen soil significantly, but 
	increases energy balance errors.
    canopy_energy_bal.c: (KAC)
	This include changes for spatial snow and frost plus the new
	canopy energy balance.  No notes indicate fixes made to v4.0.4
	that need to be incorporated.
	- converted internal time step accounting to seconds from hours
	- added check to further restrict evaporative losses if ice content
	  is less than wilting point, otherwise dry soils can evaporate too
	  much if ice is present.
    check_files:
	Added lines for lake file
    check_state_file.c: (KAC) STATE_FILE updates
        Modified to work with new binary state file, fixed in version 4.0.4.
    close_files:
	Added lines for lake file
    compress_files.c: No changes made
	Modified to nice the backgrounded file compression processes.
    dist_prcp.c: (KAC) DIST_PRCP, TREELINE and STATE_FILE updates
	This will have to be carefully merged, as there are several fixes
	and improvements which should be maintained.  Some of these include 
	fixes to the accounting storms, updates for spatial frost and snow,
	updates for the state file, updates to account for the treeline.
	** This will involve some extra effort, especially since distributed
	precipitation variables are now account for in several locations. **
	--> DRY_TIME and STILL_STORM are now computed separatly for vegetation
	    types, this needs to be incorporated into state file.
    frozen_soil.c: (KAC)
	Modified so that soil layer ice content is only calculated if the 
	frozen soil algorithm is implemented and active in the current grid 
	cell.
    full_energy.c: (KAC)
	Added lake algorithm, and canopy energy balance updates.  Nothing
	appears to be needed from v4.0.4.
    func_surf_energy_bal.c: (KAC)
	updated for spatial snow and frost as well as the canopy energy 
	balance.  No important changes in v4.0.4 found.
    get_global_param.c: (KAC)
        Nothing new in v4.0.4, all additions are in the admin code.
    initialize_atmos.c: (KAC)
        Modified to initialize atmospheric pressues in Pa rather than kPa,
	all calculations now use Pa eliminating internal conversions.  Added
	treeline calculations. Added output of forcing file statistics (this
	is in addition to the output of forcing files which was in v4.0.3).
    initialize_global.c: (KAC)
        Added initialization of lake variables, blowing snow variables and
	QUICK_SOLVE for frozen soil energy fluxes.
    initialize_model_state.c: (KAC)
        Updated for canopy energy balance, new model state file, and lakes 
	and wetlands algorithm.
    initialize_new_storm.c: (KAC)
        Modified to work with spatial frozen soil.
    initialize_snow.c: (KAC)
        Modified to initialize blowing snow and spatial snow algorithm 
        variables.
    make_dist_prcp.c: (LCB)
        Modified to include wetland variables in the new structures.
    make_dmy.c:
        Only white space changes
    make_in_and_outfiles.c: (KAC)
        Lake model file control was added
    mtclim42_vic.c: (KAC) 
        Modified to compute all pressures in Pa rather than kPa.
    mtclim42_vic.h: (KAC)
        Definition of EPS constant moved to vicNl_def.h
    mtclim42_wrapper.c: (KAC)
        Changed calls to vicrerror to calls to nrerror.  Now handles 
	pressure in Pa, rather than kPa.
    open_state_file.c: (KAC)
        Updated to handle binary state file.
    penman.c: (KAC)
        Now handles pressure in Pa, rather than kPa.
    prepare_full_energy.c: (KAC)
        modified so that ice content is set to zero unless the frozen soil 
	algorithm is implemented and active in the current grid cell. 
    put_data.c: (KAC)
        modified to incorporate the effects of treeline calculations, the
	lake and wetland algorithm, spatial snow and frost algorithms, and 
	the canopy energy balance.
    read_atmos_data.c: (KAC)
        v4.0.4 multiplies nrecs by NF to check for a long enough file,
        v.4.1.0 multiples by dt.  I think new version is probably correct.
    read_initial_model_state.c: (KAC)
        Updated to handle binary state file.  Need to address storm state
	variables!
    read_snowband.c: (KAC)
	Modified to allocate treeline variable.  Will need to update to fix
        lake model/snow band problem!
    read_soilparam.c: (KAC, JA)
        Modified to read in spatial snow and frost parameters and to read in 
	Bart's new Arno parameter set.  
    read_soilparam_arc.c: (KAC)
        Modified to read in spatial snow and frost parameters and to read in 
	Bart's new Arno parameter set.
    read_vegparam.c: (DP,KAC)
        Modified code to update Wdmax based on LAI values read in for the 
	current grid cell.  If LAI is not obtained from this function, then 
	the values cacluated in read_veglib.c are left unchanged.
    redistribute_during_storm.c: (KAC)
        Modified to work with distributed snow and frozen soil.
    runoff.c: (KAC)
        Modified to handle spatial snow and soil frost.
    snow.h: (KAC)
        Added minimum SWQ for computing snow pack energy balance.  Added
        coefficients of shortwave attenuation through the snowpack.  
	Removed minimum SWQ for full coverage, now included in the soil 
	files.
    snow_intercept.c: (KAC)
        Modified to handle new variables required to close the canopy 
	energy balance.
    snow_melt.c: (KAC)
        Modified to handle partial snow cover.  Modified to assure that 
	ground heat flux is used properly in the snow surface energy 
	balance as well as imporving the handling of energy fluxes for 
	partial snow cover.  Modified to handle blowing snow.
    snow_utility.c: (KAC)
        Moved definition of G to vicNl_def.h.
    soil_conduction.c: (KAC)
	set_node_parameters:
        Modified to correct differences between calculations to determine 
	maximum node moisture and node moisture, so that nodes on the 
	boundary between soil layers are computed the same way for both.
	distribute_node_moisture_propertes:
	Modified to check that node soil moisture is less than or equal to 
	maximum node soil moisture, otherwise an error is printed to the 
	screen and the model exits.
	estimate_layer_ice_content:
	Modified to find ice content in spatial frost bands.
    solve_snow.c: (KAC)
	- Added partial snow cover and advection of sensible heat from local 
	bare patches.
	- Modified to pass the minimum depth of full snow cover as a variable 
	in soil_con rather than a globally defined constant.
	- Fixed check of new snow accumulation for setting understory flag to 
	use snowfall[WET] not snowfall.
	- Set MELTING flag to maintain melting albedo curve even during brief 
	periods of refreezing, until a snowfall exceeds SnowThres.
	- Modified to handle the effects of blowing snow.
	- Modified to handle closed canopy energy balance.
    surface_fluxes.c: (KAC)
        - Modified to handle partial snow cover.
        - Modified to iterate a solution for the exchange of energy between 
	the snowpack and the ground surface.
        - Modified to add the effects of blowing snow.
        - Fixed indexing problem for sub-daily snow model within daily water 
	balance VIC: hour (now hidx) is incremented by 1 rather than the 
	sub-daily time step, so the atmospheric forcing data is now properly 
	indexed.
	- Indexing fix sent SNOW_STEP to calc_surf_energy_bal rather than the 
	model time step, meaning that without snow the evaporation was 
	computed for SNOW_STEP hours rather than a full day.  This was fixed 
	by introducing step_inc to index the arrays, while step_dt keeps 
	track of the correct time step.
    svp.c: (KAC)
        Modified internal calculations to use kPa rather than Pa for 
	consistancy
    user_def.h: (KAC)
        Added options for lake model, closing the canopy energy balance,
	computing the treeline, computing statistics from the input 
	forcings, and controls for spatial snow and frost.
    vicNl.c: (KAC)
        - Added controls for lake model.
	- Updated storage of lake water for water balance calculations.
	- Modifed to add AboveTreeLine to soil_con_struct so that the model 
	can make use of the computed treeline.
	- Modified to initialize storm parameters using the state file.
        - Modified to start the model by skipping records until the state 
	file date is found.  This replaces the previous method of modifying 
	the global file start date, which can change the interpolation of 
	atmospheric forcing data.
	- Modified to store wet and dry fractions when intializing water 
	balance storage.  This accounts for changes in model state 
	initialization, which now stores wet and dry fractions rather than 
	just averaged values.
    vicNl.h: (KAC)
        Added definitions for treeline and forcing stats subroutines and
	made some needed modifications to subroutine definitions due to
	other changes, however, other changes still need to be made.
    vicNl_def.h: (KAC)
        Modified to handle lake algorithm, blowing snow, spatial snow and
	frost, and treeline calculations.
    write_data.c: (KAC)
	- Made hour a variable in all output data file formats even if the 
	model is run at a daily time step.  Also modified all output files 
	to account for new variables introduced by the spatial frost and 
	snow algorithms, the lake algorithm and the PILPS 2e study.
	- Added energy fluxes to snow band output files.
	- Updated output of model for lakes and wetlands algorithm.  Added 
	output of blowing snow sublimation to LDAS and standard snow output 
	files.  ** No Lake Variables are included in the LDAS output format. 
	**
	- Modified LDAS SWQ output, so that it is multiplied by 10 instead 
	of 100 before being converted to a short integer.  This reduces 
	stored value precision to 0.1, but increases the maximum storable 
	SWQ, which was exceeded in previous LDAS simulations.
        - Eliminated different formats between energy balance and water
	balance model output.
    write_debug.c: (KAC)
        Modified to work with closed canopy energy balance.
    write_forcing_file.c: (KAC)
	Modified to output pressures, which are handled internally in kPa, 
	as Pa for backward compatability.
    write_layer.c: (KAC)
        Modified to handle spatial soil frost.
    write_model_state.c: (KAC)
	- Rewritten to handle updates to vicNl_def.h and to write the file 
	as binary to minimize write time and differences with simulations 
	started with the state file.
	- Model is now restarted with the correct values for mu and 
	LAST_STORM.
        ** Still need to account for differences with distributed 
	precipitation flags, which now differ with vegetation - also need
	to check storage of lake variables, blowing snow variables and 
	MELTING flag.

-----------------------------------------------------------------------------

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
	    at the beginning of each time step, so storing its value
	    is unnecessary.

-----------------------------------------------------------------------------

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

-----------------------------------------------------------------------------

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
	     actually LAI values used for the grid cell.  Values for Wdmax 
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
	     model state was saved.  This means that calculations to yield
	     sub-daily metrological forcings from daily forcings will produce 
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
	     

-----------------------------------------------------------------------------

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

-----------------------------------------------------------------------------

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

-----------------------------------------------------------------------------

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
             

-----------------------------------------------------------------------------

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
