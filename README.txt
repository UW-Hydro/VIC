# README.txt - Release notes
# $Id$
#-----------------------------------------------------------------------

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
***** Description of changes from VIC 4.0.5 to VIC 4.0.6 beta r1 *****
--------------------------------------------------------------------------------


New Features:
-------------

Added sample global param file to distribution

        Files affected:

        global.param.sample (new)
										TJB


Flexible output configuration

	Files affected:

	Makefile
	calc_water_energy_balance_errors.c (new)
	close_files.c
	display_current_settings.c
	dist_prec.c
	get_global_param.c
	global.param.sample
	initialize_atmos.c
	make_in_and_outfiles.c
	output.LDAS_OUTPUT.template (new)
	output.OPTIMIZE.template (new)
	output.OUTPUT_FORCE.template (new)
	output.TRADITIONAL.template (new)
	output_list_utils.c (new)
	parse_output_info.c (new)
	put_data.c
	set_output_defaults.c (new)
	user_def.h
	vicNl.c
	vicNl.h
	vicNl_def.h
	vicerror.c
	write_data.c
	write_forcing_file.c

	Description:

	In earlier versions of VIC, the set of output files and their contents
	were hard-coded.  A few settings in user_def.h (OPTIMIZE and LDAS_OUTPUT)
	allowed the user to change the contents of the output files, but this
	has not been enough to accomodate the various needs of users.  Users
	inevitably have had to modify write_data.c to produce the type of output
	they want, and when VIC is updated, they must merge their changes into
	the new version of the code.

	VIC 4.0.6 now allows the user to specify exactly which output files to
	create and which variables to store in each file.  This way, users can
	save space by only writing those variables that are useful, and will be
	less likely to need to maintain a private version of the code to do this.


	Main points:

	1. Output file names and contents can be specified in the global param
	   file (see below).

	2. If you do not specify file names and contents in the global param
	   file, VIC will produce the same set of output files that it has
	   produced in earlier versions, namely "fluxes" and "snow" files, plus
	   "fdepth" files if FROZEN_SOIL is TRUE and "snowband" files if
	   PRT_SNOW_BAND is TRUE.  These files will have the same contents and
	   format as in earlier versions.

	3. The OPTIMIZE and LDAS_OUTPUT options have been removed.  These
	   output configurations can be selected with the proper set of
	   instructions in the global param file.  (see the output.*.template
	   files included in this distribution for more information.)

	4. If you do specify the file names and contents in the global param file,
	   PRT_SNOW_BAND will have no effect.

	To specify file names and contents in the global param file, one
	should use the following format:

	  (typical global param file contents here...)

	  # Output File Contents
	  N_OUTFILES	<n_outfiles>

	  OUTFILE	<prefix>	<nvars>
	  OUTVAR	<varname>	[<format>	<type>	<multiplier>]
	  OUTVAR	<varname>	[<format>	<type>	<multiplier>]
	  OUTVAR	<varname>	[<format>	<type>	<multiplier>]

	  OUTFILE	<prefix>	<nvars>
	  OUTVAR	<varname>	[<format>	<type>	<multiplier>]
	  OUTVAR	<varname>	[<format>	<type>	<multiplier>]
	  OUTVAR	<varname>	[<format>	<type>	<multiplier>]

        where
		<n_outfiles> = number of output files
		<prefix>     = name of the output file, NOT including latitude
		               and longitude
		<nvars>      = number of variables in the output file
		<varname>    = name of the variable (this must be one of the
		               output variable names listed in vicNl_def.h.)

                <format>, <type>, and <multiplier> are optional.  For a given
		variable, you can specify either NONE of these, or ALL of
		these.  If these are omitted, the default values will be used.

		<format>     = (for ascii output files) fprintf format string,
		               e.g.
		                 %.4f = floating point with 4 decimal places
		                 %.7e = scientific notation w/ 7 decimal places
				 *    = use the default format for this variable
		<type>       = (for binary output files) data type code.
		               Must be one of:
		                 OUT_TYPE_DOUBLE = double-precision floating point
		                 OUT_TYPE_FLOAT  = single-precision floating point
		                 OUT_TYPE_INT    = integer
		                 OUT_TYPE_USINT  = unsigned short integer
		                 OUT_TYPE_SINT   = short integer
		                 OUT_TYPE_CHAR   = char
				 *               = use the default type
		<multiplier> = (for binary output files) factor to multiply
		               the data by before writing, to increase
		               precision.
				 *    = use the default multiplier for this variable

	Here's an example.  To specify 2 output files, named "wbal" and
	"ebal", and containing water balance and energy balance terms,
	respectively, you could do something like this:

	  N_OUTFILES	2

	  OUTFILE	wbal	6
	  OUTVAR	OUT_PREC
	  OUTVAR	OUT_EVAP
	  OUTVAR	OUT_RUNOFF
	  OUTVAR	OUT_BASEFLOW
	  OUTVAR	OUT_SWE
	  OUTVAR	OUT_SOIL_MOIST

	  OUTFILE	ebal	7
	  OUTVAR	OUT_NET_SHORT
	  OUTVAR	OUT_NET_LONG
	  OUTVAR	OUT_LATENT
	  OUTVAR	OUT_SENSIBLE
	  OUTVAR	OUT_GRND_FLUX
	  OUTVAR	OUT_SNOW_FLUX
	  OUTVAR	OUT_ALBEDO

	Since no format, type, or multiplier were specified for any variables, VIC will
	use the default format, type, and multiplier for the variables.

	If you wanted scientific notation with 10 significant digits for ALBEDO,
	you could do the following:

	  OUTVAR	OUT_ALBEDO		%.9e	*	*

	Note that even if you only want to specify the format, you must supply a value
	in the type and multiplier columns as well.  This can be "*" to indicate the
	default value.  Similarly, if you only want to specify the type (e.g. as a double),
	you would need to do something like:

	  OUTVAR	OUT_ALBEDO		*	OUT_TYPE_DOUBLE		*


	Date variables:

	For typical output files, the date is always written at the beginning of
	each record.  This will consist of the following columns:
	  year month day hour
	For daily output timestep, "hour" is not written.

	If BINARY_OUTPUT is TRUE, these will all be written as type int (OUT_TYPE_INT).

	If OUTPUT_FORCE is TRUE, the date will NOT be written.


	Multiple-valued variables:

	Since variables like SOIL_MOIST have 1 value per soil layer, these variables
	will be written to multiple columns in the output file, one column per soil
	layer.  Other multiple-valued variables are treated similarly.


	Snow band output:

	To specify writing the values of variables in each snow band,  append
	"_BAND" to the variable name (this only works for some variables - see
	the list in vicNl_def.h).  If you specify these variables, the value of
	the variable in each band will be written, one band per column.  For
	example, for a cell having 2 snow bands:

	  OUTVAR	OUT_SWE_BAND
	  OUTVAR	OUT_ALBEDO_BAND

	will result in an output file containing:

	  year month day (hour) swe[0] swe[1] albedo[0] albedo[1]

										TJB


ALMA-compliant input and output

	Files affected:

	SnowPackEnergyBalance.c
	display_current_settings.c
	dist_prec.c
	full_energy.c
	get_force_type.c
	get_global_param.c
	global.param.sample
	initialize_atmos.c
	initialize_global.c
	initialize_model_state.c
	output.OUTPUT_FORCE.ALMA.template (new)
	output.PILPS-2E.ALMA.template (new)
	output.TRADITIONAL.template
	output_list_utils.c
	put_data.c
	surface_fluxes.c
	vicNl.c
	vicNl.h
	vicNl_def.h
	write_forcing_file.c

	Description:

	This change allows VIC to handle input and output in a manner
	compliant the the ALMA convention used in the PILPS-2e experiment
	(http://www.lmd.jussieu.fr/~polcher/ALMA/).

	ALMA INPUT:

	VIC now accepts the following new ALMA-compliant input forcings
	in addition to the forcings that it already accepts:
		SNOWF     snowfall rate (kg/m^2s)
		RAINF     rainfall rate (kg/m^2s)
		CRAINF    convective rainfall rate (kg/m^2s)
		LSRAINF   large scale rainfall rate (kg/m^2s)
		QAIR      specific humidity (kg/kg)
		WIND_E    zonal wind speed (m/s)
		WIND_N    meridional wind speed (m/s)
		TAIR      air temperature per time step (K)
		PSURF     atmospheric pressure (Pa)

	When giving VIC ALMA-compliant input files, you must be sure to use
	the names given above in the forcing section of your global parameter
	file.

	Instead of the existing PREC (precipitation per timestep in mm), you
	can now specify SNOWF and RAINF (snowfall and rainfall rates, both in
	mm/s).  VIC will simply add these two quantities together, multiply
	by the forcing interval, and treat their sum the same way it treats
	PREC.

	An alternative to supplying RAINF is to supply CRAINF (convective
	rainfall rate, mm/s) and LSRAINF (large-scale rainfall rate, mm/s).
	VIC will add these two quantities together to get RAINF.

	Instead of the existing WIND, alternatively you can specify WIND_E
	and WIND_N (zonal and meridional wind speed, m/s).  VIC will simply
	compute WIND = sqrt(WIND_E**2+WIND_N**2).

	TAIR has units of K, while the existing AIR_TEMP is in C.  Similarly,
	PSURF is in Pa, while PRESSURE is in kPa.  VIC will convert these
	to AIR_TEMP and PRESSURE after reading them in.

	More information is available on ALMA forcing variables at:
	  http://www.lmd.jussieu.fr/~polcher/ALMA/convention_input_3.html

	ALMA OUTPUT:

	If the user sets ALMA_OUTPUT=TRUE in the global parameter file, then
	VIC will convert its output variables to ALMA-compliant forms.  The
	majority of the changes are changes of units.  Moisture fluxes are
	changed from VIC's standard (mm accumulated over the time step) to
	the average flux rate (mm/s).  Temperatures are converted from C
	to K.  More information on ALMA output variables is available at:
	  http://www.lmd.jussieu.fr/~polcher/ALMA/convention_output_3.html

	In addition, several more variables have been added to the list of
	available output variables.  See vicNl_def.h for the complete list
	of available output variables.						TJB


Aggregation of output variables

	Files affected:

	display_current_settings.c
	get_global_param.c
	global.param.sample
	output_list_utils.c
	put_data.c
	vicNl_def.h
	write_data.c

	Description:

	VIC can now aggregate the output variables to a user-defined
	output interval, via the OUT_STEP setting in the global parameter
	file.  Currently, the largest output interval allowed is 24 hours,
	so this option is only useful for simulations running at sub-daily
	time steps.								TJB


Cleanup of structures holding filenames and file pointers

	Files affected:

	check_files.c
	check_state_file.c
	close_files.c
	display_current_settings.c
	dist_prec.c
	get_global_param.c
	initialize_model_state.c
	make_in_and_outfiles.c
	open_state_file.c
	read_initial_model_state.c
	vicerror.c
	vicNl.c
	vicNl_def.h
	vicNl.h
	write_model_state.c

	Description:

	1. Merged infiles and outfiles structs into filep_struct.
	2. Merged builtnames into filenames struct.
	3. Renamed infiles.statefile to filep.init_state
	4. Moved global.statename to filenames.statefile.
	5. Added f_path_pfx[] to the filenames_struct, to store
	   the path and prefix of forcing files.  Now, forcing[]
	   only stores the full forcing file names.				TJB


More complete set of supported input variables

        Files affected:

        display_current_settings.c
        get_force_type.c
        get_global_param.c
        global.param.sample
        initialize_atmos.c
        initialize_global.c
        vicNl_def.h

        Description:

        Added REL_HUMID (relative humidity, as a fraction), CSNOWF
        (convective snowfall) and LSSNOWF (large-scale snowfall) to the
        list of supported met input variables.  Added the ALMA_INPUT
        option, which causes temperatures to be interpreted as Kelvin,
        pressures as kPa, and moisture fluxes as rates (mm/s) instead of
        accumulated fluxes (mm/timestep).  This allowed us to remove TAIR
        and PSURF from the list of supported met input variables, since
        AIR_TEMP and PRESSURE provide this ability when ALMA_INPUT is
        TRUE.  This also makes input variable specification more
        consistent with output variable specification.				TJB


Atmos_data arrays are always allocated dynamically now.

        alloc_atmos.c
        vicNl_def.h

        Description:

        Under some circumstances (OPTIMIZE and LINK_DEBUG options),
        the arrays in the atmos_data structure were hard-wired to have
        25 elements (large enough to store 24 hourly values and 1 daily
        value).  This wasted memory and has been abandoned.  Now, the
        arrays in atmos_data are always allocated dynamically to be only
        big enough to store exactly the number of necessary elements.		TJB


Optional headers for output and input files

        Files affected:

        display_current_settings.c
        get_global_param.c
        global.param.sample
        initialize_global.c
        Makefile
        read_atmos_data.c
        vicNl.c
        vicNl_def.h
        vicNl.h
        write_header.c (new)

        Description:

        Added the PRT_HEADER option to the global parameter file.  If
        this is set to TRUE, VIC will insert a short header into its
        output files, describing the time step, start date/time,
        variables and units included in the file.

        For ascii files, the output header has the following format:

                # NRECS: (nrecs)
                # DT: (dt)
                # STARTDATE: yyyy-mm-dd hh:00:00
                # ALMA_OUTPUT: (0 or 1)
                # NVARS: (Nvars)
                # VARNAME    VARNAME   VARNAME   ...

        where
           nrecs       = Number of records in the file
           dt          = Time step length in hours
           start date  = Date and time of first record of file
           ALMA_OUTPUT = Indicates units of the variables; 0 = standard VIC
                         units; 1 = ALMA units
           Nvars       = Number of variables in the file, including date
                         fields

        For binary files, the output header has the following format:

            // Data        Stored As           Comment
            //
            // Identifier  (unsigned short)*4  0xFFFF, repeated 4 times
            // Nbytes      (unsigned short)*1  Number of bytes in the header,
            //                                 INCLUDING THE IDENTIFIER
            //
            // Part 1: Global Attributes
            // Nbytes1     (unsigned short)*1  Number of bytes in part 1
            // nrecs       (int)*1             Number of records in the file
            // dt          (int)*1             Time step length in hours
            // startyear   (int)*1             Year of first record
            // startmonth  (int)*1             Month of first record
            // startday    (int)*1             Day of first record
            // starthour   (int)*1             Hour of first record
            // ALMA_OUTPUT (char)*1            0 = standard VIC units; 1 = ALMA units
            // Nvars       (char)*1            Number of variables in the file,
            // including date fields
            //
            // Part 2: Variables
            // Nbytes2     (unsigned short)*1  Number of bytes in part 2
            // For each variable, the following fields: { len varname type mult }
            //   len       (char)*1            Number of characters in varname
            //   varname   (char)*len          Variable name
            //   type      (char)*1            Code identifying variable type
            //   mult      (float)*1           Multiplier for variable


        To accomodate input forcing files that might have been produced via
        VIC's OUTPUT_FORCE option, and therefore could contain a header,
        read_atmos_data.c has been modified to detect and skip headers that
        follow the formats outlined above.					TJB


Variable TYPE specifications for binary-format output files in the global parameter
file must match the strings listed in vicNl_def.h.

        Files affected:

        parse_output_info.c

        Description:

        When listing output variables in the global parameter file, if the
        output file format is binary, the variable data TYPE must match the
	string from vicNl_def.h exactly, e.g. "OUT_TYPE_INT" rather than
	just "INT".								TJB


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
        the beginning of the new forcing file.					TJB


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
        lines are present, VIC will give an error.				TJB


Bug Fixes:
----------

Aerodynamic resistance incorrect in output fluxes file

        Files affected:
        SnowPackEnergyBalance.c, calc_surf_energy_bal.c, full_energy.c,
        put_data.c, snow_intercept.c, snow_melt.c, solve_snow.c,
        surface_fluxes.c, vicNl.h, vicNl_def.h

        Description:
        In 4.0.5 and earlier, the aerodynamic resistance written to the output
        fluxes file rarely reflected the value actually used in flux computations.
        This has been fixed.

        VIC actually computes an array of 3 different aerodynamic resistances,
        as follows:
          aero_resist[0] : over vegetation or bare soil
          aero_resist[1] : over snow-filled overstory
          aero_resist[2] : over snow pack
        VIC determines which element of the array to use depending on the current
        vegetation type and whether snow is present.  In addition, in most cases,
        VIC applies a stability correction to this aerodynamic resistance before
        using it in flux computations.  Furthermore, when the current vegetation
        tile contains overstory and snow is present on the ground, aero_resist[2]
        is used for snow pack flux computations and either aero_resist[1] or
        aero_resist[0] is used for canopy flux computations, meaning that two
        different aerodynamic resistances are in use in the same time step.

        However, VIC 4.0.5 always wrote the uncorrected value of aero_resist[0]
        to the fluxes file.

        In 4.0.6, the value written to the fluxes file is the actual value used
        in flux computations, including any corrections that were applied.  In
        the case mentioned above in which two different aerodynamic resistances
        are in use at the same time, the one used for the snow pack is
	written.								TJB


Aerodynamic resistance not correctly aggregated for output

        Files affected:
        put_data.c
        vicNl_def.h

        In previous releases, aerodynamic resistance
        (out_data->aero_resist) was aggregated by a simple
        area-weighted average over veg tiles.  This led to an
        aggregate value that was not the true effective resistance of
        the entire grid cell.  Since evaporation is proportional to
        1/aero_resist, it is (1/aero_resist), or the aerodynamic
        conductivity, that should be averaged over the grid cell.
        Therefore, a new variable, out_data.aero_cond, was created for
        the purposes of aggregation.  After aggregation,
        out_data.aero_resist is computed as 1/out_data.aero_cond.

        The effect of the change is most pronounced when there is a
        large range of values of aerodynamic resistance among the veg
        tiles in a grid cell.  The effective aerodynamic resistance,
        computed the new way, will tend to be smaller than the old way
        when the values cover a wide range.  However, the effective
        aerodynamic resistance will never be smaller than the smallest
        value among the various veg tiles in the cell.				TJB


Skipping deactivated cells in binary state file

        Files affected:
        write_model_state.c
        
        (Port from 4.1.0) Changed calculation of Nbytes in binary
        state file to account for bare soil values (extra veg class
        per grid cell). Without this fix, attempts to skip grid cells
        fail.									GCT


Fix reading/writing of state files when QUICK_FLUX is TRUE

        Files affected:

        get_global_param.c
        make_dist_prcp.c
        make_energy_bal.c
        vicNl.c
        vicNl.h

        In previous releases, if QUICK_FLUX=TRUE, the number of soil
        thermal nodes is changed after being recorded in the state
        file.  The resulting mismatch between Nnodes in the state
        file header and the actual number of node temperatures
        recorded in the state file prevents VIC from being able to
        read the state file.  This has been fixed.				GCT


Fixed bug in error trapping when INIT_STATE filename matches SAVE_STATE filename

        Files affected:
        get_global_param.c, open_state_file.c

        Previously, the checks for matches would occur in get_global_param even
        though the SAVE_STATE filename wasn't created until open_state_file.
        The output name setting was moved from open_state_file to
	get_global_param.							GCT


Added checks for range/valid month days

        Files affected:
        get_global_param.c

        In previous versions the user could set a non-valid date for STATE files
        and the model would run without writing to STATE file. Code now checks
        for valid date.								GCT


Replace %i with %d in scanf statements.

        Files affected:
        check_state_file.c
        get_global_param.c
        read_arcinfo_ascii.c
        read_initial_model_state.c
        read_snowband.c
        read_soilparam.c,
        read_veglib.c

        Having %i in fscanf statements was causing input values of
        "08" to be interpreted as octal rather than decimal.  These
        instances of %i have been replaced with %d.				GCT


ARNO_PARAMS global parameter option changed to BASEFLOW

        Files affected:
        display_current_settings.c
        get_global_param.c
        global.param.sample
        initialize_global.c
        read_soilparam.c
        read_soilparam_arc.c
        vicNl_def.h

        Changed the name of the ARNO_PARAMS global parameter option to
        BASEFLOW.  The meaning of the ARNO_PARAMS option was actually
        opposite to its name: when ARNO_PARAMS was FALSE, VIC would
        interpret the first four parameters in the soil parameter file
        to be the standard ARNO soil parameters Ds, Dsmax, Ws, and c,
        while when ARNO_PARAMS was TRUE, VIC would interpret the first
        four parameters to be d1, d2, d3, and d4, the soil parameters
        used in Nijssen et al. (2001).  The new option now can take
        values of "ARNO" and "NIJSSEN2001".  When BASEFLOW == NIJSSEN2001,
        VIC assumes the soil parameter file contains d1, d2, d3, and d4.
        When BASEFLOW == ARNO, VIC assumes the soil parameter file
        contains Ds, Dsmax, Ws, and c.						TJB


Allow NO_FLUX in addition to NOFLUX in global.param.file

        Files affected:
        get_global_param.c

        The option NOFLUX has a syntax (ie, the missing underscore) that is
        inconsistent with other FLUX options. The change will allow users to
        enter either string.							GCT


Skip reading/writing of snow band for areafract <= 0 

        Files affected:
        read_initial_model_state.c
        write_model_state.c

        This will reduce the size of the statefile.				GCT


Changed argument order in fread, fwrite statements.

        Files affected:
        check_state_file.c
        open_state_file.c
        read_atmos_data.c
        read_initial_model_state.c
        write_data.c
        write_forcing_file.c
        write_model_state.c

        Statements had arguments with ...1, sizeof()....Those were changed to
        ...sizeof(), 1, ...							GCT


OUTPUT_FORCE option does not close output files properly

        Files affected:
        check_files.c
        close_files.c
        get_global_param.c
        initialize_atmos.c
        read_soilparam.c
        read_soilparam_arc.c
        vicNl.c

        The OUTPUT_FORCE compile-time option excluded the code that closed the
        output files, causing the output disaggregated forcing files to end
        abruptly before the end of the simulation period (due to the final
        chunk of the buffer not being flushed).  In addition, the OUTPUT_FORCE
        option included the validation code for the entire soil parameter
        file, whose validation is not necessary for this option (and can prevent
        this option from working.  This fix remedies these problems, by
        including the necessary file-closing code and excluding the unnecessary
        soil-parameter-file-checking code.					TJB


Bus error in cells that have bare soil

	Files affected:
	initialize_model_state.c
        
	An error in the flexible output configuration feature resulted in a bus
	error in cells that have non-zero bare soil fractions.  This has been
	fixed.									TJB


Incorrect sub-daily temperature interpolation when referencing GMT instead of
local time

	Files affected:
	calc_air_temperature.c
        
	Temperature interpolation didn't account for case in which min or max
	temperature could cross the boundary of the current day.  This can
	happen when referencing GMT instead of local time, for cells far away
	from 0 E longitude.  This has been fixed.				TJB


AIR_TEMP was not being allowed as an input forcing variable

	Files affected:
	initialize_atmos.c

	Description:
	A typo in an "if" statement prevented AIR_TEMP from being allowed as
	a valid input forcing variable at any time step.  This has been
	fixed.									TJB


Bug: undeclared variable i in write_forcing_file.c

        Files affected:

        write_forcing_file.c

        Description:

        Index i was not declared.  This has been fixed.				TJB


If() statements in get_force_type() fail for some global parameter files

        Files Affected:

        get_force_type.c

        Description:

        Removed all of the if statements
          if(param_set.FORCE_FORMAT[file_num]==BINARY)
        since this ended up requiring that the "FORCE_FORMAT BINARY"
        line appear in the global parameter file before the list of
        forcing variables in order to work.  Since the sscanf()
        performs proper parsing regardless of ASCII (which doesn't
        have SIGNED or MULTIPLIER fields) vs. BINARY, we have removed
        the if() statements altogether.					TJB


Aggregation methods of some variables not set properly.

        Files affected:

        output_list_utils.c

        Description:

        Corrected AGG_TYPE definitions for miscellaneous
        output variables; re-organized the code to make
        it easier to debug.						TJB


Fixes for memory leaks and variable initialization.

        Files affected:

        free_dist_prcp.c
        get_global_param.c
        make_dmy.c
        mtclim42_vic.c
        output_list_utils.c
        parse_output_info.c
        read_veglib.c
        vicNl.c
        vicNl.h

        Description:

        Miscellaneous fixes for memory leaks and variable
	initialization.							TJB


Memory errors for ARC_SOIL=TRUE and OUTPUT_FORCE=TRUE

	Files Affected:

	get_force_type.c
	get_global_param.c
	initialize_global.c
	read_soilparam.c
	read_soilparam_arc.c
	vicNl.c

	Description:

	Memory errors would occur when ARC_SOIL=TRUE and
	OUTPUT_FORCE=TRUE.  In addition, the output files
	would not contain sufficient contents due to not
	closing properly.						TJB


Fixed fread checks

	Files affected:

	read_initial_model_state.c

	Fixed fread checks to make sure correct number of
	items were read in rather than the size of the item
	read in.							JCA
	port from 4.1.0_r4						GCT


Bug fix for previous bug fix to dt_baseflow calculation.

	Files Affected:

	runoff.c

	Description:

	Fixed bug arising from earlier fix to dt_baseflow
	calculation.  Earlier fix took residual moisture
	into account in the linear part of the baseflow eqn,
	but not in the non-linear part.  Now we take residual
	moisture into account correctly throughout the whole
	equation.							TJB


Removed logic that reset resid_moist[i].

	Files Affected:

	runoff.c

	Description:

	Removed logic that reset resid_moist[i].  Previously,
	resid_moist[i] was reset to 0 for i > 0 when
	resid_moist[0] == 0.  Such resetting of soil properties
	was deemed unnecessary and confusing, since VIC would end
	up using different residual moisture values than those
	specified by the user.  If a user truly wants to specify
	residual moisture in all layers to be 0, the user should
	set these explicitly in the soil parameter file.		TJB


Seg fault when trying to read initial state file

	Files Affected:

	vicNl.c

	Description:

	Fixed typo in call to check_state_file().  Was assigning
	init_state file pointer to filep.statefile; now assigns
	pointer to filep.init_state.					TJB


Model aborts when TIME_STEP = 24 and STARTHOUR not specified.

	Files Affected:

	get_global_param.c

	Description:

	Added validation of dt, start date, end date, and nrecs.	TJB


Output file headers contain "hour" field despite output dt == 24 hours.

	Files affected:

	write_header.c

	Description:

	Replaced all instances of global.dt with global.out_dt,
	since out_dt is the time interval used in the output files.	TJB


Liquid soil moisture sometimes falls below residual.

	Files affected:

	runoff.c

	Description:

	Fixed the checks on the lower bound of soil moisture.
	Previously, the condition was
	  (moist[lindex]+ice[lindex]) < resid_moist[lindex]
	which led to liquid soil moisture falling below residual
	during winter conditions.  This has been changed to
	  moist[lindex] < resid_moist[lindex]
	to eliminate these errors and make the logic consistent
	with the rest of the code.					TJB


Variable "moist" in runoff() has different meaning than in other functions.

	Files affected:

	runoff.c

	Description:

	Renamed all *moist* variables to *liq* if they only refer
	to liquid soil moisture.  This makes the logic much easier
	to understand.							TJB


Soil moisture drops below residual for sub-daily time step interval and
FULL_ENERGY or FROZEN_SOIL = TRUE.

	Files Affected:

	runoff.c

	Description:

	Fixed the checks on the lower bound of soil moisture.  Previously,
	the lower bound on soil moisture was
	  (liquid + ice >= residual moisture)
	and the way this bound was enforced was to reduce baseflow to
	bring total liquid + ice content back up to residual, but this
	compensation was limited so that the compensation would never
	be larger than baseflow (i.e. the compensation would never create
	negative baseflow).

	However, this condition had two main problems: 1. it allowed liquid
	moisture to fall to very low values (due to evap being over-estimated
	in arno_evap() and transpiration(), and/or Q12 being over-estimated
	earlier in runoff() due to bad numerics) in the presence of ice, and
	2.  it had limited ability to recover from these low values because
	baseflow (already small) wasn't allowed to be reduced below 0.

	This behavior has been replaced with the following conditions: For
	unfrozen soil, the new lower bound is
	  (liquid >= resid_moist)
	while for frozen soil, the new error condition is
	  (liquid >= min_liq_fraction * resid_moist)
	where
	  min_liq_fraction = the value returned by maximum_unfrozen_water()
	                     for a unit maximum moisture content.

	If overestimates of evap or drainage do occur, baseflow is allowed
	to be reduced below 0 to bring liquid water back up to the lower
	limit.  Then, further down in the code, if baseflow is negative,
	bottom-layer evap is reduced by (-baseflow) and baseflow is set to 0.	TJB



-------------------------------------------------------------------------------
***** Description of changes from VIC 4.0.4 to VIC 4.0.5 *****
--------------------------------------------------------------------------------


New Features:
-------------

"-o" command-line option and display of run/compile-time options

	Files affected:
	cmd_proc.c, display_current_settings.c, get_global_param.c, global.h,
	vicNl.c, vicNl.h, vicNl_def.h

	Description:
	In VIC 4.0.5, if VERBOSE is TRUE, all compile-time options (from
	user_def.h) and run-time options (from your global parameter file) are
	displayed for you at the beginning of a model run.  If you are saving
	your model output to a log file, this information will be stored with
	your log, so you will have a record of the option settings that produced
	the results of that particular model run.
	In addition, added the "-o" option, to display current compile-time
	options.  One benefit of this option is that you can see what the
	options in user_def.h were set to when vicNl was compiled, without
	having to start a model run.  Since your version of user_def.h may have
	changed since you compiled vicNl, this is the most reliable way to see
	what these options are set to.

	Example:

	gen.hydro.washington.edu 41: vicNl -o

	***** VIC Version 4.0.5 - Current Model Settings *****

	COMPILE-TIME OPTIONS (set in user_def.h)
	----------------------------------------

	Output to Screen:
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
	LOW_RES_MOIST           FALSE
	QUICK_FS                FALSE

	Maximum Array Sizes:
	MAX_BANDS               10
	MAX_FRONTS               3
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
	In VIC 4.0.4 and earlier, updating a .h file and recompiling VIC would
	not result in recompilation of files that depend on the .h file, unless
	a "make clean" command was issued first.  Now, if any .h files are
	updated, all dependent .c files are recompiled on the next "make".


Optional July average temperature field in soil parameter file

	Files affected:
	check_files.c, compute_treeline.c, initialize_atmos.c,
	initialize_global.c, read_soilparam.c, read_soilparam_arc.c, vicNl.c,
	vicNl.h, vicNl_def.h

	Description:
	In 4.0.4, when using the COMPUTE_TREELINE option, if the forcing files
	did not contain data for at least one July, the computed July average
	temperature would be undefined, resulting in incorrect location of the
	treeline.  In addition, runs covering different periods of the same
	forcing data would have different July average temperatures, causing the
	location of the treeline to change from run to run.  To allow more user
	control over the location of the treeline, added an optional field to
	the soil file, containing average July air temperature.  If the soil
	file contains 54 fields (or if the arc/info soil file list has 1 more
	line), this final entry is assumed to be the average July air
	temperature.  If COMPUTE_TREELINE is TRUE, this temperature is used in
	computing the treeline rather than calculating the average July
	temperature from the forcing data.  If COMPUTE_TREELINE is FALSE, this
	field is ignored.  If this field is not present in the soil file, vic
	behaves as in 4.0.4.


Bug Fixes:
----------

Spurious condensation at low temperatures

	Files affected:
	arno_evap.c

	Description:
	Changed logic of evap limit check to avoid creating spurious
	condensation.  In VIC 4.0.4 and earlier, whenever evaporation > (liquid
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
	In 4.0.4, runoff.c checked for the wrong bounds on baseflow, allowing
	baseflow to become negative when liquid soil moisture < residual
	moisture.  These bounds have been fixed in 4.0.5, as follows: baseflow
	is not allowed to exceed (liquid soil moisture - residual moisture); when
	baseflow < 0, baseflow is set to 0; when baseflow > 0 and the resulting
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
	In 4.0.4, check_state_file.c would increment the index of the forcing
	data array until it reached the record corresponding to the date
	stored in the state file.  This caused the simulation to start at the
	date at which the state file was saved rather than the start date
	specified in the global parameter file.  If the state file's date was
	earlier than the start date in the global parameter file, the index
	would be incremented until a segmentation fault occurred.  This has
	been fixed in 4.0.5 so that the start date in the global parameter
	file is always the start date of the simulation.  The date stored in
	the initial state file is ignored.

	NOTE: If you have been relying on the state file to dictate when your
	simulations start, this fix may require you to change your global
	parameter file so that STARTYEAR, STARTMONTH, etc. reflect the start
	date/time you want.


Negative incoming shortwave radiation at high latitudes

	Files affected:
	mtclim42_vic.c

	Description:
	In 4.0.4, when sub-daily shortwave radiation is estimated from daily
	min/max temperatures, negative values occasionally are calculated in
	coastal areas above the Arctic Circle in winter.  Now, if estimated sub-
	daily incident shortwave is negative, it is set to 0.0.


Undefined daily precipitation for deserts

	Files affected:
	mtclim42_vic.c

	Description:
	In 4.0.4, if a grid cell's annual precipitation (specified in the soil
	parameter file) is 0, then the adjusted daily precipitation calculated
	in mtclim42_vic.c ends up being undefined.  In 4.0.5 this has been
	fixed.

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


100% snow when air_temp = MAX_SNOW_TEMP

	Files affected:
	calc_rainonly.c

	Description:
	In 4.0.4, when air_temp = MAX_SNOW_TEMP, the portion of precipitation
	that is snow was set to 100%.  This has been fixed.  (found by Justin
	Sheffield at Princeton)


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


Forcing data validation

	Files affected:
	read_atmos_data.c

        Description:
	Replaced NF with global_param.dt in condition checking whether forcing
	file contains enough records to cover the time range of the simulation.
	(port from 4.1.0 beta)


Model state validation

	Files affected:
	read_initial_model_state.c

	Description:
	Added check to verify that the sum of the defined nodes equals the damping depth. (port from 4.1.0 beta)


Arc/info soil parameter files and ARNO soil parameters

	Files affected:
	read_soilparam_arc.c

	Description:
	VIC 4.0.4 does not handle ARNO soil parameters in arc/info soil files.
	VIC 4.0.5 now handles them in all soil files (assuming ARC_SOIL is
	TRUE).


Soil thermal node calculations

	Files affected:
	soil_conduction.c

	Description:
	set_node_parameters(): Modified to correct differences between
	calculations to determine maximum node moisture and node moisture, so
	that nodes on the boundary between soil layers are computed the same
	way for both. (port from 4.1.0 beta)
	distribute_node_moisture_properties(): Modified to check that node soil
	moisture is less than or equal to maximum node soil moisture, otherwise
	an error is printed to the screen and the model exits. (port from 4.1.0
	beta)


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
	calc_surf_energy_bal.c, frozen_soil.c, root_brent.c, snow_melt.c,
	vicNl.h

	Description:
	Instead of printing error messages, root_brent.c passes error messages
	to the functions that called it, leaving it to them to describe the
	specific error and its consequences.  In 4.0.4, root_brent's messages
	sometimes were wrong.


Statefile name change

	Files affected:
	open_state_file.c

	Description:
	Modified the statefile name to contain year, month, day rather than
	day, month, year.  This makes it consistent with the planned release
	of 4.1.0. (port from 4.1.0 beta)


Display current grid cell number for arc/info soil files

	Files affected:
	read_soilparam_arc.c

	Description:
	In 4.0.4, for arc/info-format soil files, the current grid cell is
	not displayed (for regular-format soil files, it is displayed).  Now
	the current grid cell number is always displayed.


Fluxes file output format (ascii)

	Files affected:
	write_data.c

	Description:
	In 4.0.4, some fields in the ascii-format fluxes output file were not
	separated by white space when values became large.  This has been
	fixed.


Binary file opening message

	Files affected:
	open_file.c

	Description:
	In 4.0.4, the message announcing the opening of a binary input file
	for reading was truncated.  This has been fixed.

Initial state reading changed to read initial state not final state

	Files affected:
	initialize_model_state.c

	Description:
	Changed read statement to read init_state not statefile.


--------------------------------------------------------------------------------

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

--------------------------------------------------------------------------------

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

--------------------------------------------------------------------------------

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

--------------------------------------------------------------------------------

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

--------------------------------------------------------------------------------

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
	     

--------------------------------------------------------------------------------

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

--------------------------------------------------------------------------------

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

--------------------------------------------------------------------------------

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
