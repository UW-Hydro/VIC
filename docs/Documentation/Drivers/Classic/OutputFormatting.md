# VIC Model Output Formatting

## How to control the contents of VIC output files

These instructions apply to VIC model versions 4.0.6 and later. In earlier versions of VIC, the user has little control over the contents of VIC output files without changing the VIC source code.

## Specifying Output Files and Variables

In earlier versions of VIC, the set of output files and their contents were hard-coded. Users inevitably had to modify the code to produce the type of output they wanted, and had to merge their changes into any new versions of the code as they were released.

VIC (4.0.6 and later) now allows the user to specify exactly which output files to create and which variables to store in each file. This way, users can save space by only writing those variables that are useful, and will be less likely to need to maintain a private version of the code to do this.

**Main points:**

1.  Output file names and contents can be specified in the [global parameter file](GlobalParam.md) (see below).
2.  If you do not specify file names and contents in the [global parameter file](GlobalParam.md), VIC will produce the same set of output files that it has produced in earlier versions, namely "fluxes" and "snow" files, plus "fdepth" files if FROZEN_SOIL is TRUE and "snowband" files if PRT_SNOW_BAND is TRUE. These files will have the same contents and format as in earlier versions.
3.  The OPTIMIZE and LDAS_OUTPUT options have been removed. These output configurations can be selected with the proper set of instructions in the [global parameter file](GlobalParam.md). (see the `output.*.template` files included in this distribution for more information.)
4.  If you do specify the file names and contents in the [global parameter file](GlobalParam.md), PRT_SNOW_BAND will have no effect.

**To specify file names and contents in the [global parameter file](GlobalParam.md):**

1.  Find the names of the desired variables in the [output variable list](../../OutputVarList.md)
2.  Decide how many output files you would like, what to name them, and which output variables will appear in each of these output files
3.  Add this information to the [global parameter file](GlobalParam.md) in the following format:

```
# Output File Contents
OUTFILE	_prefix_
OUTVAR	_varname_	[_format_	_type_	_multiplier_]
OUTVAR	_varname_	[_format_	_type_	_multiplier_]
OUTVAR	_varname_	[_format_	_type_	_multiplier_]

OUTFILE	_prefix_
OUTVAR	_varname_	[_format_	_type_	_multiplier_]
OUTVAR	_varname_	[_format_	_type_	_multiplier_]
OUTVAR	_varname_	[_format_	_type_	_multiplier_]
```
where

_prefix_ = name of the output file, NOT including latitude and longitude

_varname_ = name of the variable (this must be one of the output variable names listed in `vic_driver_shared.h`.)

_format_, _type_, and _multiplier_ are optional.  For a given variable,
you can specify either NONE of these, or ALL of these.  If these
are omitted, the default values will be used.

_format_ = (for ascii output files) `fprintf` format string, e.g.
  - `%.4f` = floating point with 4 decimal places
  - `%.7e` = scientific notation w/ 7 decimal places
  - `*` = use the default format for this variable

_type_ = (for binary output files) data type code. Must be one of:
  - `OUT_TYPE_DOUBLE` = double-precision floating point
  - `OUT_TYPE_FLOAT` = single-precision floating point
  - `OUT_TYPE_INT` = integer
  - `OUT_TYPE_USINT` = unsigned short integer
  - `OUT_TYPE_SINT` = short integer
  - `OUT_TYPE_CHAR` = char
  - `*` = use the default type

_multiplier_ = (for binary output files) factor to multiply the data by before writing, to increase precision.
  - `*` = use the default multiplier for this variable

Here's an example. To specify 2 output files, named "wbal" and "ebal", and containing water balance and energy balance terms, respectively, you could do something like this:

```
OUTFILE	wbal
OUTVAR	OUT_PREC
OUTVAR	OUT_EVAP
OUTVAR	OUT_RUNOFF
OUTVAR	OUT_BASEFLOW
OUTVAR	OUT_SWE
OUTVAR	OUT_SOIL_MOIST

OUTFILE	ebal
OUTVAR	OUT_NET_SHORT
OUTVAR	OUT_NET_LONG
OUTVAR	OUT_LATENT
OUTVAR	OUT_SENSIBLE
OUTVAR	OUT_GRND_FLUX
OUTVAR	OUT_SNOW_FLUX
OUTVAR	OUT_ALBEDO
```

Since no format, type, or multiplier were specified for any variables, VIC will use the default format, type, and multiplier for the variables.

If you wanted scientific notation with 10 significant digits for ALBEDO, you could do the following:

```OUTVAR	OUT_ALBEDO	%.9e	*	*```

Note that even if you only want to specify the format, you must supply a value in the type and multiplier columns as well. This can be `*` to indicate the default value. Similarly, if you only want to specify the type (e.g. as a double), you would need to do something like:

```OUTVAR	OUT_ALBEDO	*	OUT_TYPE_DOUBLE	*```

**Date variables:**

For typical output files, the date is always written at the beginning of each record. This will consist of the following columns:

year month day hour

For daily output timestep, "hour" is not written.

If BINARY_OUTPUT is TRUE, these will all be written as type int (OUT_TYPE_INT).

**Multiple-valued variables:**

Since variables like SOIL_MOIST have 1 value per soil layer, these variables will be written to multiple columns in the output file, one column per soil layer. Other multiple-valued variables are treated similarly.

**Snow band output:**

To specify writing the values of variables in each snow band, append "BAND" to the variable name (this only works for some variables - see the list in vic_driver_shared.h). If you specify these variables, the value of the variable in each band will be written, one band per column. For example, for a cell having 2 snow bands:

```
OUTVAR	OUT_SWE_BAND
OUTVAR	OUT_ALBEDO_BAND
```

will result in an output file containing:

```year month day (hour) swe[0] swe[1] albedo[0] albedo[1]```

## Specifying Units

The user now has some control over the units of the input and output variables. The standard VIC units for moisture fluxes are total mm over the output time interval, and degrees C for temperatures. However, other land surface schemes and circulation or climate models tend to use mm/s for moisture fluxes and degrees K for temperatures.

Now there are options in the [global parameter file](GlobalParam.md) that allow you to specify whether to use traditional VIC units or the mm/s and K convention for input or output variables. The option names are "ALMA_INPUT" and "ALMA_OUTPUT", named after the [ALMA convention](http://www.lmd.jussieu.fr/~polcher/ALMA/) used in the PILPS-2e experiment.

**ALMA INPUT:**

VIC now accepts the following new ALMA-compliant input forcings in addition to the forcings that it already accepts:

```
SNOWF     snowfall rate (kg/m^2s)
RAINF     rainfall rate (kg/m^2s)
CRAINF    convective rainfall rate (kg/m^2s)
LSRAINF   large scale rainfall rate (kg/m^2s)
QAIR      specific humidity (kg/kg)
WIND_E    zonal wind speed (m/s)
WIND_N    meridional wind speed (m/s)
TAIR      air temperature per time step (K)
PSURF     atmospheric pressure (Pa)
```

When giving VIC ALMA-compliant input files, you must be sure to use the names given above in the forcing section of your [global parameter file](GlobalParam.md).

Instead of the existing PREC (precipitation per timestep in mm), you can now specify SNOWF and RAINF (snowfall and rainfall rates, both in mm/s). VIC will simply add these two quantities together, multiply by the forcing interval, and treat their sum the same way it treats PREC.

An alternative to supplying RAINF is to supply CRAINF (convective rainfall rate, mm/s) and LSRAINF (large-scale rainfall rate, mm/s). VIC will add these two quantities together to get RAINF.

Instead of the existing WIND, alternatively you can specify WIND_E and WIND_N (zonal and meridional wind speed, m/s). VIC will simply compute `WIND = sqrt(WIND_E**2+WIND_N**2)`.

TAIR has units of K, while the existing AIR_TEMP is in C. Similarly, PSURF is in Pa, while PRESSURE is in kPa. VIC will convert these to AIR_TEMP and PRESSURE after reading them in.

More information is available on ALMA forcing variables at: [http://www.lmd.jussieu.fr/~polcher/ALMA/convention_input_3.html](http://www.lmd.jussieu.fr/~polcher/ALMA/convention_input_3.html)

**ALMA OUTPUT:**

If the user sets ALMA_OUTPUT=TRUE in the global parameter file, then VIC will convert its output variables to ALMA-compliant forms. The majority of the changes are changes of units. Moisture fluxes are changed from VIC's standard (mm accumulated over the time step) to the average flux rate (mm/s). Temperatures are converted from C to K.

More information on ALMA output variables is available at: [http://www.lmd.jussieu.fr/~polcher/ALMA/convention_output_3.html](http://www.lmd.jussieu.fr/~polcher/ALMA/convention_output_3.html).

In addition, several more variables have been added to the list of available output variables. See `vic_driver_shared.h` for the complete list of available output variables.

## Specifying Output Time Step

VIC can now aggregate the output variables to a user-defined output interval, via the OUT_STEP setting in the [global parameter file](GlobalParam.md). Currently, the largest output interval allowed is 24 hours, so this option is only useful for simulations running at sub-daily time steps.

## Optional Output File Headers

Now VIC provides an option to insert descriptive headers into its output files, via the PRT_HEADER option in the [global parameter file](GlobalParam.md). If this is set to TRUE, VIC will insert a short header into its output files, describing the time step, start date/time, variables and units included in the file.

For ascii files, the output header has the following format:

```
# NRECS: (nrecs)
# DT: (dt)
# STARTDATE: yyyy-mm-dd hh:00:00
# ALMA_OUTPUT: (0 or 1)
# NVARS: (Nvars)
# VARNAME    VARNAME   VARNAME   ...
```
where

- nrecs = Number of records in the file
- dt = Time step length in hours
- start date = Date and time of first record of file
- ALMA_OUTPUT = Indicates units of the variables; 0 = standard VIC units; 1 = ALMA units
- Nvars = Number of variables in the file, including date fields

For binary files, the output header has the following format:

```
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
```
