# VIC Model Output Formatting

## Specifying Output Files and Variables

VIC allows the user to specify exactly which output files to create and which variables to store in each file. This way, users can save space by only writing those variables that are useful, and will be less likely to need to maintain a private version of the code to do this.

**Main points:**

1.  Output file names and contents can be specified in the [global parameter file](GlobalParam.md) (see below).
2.  If you do not specify file names and contents in the [global parameter file](GlobalParam.md), VIC will produce the same set of output files that it has produced in earlier versions, namely `fluxes_` and `snow_` files, plus `fdepth_` files if `FROZEN_SOIL` is TRUE. These files will have the same contents and format as in earlier versions.

**To specify file names and contents in the [global parameter file](GlobalParam.md):**

1.  Find the names of the desired variables in the [output variable list](../../OutputVarList.md)
2.  Decide how many output files you would like, what to name them, and which output variables will appear in each of these output files
3.  Add this information to the [global parameter file](GlobalParam.md) in the following format:

```
# Output File Contents
OUTFILE	_prefix_
OUTFREQ         _freq_          _VALUE_
COMPRESS        _compress_
OUTVAR	_varname_	[_format_	[_type_	[_multiplier_ [_aggtype_]]]]
OUTVAR	_varname_	[_format_	[_type_	[_multiplier_ [_aggtype_]]]]
OUTVAR	_varname_	[_format_	[_type_	[_multiplier_ [_aggtype_]]]]

OUTFILE	_prefix_
OUTVAR	_varname_	[_format_	[_type_	[_multiplier_ [_aggtype_]]]]
OUTVAR	_varname_	[_format_	[_type_	[_multiplier_ [_aggtype_]]]]
OUTVAR	_varname_	[_format_	[_type_	[_multiplier_ [_aggtype_]]]]
```
where

```
 _prefix_     = name of the output file, NOT including latitude
                and longitude
 _freq_       = Describes aggregation frequency for output stream. Valid
                options for frequency are:
                  NEVER     = never write to history file
                  NSTEPS    = write to history every _value_ steps
                  NSECONDS  = write to history every _value_ seconds
                  NMINUTES  = write to history every _value_ minutes
                  NHOURS    = write to history every _value_ hours
                  NDAYS     = write to history every _value_ days
                  NMONTHS   = write to history every _value_ months
                  NYEARS    = write to history every _value_ years
                  DATE      = write to history on the date: _value_
                  END       = write to history at the end of the simulation
 _value_      = integer describing the number of _freq_ intervals to pass
                before writing to the history file.
 _compress_   = gzip compression option.  TRUE, FALSE, or integer between 1-9.
 _varname_    = name of the variable (this must be one of the
                output variable names listed in vic_driver_shared_all.h.)
 _format_     = (for ascii output files) fprintf format string,
                e.g.
                  %.4f = floating point with 4 decimal places
                  %.7e = scientific notation w/ 7 decimal places
                  *    = use the default format for this variable

 _format_, _type_, _multiplier_, and _aggtype_ are optional.
 these.  If these are omitted, the default values will be used.

 _type_       = (for binary output files) data type code.
                Must be one of:
                  OUT_TYPE_DOUBLE = double-precision floating point
                  OUT_TYPE_FLOAT  = single-precision floating point
                  OUT_TYPE_INT    = integer
                  OUT_TYPE_USINT  = unsigned short integer
                  OUT_TYPE_SINT   = short integer
                  OUT_TYPE_CHAR   = char
                  *               = use the default type
 _multiplier_ = (for binary output files) factor to multiply
                the data by before writing, to increase precision.
                  *    = use the default multiplier for this variable
 _aggtype_    = Aggregation method to use for temporal aggregation. Valid
                options for aggtype are:
                  AGG_TYPE_DEFAULT = default aggregation type for variable
                  AGG_TYPE_AVG     = average over aggregation window
                  AGG_TYPE_BEG     = beginning of aggregation window
                  AGG_TYPE_END     = end of aggregation window
                  AGG_TYPE_MAX     = maximum in aggregation window
                  AGG_TYPE_MIN     = minimum in aggregation window
                  AGG_TYPE_SUM     = sum over aggregation window
```

Here's an example. To specify 2 output files, named `wbal` and `ebal`, and containing water balance and energy balance terms, respectively, you could do something like this:

```
OUTFILE	wbal
AGGFREQ         NDAYS           1
COMPRESS        FALSE
OUT_FORMAT      ASCII
OUTVAR	OUT_PREC        %.7g * * AGG_TYPE_AVG
OUTVAR	OUT_EVAP        %.7g * * AGG_TYPE_AVG
OUTVAR	OUT_RUNOFF      %.7g * * AGG_TYPE_AVG
OUTVAR	OUT_BASEFLOW    %.7g * * AGG_TYPE_AVG
OUTVAR	OUT_SWE         %.7g * * AGG_TYPE_END
OUTVAR	OUT_SOIL_MOIST  %.7g * * AGG_TYPE_AVG

OUTFILE	ebal
AGGFREQ         NHOURS           3
COMPRESS        TRUE
OUT_FORMAT      BINARY
OUTVAR	OUT_NET_SHORT
OUTVAR	OUT_NET_LONG
OUTVAR	OUT_LATENT
OUTVAR	OUT_SENSIBLE
OUTVAR	OUT_GRND_FLUX
OUTVAR	OUT_SNOW_FLUX
OUTVAR	OUT_ALBEDO
```

In the second file, none of the _format_, _type, _multiplier_, or _aggtype_ parameters were specified for any variables, VIC will use the default _format_, _type, _multiplier_, or _aggtype_ for the variables.

For example, to specify scientific notation with 10 significant digits, you could do the following:

```
OUTVAR	OUT_ALBEDO	%.9e
```

Note that even if you only want to specify the format, you must supply a value in the type and multiplier columns as well. This can be `*` to indicate the default value. Similarly, if you only want to specify the type (e.g. as a double), you would need to do something like:

```
OUTVAR	OUT_ALBEDO	*	OUT_TYPE_DOUBLE
```

**Date variables:**

For typical output files, the date is always written at the beginning of each record. This will consist of the following columns:

year month day seconds

For daily output timestep, "seconds" is not written.

If `OUT_FORMAT` is `BINARY`, these will all be written as type int (OUT_TYPE_INT).

**Multiple-valued variables:**

Since variables like SOIL_MOIST have 1 value per soil layer, these variables will be written to multiple columns in the output file, one column per soil layer. Other multiple-valued variables are treated similarly.

**Snow band output:**

To specify writing the values of variables in each snow band, append "BAND" to the variable name (this only works for some variables - see the list in `vic_driver_shared_all.h`). If you specify these variables, the value of the variable in each band will be written, one band per column. For example, for a cell having 2 snow bands:

```
OUTVAR	OUT_SWE_BAND
OUTVAR	OUT_ALBEDO_BAND
```

will result in an output file containing:

```
year month day (seconds) swe[0] swe[1] albedo[0] albedo[1]
```

## Specifying Output Time Step

VIC can now aggregate the output variables to a user-defined output interval, via the `OUTFREQ` setting in the [global parameter file](GlobalParam.md). When  `OUTFREQ` is set, it describes aggregation frequency for an output stream. Valid options for frequency are: NEVER, NSTEPS, NSECONDS, NMINUTES, NHOURS, NDAYS, NMONTHS, NYEARS, DATE, END. Count may be a positive integer or a string with date format YYYY-MM-DD[-SSSSS] in the case of DATE. Default `frequency` is `NDAYS`. Default `count` is 1.

## Output File Headers

Now VIC provides descriptive headers into its output files. VIC will insert a short header into its output files, describing the simulation and variables included in the file. See our [Best Practices](../../best_practices.md) page for further information on how to add descriptive metadata to VIC output.

For ascii files, the output header has the following format:

```
# SIMULATION: (OUTFILE prefix)
# MODEL_VERSION: (Version String)
VARNAME    VARNAME   VARNAME   ...
```
where

- SIMULATION: OUTFILE prefix is taken from the global parameter file
- MODEL_VERSION: VIC Version String (e.g VIC version 5.0.0)

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
