# VIC Model Output Formatting

The VIC Image Driver writes output files using the [netCDF](http://www.unidata.ucar.edu/software/netcdf/) file format. Each output stream may include 1 or more variables and may aggregate model history at a user specified output frequency.

## Specifying Output Files and Variables

VIC allows the user to specify exactly which output files to create and which variables to store in each file.

**Main points:**

1.  Output file names and contents can be specified in the [global parameter file](GlobalParam.md) (see below).
2.  If you do not specify file names and contents in the [global parameter file](GlobalParam.md), VIC will produce the same set of output files that it has produced in earlier versions, namely `fluxes_` and `snow_` files, plus `fdepth_` files if `FROZEN_SOIL` is TRUE. See the [documentation](DefaultOutputs.md) on the default output files for more information.

**To specify file names and contents in the [global parameter file](GlobalParam.md):**

1.  Find the names of the desired variables in the [output variable list](../../OutputVarList.md)
2.  Decide how many output files you would like, what to name them, and which output variables will appear in each of these output files
3.  Add this information to the [global parameter file](GlobalParam.md) in the following format:

```
# Output File Contents
OUTFILE	_prefix_
OUTFREQ         _freq_          _VALUE_
HISTFREQ        _freq_          _VALUE_
COMPRESS        _compress_
OUT_FORMAT      _nc_format_
OUTVAR	_varname_	[_format_  [_type_ [_multiplier_ [_aggtype_]]]]
OUTVAR	_varname_	[_format_  [_type_ [_multiplier_ [_aggtype_]]]]
OUTVAR	_varname_	[_format_  [_type_ [_multiplier_ [_aggtype_]]]]

OUTFILE	_prefix_
OUTFREQ         _freq_          _VALUE_
OUTVAR	_varname_	[_format_  [_type_ [_multiplier_ [_aggtype_]]]]
OUTVAR	_varname_	[_format_  [_type_ [_multiplier_ [_aggtype_]]]]
OUTVAR	_varname_	[_format_  [_type_ [_multiplier_ [_aggtype_]]]]
```

where

```
 _prefix_     = name of the output file, NOT including the date stamp or the suffix
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
 _compress_   = netCDF gzip compression option.  TRUE, FALSE, or integer between 1-9.
 _nc_format_  = netCDF format. NETCDF3_CLASSIC, NETCDF3_64BIT_OFFSET,
                NETCDF4_CLASSIC, or NETCDF4
 _varname_    = name of the variable (this must be one of the
                output variable names listed in vic_driver_shared_all.h.)

_format_     = not used in image driver, replace with *

_type_, and _multiplier_, and _aggtype_ are optional.
If these are omitted, the default values will be used.

 _type_       = data type code. Must be one of:
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
HISTFREQ        NYEARS          1
COMPRESS        FALSE
OUT_FORMAT      NETCDF3_CLASSIC
OUTVAR	OUT_PREC        * * AGG_TYPE_AVG
OUTVAR	OUT_EVAP        * * AGG_TYPE_AVG
OUTVAR	OUT_RUNOFF      * * AGG_TYPE_AVG
OUTVAR	OUT_BASEFLOW    * * AGG_TYPE_AVG
OUTVAR	OUT_SWE         * * AGG_TYPE_END
OUTVAR	OUT_SOIL_MOIST  * * AGG_TYPE_AVG

OUTFILE	ebal
AGGFREQ         NHOURS           3
COMPRESS        TRUE
OUT_FORMAT      NETCDF4
OUTVAR	OUT_NET_SHORT
OUTVAR	OUT_NET_LONG
OUTVAR	OUT_LATENT
OUTVAR	OUT_SENSIBLE
OUTVAR	OUT_GRND_FLUX
OUTVAR	OUT_SNOW_FLUX
OUTVAR	OUT_ALBEDO
```

In the second file, none of the _type, _multiplier_, or _aggtype_ parameters were specified for any variables, VIC will use the default _type, _multiplier_, or _aggtype_ for the variables.

**Multiple-valued variables:**

Since variables like SOIL_MOIST have 1 value per soil layer, these variables will be written along a fourth netCDF dimension. Other multiple-valued variables are treated similarly.

**Snow band output:**

To specify writing the values of variables in each snow band, append "BAND" to the variable name (this only works for some variables - see the list in `vic_driver_shared_all.h`). If you specify these variables, the values of the variable in each band will be written in an additional netCDF dimension.

```
OUTVAR	OUT_SWE_BAND
OUTVAR	OUT_ALBEDO_BAND
```

will result in an output file containing:

```
float OUT_SWE_BAND(time, snow_band, lat, lon) ;                                     
        OUT_SWE_BAND:_FillValue = 9.96921e+36f ;                                    
        OUT_SWE_BAND:long_name = "swe_band" ;                                       
        OUT_SWE_BAND:standard_name = "lwe_thickness_of_snow" ;                      
        OUT_SWE_BAND:units = "mm" ;                                                 
        OUT_SWE_BAND:description = "snow water equivalent in snow pack" ;
        OUT_SWE_BAND:cell_methods = "time: mean" ;
float OUT_ALBEDO_BAND(time, snow_band, lat, lon) ;                                  
        OUT_ALBEDO_BAND:_FillValue = 9.96921e+36f ;                                 
        OUT_ALBEDO_BAND:long_name = "albedo_band" ;                                 
        OUT_ALBEDO_BAND:standard_name = "surface_albedo" ;                          
        OUT_ALBEDO_BAND:units = "1" ;                                               
        OUT_ALBEDO_BAND:description = "albedo" ;                                    
        OUT_ALBEDO_BAND:cell_methods = "time: mean" ;

```

## Specifying Output Time Step

VIC can now aggregate the output variables to a user-defined output interval, via the `OUTFREQ` setting in the [global parameter file](GlobalParam.md). When  `OUTFREQ` is set, it describes aggregation frequency for an output stream. Valid options for frequency are: NEVER, NSTEPS, NSECONDS, NMINUTES, NHOURS, NDAYS, NMONTHS, NYEARS, DATE, END. Count may be a positive integer or a string with date format YYYY-MM-DD[-SSSSS] in the case of DATE. Default `frequency` is `NDAYS`. Default `count` is 1.

The number of output records per output file is controlled by the `HISTFREQ` option. Valid options for the`HISTFREQ` option are NEVER, NSTEPS, NSECONDS, NMINUTES, NHOURS, NDAYS, NMONTHS, NYEARS, DATE, END. Count may be a positive integer or a string with date format YYYY-MM-DD[-SSSSS] in the case of DATE. Default `frequency` is `NDAYS`. Default `count` is END.

**Output file format:**

For each `OUTVAR` specified in the [global parameter file](GlobalParam.md), one netCDF file will be output. In each output netCDF file, each specified output variable has the following attributes:

  - _FillValue_ = value for inactive cells
  - _long_name_ = variable long name
  - _standard_name_ = variable standard name
  - _units_ = variable unit
  - _description_ = variable description
  - _cell_methods_ = variable aggregation method

Here is an example of the structure and content of an output netCDF file:

```
dimensions:                                         
        snow_band = 1 ;                             
        front = 3 ;                                 
        frost_area = 1 ;                            
        nlayer = 3 ;                                
        lon = 5 ;                                   
        lat = 4 ;                                   
        node = 3 ;                                  
        root_zone = 3 ;                             
        veg_class = 12 ;                            
        time = UNLIMITED ; // (5 currently)
        nv = 2 ;
variables:                                          
        double time(time) ;                         
                time:standard_name = "time" ;       
                time:units = "days since 0001-01-01 00:00:00" ;
                time:calendar = "proleptic_gregorian" ;
        double time_bnds(time, nv) ;
                time_bnds:standard_name = "time_bounds" ;
                time_bnds:units = "days since 0001-01-01 00:00:00" ;
                time:calendar = "proleptic_gregorian" ;
        double lon(lon) ;                                      
                lon:long_name = "longitude" ;                  
                lon:units = "degrees_east" ;                   
                lon:standard_name = "longitude" ;              
        double lat(lat) ;                                      
                lat:long_name = "latitude" ;                   
                lat:units = "degrees_north" ;                  
                lat:standard_name = "latitude" ;               
        float OUT_PREC(time, lat, lon) ;                       
                OUT_PREC:_FillValue = 9.96921e+36f ;           
                OUT_PREC:long_name = "prec" ;                  
                OUT_PREC:standard_name = "precipitation_amount" ;
                OUT_PREC:units = "mm" ;                          
                OUT_PREC:description = "incoming precipitation" ;
                OUT_PREC:cell_methods = "time: sum" ;            
        float OUT_EVAP(time, lat, lon) ;                         
                OUT_EVAP:_FillValue = 9.96921e+36f ;             
                OUT_EVAP:long_name = "evap" ;                    
                OUT_EVAP:standard_name = "water_evaporation_flux_net" ;
                OUT_EVAP:units = "mm" ;                                
                OUT_EVAP:description = "total net evaporation" ;       
                OUT_EVAP:cell_methods = "time: sum" ;    
        float OUT_RUNOFF(time, lat, lon) ;                             
                OUT_RUNOFF:_FillValue = 9.96921e+36f ;                 
                OUT_RUNOFF:long_name = "runoff" ;                      
                OUT_RUNOFF:standard_name = "runoff_amount" ;           
                OUT_RUNOFF:units = "mm" ;                              
                OUT_RUNOFF:description = "surface runoff" ;            
                OUT_RUNOFF:cell_methods = "time: sum" ;                
        float OUT_BASEFLOW(time, lat, lon) ;                           
                OUT_BASEFLOW:_FillValue = 9.96921e+36f ;               
                OUT_BASEFLOW:long_name = "baseflow" ;                  
                OUT_BASEFLOW:standard_name = "baseflow_amount" ;       
                OUT_BASEFLOW:units = "mm" ;                            
                OUT_BASEFLOW:description = "baseflow out of the bottom layer" ;
                OUT_BASEFLOW:cell_methods = "time: sum" ;
        float OUT_SOIL_LIQ(time, nlayer, lat, lon) ;                                    
               OUT_SOIL_LIQ:_FillValue = 9.96921e+36f ;                                
               OUT_SOIL_LIQ:long_name = "soil_liq" ;                                   
               OUT_SOIL_LIQ:standard_name = "soil_moisture_liquid_depth" ;             
               OUT_SOIL_LIQ:units = "mm" ;                                             
               OUT_SOIL_LIQ:description = "soil liquid moisture content for each soil layer" ;
               OUT_SOIL_LIQ:cell_methods = "time: mean" ;
```
