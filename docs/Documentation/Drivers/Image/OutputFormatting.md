# VIC Model Output Formatting - Image Driver

## Specifying Output Files and Variables

VIC allows the user to specify exactly which output files to create and which variables to store in each file. This way, users can save space by only writing those variables that are useful, and will be less likely to need to maintain a private version of the code to do this.

**Main points:**

1.  Output file names and contents can be specified in the [global parameter file](GlobalParam.md) (see below).
2.  If you do not specify file names and contents in the [global parameter file](GlobalParam.md), VIC will produce the same set of output files that it has produced in earlier versions, namely "fluxes" and "snow" files, plus "fdepth" files if FROZEN_SOIL is TRUE and "snowband" files if PRT_SNOW_BAND is TRUE. These files will all be in netCDF format.
3.  The OPTIMIZE and LDAS_OUTPUT options have been removed. These output configurations can be selected with the proper set of instructions in the [global parameter file](GlobalParam.md). (see the `output.*.template` files included in this distribution for more information.)
4.  If you do specify the file names and contents in the [global parameter file](GlobalParam.md), PRT_SNOW_BAND will have no effect.

**To specify file names and contents in the [global parameter file](GlobalParam.md):**

1.  Find the names of the desired variables in the [output variable list](../../OutputVarList.md)
2.  Decide how many output files you would like, what to name them, and which output variables will appear in each of these output files
3.  Add this information to the [global parameter file](GlobalParam.md) in the following format:

```
# Output File Contents
OUTFILE	_prefix_
OUTVAR	_varname_	[_aggtype_]
OUTVAR	_varname_	[_aggtype_]
OUTVAR	_varname_	[_aggtype_]

OUTFILE	_prefix_
OUTVAR	_varname_	[_aggtype_]
OUTVAR	_varname_	[_aggtype_]
OUTVAR	_varname_	[_aggtype_]
```
where

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; _varname_ = name of the variable (this must be one of the output variable names listed in `vic_driver_shared.h`.)

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;  _aggtype_ = Aggregation method to use for temporal aggregation. Valid options for aggtype are: <br />
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; `AGG_TYPE_DEFAULT` = default aggregation type for variable <br />
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; `AGG_TYPE_AVG` = average over aggregation window <br />
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; `AGG_TYPE_BEG` = beginning of aggregation window <br />
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; `AGG_TYPE_END` = end of aggregation window <br />
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; `AGG_TYPE_MAX` = maximum in aggregation window <br />
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; `AGG_TYPE_MIN` = minimum in aggregation window <br />
  &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; `AGG_TYPE_SUM` = sum over aggregation window <br />

NOTE: currently, the output of the whole domain and run period is saved to a single netCDF file; in later version, we hope to let the user specify the length of simulation results to be saved in one netCDF file, so that splitting a long-period simulation results to multiple output files will be easy.

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

Since no _aggtype_ were specified for any variables, VIC will use the default aggregate type for the variables.


**Output file format:**

For each OUTVAR specified in the [global parameter file](GlobalParam.md), one netCDF file will be output. In each output netCDF file, each specified output variable has the following attributes:

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; _FillValue_ = value for inactive cells <br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; _long_name_ = variable long name <br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; _standard_name_ = variable standard name <br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; _units_ = variable unit <br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; _description_ = variable description <br />
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; _cell_methods_ = variable aggregation method <br />

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
variables:                                          
        double time(time) ;                         
                time:standard_name = "time" ;       
                time:units = "days since 0001-01-01 00:00:00" ;
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



**Snow band output:**

To specify writing the values of variables in each snow band, append "BAND" to the variable name (this only works for some variables - see the list in vic_driver_shared.h). If you specify these variables, the value of the variable in each band will be written, in a similar data structure as shown above (but there will be an extra "snow_band" dimension for each variable). For example, for a cell having 2 snow bands, you can specify the following in the [global parameter file](GlobalParam.md):

```
OUTVAR	OUT_SWE_BAND
OUTVAR	OUT_ALBEDO_BAND
```

will result in an output netCDF file containing:

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
