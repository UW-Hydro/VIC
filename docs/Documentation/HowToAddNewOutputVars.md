# How to Add New Output Variables

VIC allows users to select which output variables to write to their output files. Variables that can be selected are listed [here](OutputVarList.md) as well as in the VIC source code file `vic_def.h` (listed as `OUT_*`). However, not all variables that VIC simulates are available for output. And of course if you are changing VIC's physics, any new variables you add are also not available for output (yet).

You can enable VIC to output new variables by making a few changes to the code, as follows:

## 1\. Make sure VIC sends the variable to put_data()

If VIC doesn't already compute the variable you need, or computes it but doesn't store it in one of the structures that is accessible to put_data(), you'll have to modify VIC to do this. The structures that put_data() has access to are: dmy, atmos_data, soil_con, cell_data, veg_con, veg_var, energy, snow_data, and all_vars.

## 2\. Define the variable's name in `vic_def.h`.

The list of output variables is a C `enum`.  To add a new variable, begin by appending the new variable to the ***next to*** last position in the `enum`. See the example of adding `OUT_NEW_VAR_NAME` below:

```C
/******************************************************************************
 * @brief   Output Variable Types
 *****************************************************************************/
enum
{
    // Water Balance Terms - state variables
    OUT_ASAT,             /**< Saturated Area Fraction */
    ...
    OUT_CSLOW,            /**< Carbon density in slow pool [g C/m2] */
    OUT_NEW_VAR_NAME      /**< Description of new output variable [units] */
    // Last value of enum - DO NOT ADD ANYTHING BELOW THIS LINE!!
    // used as a loop counter and must be >= the largest value in this enum
    N_OUTVAR_TYPES        /**< used as a loop counter*/
};
```

Similarly, for input meteorological variables, we have:

```C
/******************************************************************************
 * @brief   Forcing Variable Types
 *****************************************************************************/
enum
{
    AIR_TEMP,    /**< air temperature per time step [C] (ALMA_INPUT: [K]) */
    ...
    WIND_N,      /**< meridional component of wind speed [m/s] */
    NEW_FORCE_VAR, /**< Description of new forcing variable [units] */
    SKIP,        /**< place holder for unused data columns */
    // Last value of enum - DO NOT ADD ANYTHING BELOW THIS LINE!!
    // used as a loop counter and must be >= the largest value in this enum
    N_FORCING_TYPES  /**< Number of forcing types */
};
```

## 3\. For output variables, populate the variable's metadata in output_list_utils.c.

In the function create_output_list(), there is a series of calls to strcpy, one per variable:

```C
strcpy(out_data[OUT_NEW_VAR_NAME].varname,"OUT_NEW_VAR_NAME");
                     ^                       ^
                array index             variable name
```
You must add an entry for your variable. The safest thing to do is to copy an existing entry and modify it. Make sure that the array index and variable name both match the name you entered in `vic_def.h` EXACTLY.

Next, if your variable consists of multiple values per grid cell (for example, per-layer soil moisture OUT_SOIL_MOIST has options.Nlayer values, while OUT_RUNOFF has only 1 value) you must define the number of elements for the variable. There is a series of lines defining the number of elements, as:

```C
out_data[OUT_NEW_VAR_NAME].nelem = options.Nlayer;
```

To add an entry, copy an existing entry and modify it.

Next, define your variable's aggregation method. This refers to the method used to aggregate up from smaller time steps to larger time steps, e.g. from hourly to daily. By default, variables are aggregated via averaging the original values over the output time step; if this is acceptable for your variable, then you do not need to specify anything here. Otherwise, you need to specify a method, e.g.:

```C
out_data[OUT_NEW_VAR_NAME].aggtype = AGG_TYPE_END;
```

Once again, to add an entry, copy an existing entry and modify it. Possible aggregation methods are:

```
AGG_TYPE_AVG : Aggregated value = average of the values over the interval
                                  (default; used for energy balance terms)
AGG_TYPE_END : Aggregated value = final value over the interval (used for
                                  moisture storage terms)
AGG_TYPE_SUM : Aggregated value = sum of the values over the interval (used
                                  for moisture fluxes)
AGG_TYPE_BEG : Aggregated value = first value over the interval
AGG_TYPE_MIN : Aggregated value = minimum of the values over the interval
AGG_TYPE_MAX : Aggregated value = maximum of the values over the interval
```

## 4\. For output variables, add logic to put_data.c to set the variable in the out_data[] array

Assuming that at this point, your variable is computed somewhere in VIC and stored in one of the data structures that put_data() has access to, you now need to assign this to the appropriate part of the out_data structure. In put_data(), there is a loop over elevation bands and veg tiles. The contribution of each band/tile combination is added to the running total (weighted by the band/tile's area fraction) in the out_data structure. For example, for single-element variables:

```C
/** record canopy interception **/
if ( veg < veg_con[0].vegetat_type_num )
  out_data[OUT_WDEW].data[0] += veg_var[veg][band].Wdew
    * Cv * mu * AreaFract[band] * TreeAdjustFactor[band];
```

For multiple-element variables, you must loop over the variable's index, e.g.:

```C
for(index = 0; index < MAX_FRONTS; index++) {
  if(energy[veg][band].fdepth[index] != MISSING)
    out_data[OUT_FDEPTH].data[index] += energy[veg][band].fdepth[index]
      * Cv * 100\. * AreaFract[band] * TreeAdjustFactor[band];
  if(energy[veg][band].tdepth[index] != MISSING)
    out_data[OUT_TDEPTH].data[index] += energy[veg][band].tdepth[index]
      * Cv * 100\. * AreaFract[band] * TreeAdjustFactor[band];
}
```

## 5\. Add the relevant information to the global parameter file.

Input variables should be specified in the "Forcing Files" section.

Output variables should be specified in the "Output Files" section, unless using the default output file format.

Note: if you simply want to add an output variable to the default output files, you can either modify the default output file format in set_output_defaults.c so that it includes your variable (and not add any information to your global parameter file), or you can paste the contents of the appropriate `output.*.template` file into your global parameter file, and then add your variable to that description. See the documentation on the flexible output configuration for more details on specifying output file contents in the global parameter file.

## 6\. If you have added new state variables to VIC, you must add them to state files

State variables are any type of storage term that must be "remembered" from one step to the next to preserve the model's internal state. Examples include snow water equivalent, soil moisture, soil temperature, etc. If these are not stored in the state file and read from the state file upon restart, the restarted simulation will not be identical to the same point in time of a continuous simulation, i.e. the missing states (e.g., snow water equivalent) will be reset to 0.

To do this, you must add writing of the new variable to write_model_state() and add reading of the new variable to read_initial_model_state().

In addition, you must provide default initial values for the new variable in the case of no initial state file. These initial values should be assigned either in initialize_model_state() or one of the functions it calls, e.g. initialize_soil(), initialize_snow(), etc.
