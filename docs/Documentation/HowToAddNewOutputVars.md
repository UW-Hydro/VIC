# How to Add New Output Variables

As discussed [here](OutputFormatting.md), VIC allows users to select which output variables to write to their output files. Variables that can be selected are listed [here](OutputVarList.md) as well as in the VIC source code file **vicNl_def.h** (listed as `OUT_*`). However, not all variables that VIC simulates are available for output. And of course if you are changing VIC's physics, any new variables you add are also not available for output (yet).

You can enable VIC to output new variables by making a few changes to the code, as follows:

## 1\. Make sure VIC sends the variable to put_data()

If VIC doesn't already compute the variable you need, or computes it but doesn't store it in one of the structures that is accessible to put_data(), you'll have to modify VIC to do this. The structures that put_data() has access to are: dmy, atmos_data, soil_con, cell_data, veg_con, veg_var, energy, snow_data, and all_vars.

## 2\. Define the variable's name in vicNl_def.h.

The list of output variables begins with a definition of the number of output variables in the list:

```C
/***** Output Variable Types *****/
#define N_OUTVAR_TYPES 100
                        ^
                        | number of output variables (can be larger than the
                          actual number of variables)
```

This is followed by a list of output variable names and ID codes:

```C
// Water Balance Terms - state variables
#define OUT_ROOTMOIST        0  /* root zone soil moisture  [mm] */
#define OUT_SMFROZFRAC       1  /* fraction of soil moisture (by mass) that is ice, for each soil layer */
             ^               ^
           name            ID code

(etc)
```

Similarly, for input meteorological variables, we have:

```C
/***** Forcing Variable Types *****/
#define N_FORCING_TYPES 23
```
followed by:

```
#define AIR_TEMP   0 /* air temperature per time step [C] (ALMA_INPUT: [K]) */
#define ALBEDO     1 /* surface albedo [fraction] */
```

To add a variable to either list, copy an existing entry and modify it to have your variable's name and ID code. If you're creating an output variable, the name should begin with "OUT_", similar to the other output variables. Variable names must be no longer than 19 characters.

NOTE: The variable's ID code MUST BE UNIQUE!!! The safest thing to do is to create your new entry at the end of the list. Your new variable's ID code should be the previous variable's ID code + 1\. If you would like to insert your new entry elsewhere in the list (e.g. to group it with similar variables), then you will need to update all the ID codes in the list so that each code = the previous variable's ID code + 1.

Next, if necessary, increase the number of variables defined at the top of the list (N_OUTVAR_TYPES for output variables, or N_FORCING_TYPES for input variables) so that it is >= to the number of variables in the list.

## 3\. For output variables, populate the variable's metadata in output_list_utils.c.

In the function create_output_list(), there is a series of calls to strcpy, one per variable:

```C
strcpy(out_data[OUT_ROOTMOIST].varname,"OUT_ROOTMOIST");
                     ^                       ^
                array index             variable name
```
You must add an entry for your variable. The safest thing to do is to copy an existing entry and modify it. Make sure that the array index and variable name both match the name you entered in vicNl_def.h EXACTLY.

Next, if your variable consists of multiple values per grid cell (for example, per-layer soil moisture OUT_SOIL_MOIST has options.Nlayer values, while OUT_RUNOFF has only 1 value) you must define the number of elements for the variable. There is a series of lines defining the number of elements, as:

```C
out_data[OUT_SMLIQFRAC].nelem = options.Nlayer;
```

To add an entry, copy an existing entry and modify it.

Next, define your variable's aggregation method. This refers to the method used to aggregate up from smaller time steps to larger time steps, e.g. from hourly to daily. By default, variables are aggregated via averaging the original values over the output time step; if this is acceptable for your variable, then you do not need to specify anything here. Otherwise, you need to specify a method, e.g.:

```C
out_data[OUT_ROOTMOIST].aggtype = AGG_TYPE_END;
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
