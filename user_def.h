/**********************************************************************
  This header file contains model parameters that can be modified by
  the user to control model performance.  When this file is modified 
  the model needs to be recompiled for the changes to take effect.

  NOTE (to those who add or remove parameters from VIC):
  Any additions or removals of parameters in this file must also
  be made in display_current_settings.c.
**********************************************************************/

/***** If TRUE include all model messages to stdout, and stderr *****/
#define VERBOSE TRUE

/***** If TRUE limit output data to runoff and baseflow for optimization *****/
#define OPTIMIZE FALSE

/***** If TRUE all energy balance calculations are iterated to minimize
       the total column (air, canopy, snow and ground) error.  Otherwise
       no iteration is used and the model estimates the new fluxes
       based on those from the previous time step, results should
       be similar, however, the model will report energy balance errors *****/
#define CLOSE_ENERGY FALSE

/***** If TRUE include all debugging code - debugging options still
       have to be activated to get extra output.  When set to FALSE
       all debugging if-then statements are removed from the compiled 
       code *****/
#define LINK_DEBUG FALSE

/***** If TRUE output will be in LDAS binary format, which is a single
       file with limited variables, most of which are truncated to
       conserve disk space *****/
#define LDAS_OUTPUT FALSE

/***** If TRUE VIC uses a system of linear equations defined in global.h
       to estimate the maximum unfrozen water content equation.  This 
       significantly reduces the run time with frozen soil, but may
       introduce new errors (STILL UNDER TESTING) *****/
#define QUICK_FS FALSE
#define QUICK_FS_TEMPS 7

/***** If TRUE VIC binary includes the lake/wetland algorithm.  Lakes
       and wetlands do not need to be defined but code will be 
       accessed. *****/
#define LAKE_MODEL FALSE

/***** If TRUE VIC uses the linear interpolation of the logarithm of the
       matric potential from the two surrounding layers to estimate the 
       soil moisture drainage from each layer (Boone and Wetzel, 1996).
       This should improve the soil moisture drainage predicted by the
       low resolution solution computed by VIC. *****/
#define LOW_RES_MOIST FALSE

/***** If TRUE VIC code to save the model state is included in the 
       compiled code.  If STATEYEAR, STATEMONTH and STATEDAY are
       defined in the global control file the model state will
       be written to a file. *****/
#define SAVE_STATE FALSE

/***** If TRUE VIC does not rewind the vegetation, state, and snow
       band files before read data for each cell.  This saves time
       but requires that all grid cells are listed in the same
       order as the soil parameter file *****/
#define NO_REWIND TRUE

/***** If TRUE VIC reads the model forcing files, and creates the full
       internal forcing dataset (longwave, shortwave, humidity, etc.)
       which is then written to a series of gridded output files for
       later use.  Gridded forcing files are written to the RESULTS
       directory defined in the global control file, and are binary
       or ASCII based on the BINARY_OUTPUT flag. *****/
#define OUTPUT_FORCE FALSE

/***** Compute the treeline elevation.  If set to TRUE this flag will 
       force the VIC model to compute the elevation of the tree line, 
       based on elevation at which the average annual July air temperature
       is at or below 10C.  All snowbands above this evelation are then 
       assumed to be above the treeline, and vegetation types with 
       overstory are removed from the snow band average variables. *****/
#define COMPUTE_TREELINE TRUE

/***** If TRUE VIC computes the mean, standard deviation, and sum
       and finds the minimum and maximum values of the forcing 
       variables for each grid cell and outputs the results to 
       stdout.  These values are meant to be a quick check for
       obvious errors in the forcing data, more thorough checks 
       of the data should be conducted outside of the model by the 
       user. *****/
#define OUTPUT_FORCE_STATS FALSE

/***** If TRUE VIC uses a uniform distribution function to simulate
       the spatial distribution of soil frost, if FALSE VIC assumes
       that the entire grid cell is frozen uniformly *****/
#define SPATIAL_FROST FALSE
#define FROST_SUBAREAS 10

/***** If TRUE VIC uses a uniform distribution to simulate the partial
       coverage of the surface by a thin snowpack.  Coverage is 
       assumed to be uniform after snowfall until the pack begins to 
       melt. SiB uses 0.076, from Rosemount I want 0.155cm depth ~ 0.028mm swq *****/
#define SPATIAL_SNOW FALSE


/***** Define maximum array sizes for model source code *****/
#define MAX_VEG        12      /* maximum number of vegetation types per 
				  cell */
#define MAX_LAYERS     3       /* maximum number of soil moisture layers */
#define MAX_NODES      18      /* maximum number of soil thermal nodes */
#define MAX_BANDS      10      /* maximum number of snow bands */
#define MAX_FRONTS     3       /* maximum number of freezing and thawing 
				  front depths to store */
#define MAX_LAKE_NODES 20      /* maximum number of lake thermal nodes */

/***** Number of iterations to use in solving the surface energy balance.
       The original VIC model uses only 1 iteration for speed.  Increasing
       the number of iterations improves precision, and is recommended
       for single point comparisons with frozen soils *****/
#define MAXIT_FE        25
 
/***** Coefficient multiplied by the LAI to determine the amount of
       water that can be storde in the canopy *****/
#define LAI_WATER_FACTOR 0.2

/***** Longwave correction factor, used to correct estimated incoming
       longwave radiation (use 1, unless measured longwave available for
       calibration) *****/
#define LWAVE_COR	1.

/***** Snow albedo curve parameters.  Defaults are from Bras p263.
       Should not be changed except for serious problems with snow melt *****/
#define NEW_SNOW_ALB		0.85
#define SNOW_ALB_ACCUM_A	0.94
#define SNOW_ALB_ACCUM_B	0.58
#define SNOW_ALB_THAW_A		0.82
#define SNOW_ALB_THAW_B		0.46

/***** Defines the minimum amount of new snow (mm) which will reset the
       snowpack albedo to new snow *****/
#define TraceSnow 0.03
