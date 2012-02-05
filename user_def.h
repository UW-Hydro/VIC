/**********************************************************************
  This header file contains model parameters that can be modified by
  the user to control model performance.  When this file is modified 
  the model needs to be recompiled for the changes to take effect.

  NOTE (to those who add or remove parameters from VIC):
  Any additions or removals of parameters in this file must also
  be made in display_current_settings.c.

  $Id$

  Modifications:
  2006-Sep-23 Implemented flexible output configuration; removed the
              OPTIMIZE and LDAS_OUTPUT options.				TJB
  2006-Nov-07 Removed LAKE_MODEL option.				TJB
  2008-Feb-17 Moved constants related to snow albedo and trace snow to
	      snow.h.							TJB
  2009-Jan-12 Removed COMPUTE_TREELINE option (moved into options
	      struct).							TJB
  2011-Jan-04 Added MAX_ZWTVMOIST for storing array of soil moisture vs
	      water table position.					TJB
  2012-Jan-16 Removed LINK_DEBUG code					BN

**********************************************************************/

/***** If TRUE include all model messages to stdout, and stderr *****/
#define VERBOSE TRUE

/***** If TRUE all energy balance calculations are iterated to minimize
       the total column (air, canopy, snow and ground) error.  Otherwise
       no iteration is used and the model estimates the new fluxes
       based on those from the previous time step, results should
       be similar, however, the model will report energy balance errors *****/
#define CLOSE_ENERGY FALSE

/***** If TRUE VIC uses a system of linear equations defined in global.h
       to estimate the maximum unfrozen water content equation.  This 
       significantly reduces the run time with frozen soil, but may
       introduce new errors (STILL UNDER TESTING, ALSO NEEDS DEBUGGING) *****/
#define QUICK_FS FALSE
#define QUICK_FS_TEMPS 7

/***** If TRUE VIC uses the linear interpolation of the logarithm of the
       matric potential from the two surrounding layers to estimate the 
       soil moisture drainage from each layer (Boone and Wetzel, 1996).
       This should improve the soil moisture drainage predicted by the
       low resolution solution computed by VIC. *****/
#define LOW_RES_MOIST FALSE


/***** If TRUE VIC does not rewind the vegetation, state, and snow
       band files before read data for each cell.  This saves time
       but requires that all grid cells are listed in the same
       order as the soil parameter file *****/
#define NO_REWIND FALSE

/***** If TRUE VIC reads the model forcing files, and creates the full
       internal forcing dataset (longwave, shortwave, humidity, etc.)
       which is then written to a series of gridded output files for
       later use.  Gridded forcing files are written to the RESULTS
       directory defined in the global control file, and are binary
       or ASCII based on the BINARY_OUTPUT flag. *****/
#define OUTPUT_FORCE FALSE

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

/***** If TRUE VIC allows for excess ground ice, i.e. an expanded porosity
       to account for an initial volumetric ice fraction larger than
       soil porosity.  The porosity decreases as the excess ice melts.  
       Once porosity reaches the soil porosity (1-bulk density/soil density), 
       it does not change.  Initial volumetric ice fraction must be 
       defined in the soil file for each soil layer. *****/
#define EXCESS_ICE FALSE

/***** Define maximum array sizes for model source code *****/
#define MAX_VEG        12      /* maximum number of vegetation types per 
				  cell */
#define MAX_LAYERS     3       /* maximum number of soil moisture layers */
#define MAX_NODES      50      /* maximum number of soil thermal nodes */
#define MAX_BANDS      10      /* maximum number of snow bands */
#define MAX_FRONTS     3       /* maximum number of freezing and thawing 
				  front depths to store */
#define MAX_LAKE_NODES 20      /* maximum number of lake thermal nodes */
#define MAX_ZWTVMOIST  11      /* maximum number of points in water table vs moisture curve for each soil layer; should include points at lower and upper boundaries of the layer */

/***** Number of iterations to use in solving the surface energy balance.
       The original VIC model uses only 1 iteration for speed.  Increasing
       the number of iterations improves precision, and is recommended
       for single point comparisons with frozen soils *****/
#define MAXIT_FE        25
 
/***** Coefficient multiplied by the LAI to determine the amount of
       water that can be stored in the canopy *****/
#define LAI_WATER_FACTOR 0.2

/***** Longwave correction factor, used to correct estimated incoming
       longwave radiation (use 1, unless measured longwave available for
       calibration) *****/
#define LWAVE_COR	1.
