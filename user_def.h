/**********************************************************************
  This header file contains model parameters that can be modified by
  the user to control model performance.  When this file is modified 
  the model needs to be recompiled for the changes to take effect.
**********************************************************************/

/*#define EMULATE_XU*/	/* uses step function baseflow curve, remove
				after ark-red tests */

/***** Number of layers to use in energy balance and frozen soil models.
       Water balance model runs with only two layers, while the full
       energy model requires at least 3 layers *****/
#define MAXlayer	3

/***** Number of iterations to use in solving the surface energy balance.
       The original VIC model uses only 1 iteration for speed.  Increasing
       the number of iterations improves precision, and is recommended
       for single point comparisons with frozen soils *****/
#define MAXIT_FE        25
 
/***** Coefficient multiplied by the LAI to determine the amount of
       water that can be storde in the canopy *****/
#define LAI_WATER_FACTOR 0.2


/***** Offset of shortwave measurement from model hour, used for 
       computation of longwave radiation from hourly shortwave
       measurements *****/
#define SOLARTIMEOFFSET	0.0

/***** Longwave correction factor, used to correct estimated incoming
       longwave radiation (use 1, unless measured longwave available for
       calibration) *****/
#define LWAVE_COR	1.08

/***** Residual moisture content of soil column (water content mm/mm) *****/
#define MOIST_RESID 0.055

/***** Snow albedo curve parameters.  Defaults are from Bras p263.
       Should not be changed except for serious problems with snow melt *****/
#define NEW_SNOW_ALB	0.85		/** 0.85 **/
#define SNOW_ALB_ACCUM_A	0.94	/** 0.94 **/
#define SNOW_ALB_ACCUM_B	0.58	/** 0.58 **/
#define SNOW_ALB_THAW_A	0.82		/** 0.82 **/
#define SNOW_ALB_THAW_B	0.46		/** 0.46 **/
