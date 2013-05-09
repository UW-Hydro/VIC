/*
 * SUMMARY:      snow.h - header file for DHSVM snow routines
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:              nijssen@u.washington.edu
 * ORIG-DATE:    29-Aug-1996 at 16:03:11
 * LAST-MOD: Thu Feb 22 14:28:53 2001 by Keith Cherkauer <cherkaue@u.washington.edu>
 * $Id$
 * DESCRIPTION:  header file for DHSVM snow routines
 * DESCRIP-END.
 * FUNCTIONS:    
 * COMMENTS:     
 * MODIFICATIONS:     
  2008-Feb-17 Moved parameters related to snow albedo and trace snow
	      from user_def.h to here.					TJB
  2008-Feb-17 Moved parameters related to snow densification from
	      snow_utility.c to here.					TJB
  2008-Feb-17 Added parameters for new snow densification algorithm
	      based on Lundberg and Pomeroy (1998);
	      removed NEW_SNOW_DENSITY.					KMA via TJB
  2008-Apr-21 Re-inserted NEW_SNOW_DENSITY, for backwards compatibility
	      with previous snow density algorithm.			TJB
 */


#ifndef SNOW_H
#define SNOW_H

#include <stdarg.h>

/* water holding capacity of snow as a fraction of snow-water-equivalent */ 
#define LIQUID_WATER_CAPACITY		0.035

/* multiplier to calculate the amount of available snow interception as
   a function of LAI (m) */
/* #define LAI_SNOW_MULTIPLIER		0.0005 */
#define LAI_SNOW_MULTIPLIER		0.0005

/* the amount of snow on the canopy that can only be melted off. (m) */
#define MIN_INTERCEPTION_STORAGE	0.005

/* maximum depth of the surface layer in water equivalent (m)
   [default 0.125] */
#define MAX_SURFACE_SWE			0.125

/* density of new fallen snow [50] */
#define NEW_SNOW_DENSITY        50.

/* Density limit used in calculation of destructive metamorphism */
#define SNDENS_DMLIMIT			100. // (kg/m^3)

/* Constants in snow density computation */
#define SNDENS_ETA0	(3.6e6) /* viscosity of snow at T = 0C and density = 0
				   used in calculation of true viscosity (Ns/m2) */
#define SNDENS_C1	0.04
#define SNDENS_C2	(2.778e-6)
#define SNDENS_C5	0.08	/* constant used in snow viscosity calculation,
				   taken from SNTHRM.89 (/C)*/
#define SNDENS_C6	0.021	/* constant used in snow viscosity calculation,
				   taken from SNTHRM.89 (kg/m3) */
#define SNDENS_F	0.6	/* internal compaction rate coefficient */

/* Minimum SWQ for which the snowpack energy balance is computed 
   independent of the soil surface temperature */
/* #define MIN_SWQ_EB_THRES 0.0060 */
#define MIN_SWQ_EB_THRES		0.0010
 
/* Attenuation coefficients for shortwave in a snowpack.  Values and
   equation taken from Patterson and Hamblin, 1988 */
#define SNOW_A1 0.7
#define SNOW_A2 0.3  
#define SNOW_L1 6    // (1/m)
#define SNOW_L2 20   // (1/m)
 
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

#endif 
