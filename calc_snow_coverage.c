#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

double calc_snow_coverage(int    *store_snow,
			  double  depth_full_snow_cover,
			  double  old_coverage, 
			  double  swq,
			  double  old_swq,
			  double  depth,
			  double  old_depth,
			  double  melt,
			  double *max_swq, 
			  double *snowfall, 
			  double *store_swq,
			  double *swq_slope,
			  double *store_coverage) 
/**********************************************************************
  calc_snow_coverage.c      Keith Cherkauer          November 1, 2000

  This routine computes the current fraction of the vegetation band
  that is covered with snow.  The snow distribution is assumed to 
  be uniform with a slope based on the value of depth_full_snow_cover
  set in user_def.h.  The original value was based on field observations
  from the University of Minnesota's Rosemount Agricultural Experiment
  station (see dissertation by Keith Cherkauer, 2001).

  030701 Modified to accept minimum snow depth for full coverage as
         a passed variable instead of a globally defined value.   KAC

*************************************************************************/
{

  double old_max_swq;
  double coverage;

  /*************************************
    New snow falls on partial snowpack
  *************************************/

  if ( snowfall[0] > 0 ) {

    /*************************
      Continued accumulation
    *************************/

    coverage = 1;
    
    if ( *store_snow ) {

      // store coverage fraction before it is buried
      if ( *store_swq == 0 && old_coverage < 1 ) 
	*store_coverage = old_coverage;
      else if ( *store_swq == 0 ) *store_coverage = 1;

      /* store snow falling over partial snowpack */
      *store_swq += swq - old_swq;
      
      if ( depth >= depth_full_snow_cover / 2. ) {
	/* snow accumulation deep enough to remove memory of previous
	   melt distribution */
	*store_snow     = FALSE;
	*store_swq      = 0;
	*swq_slope      = 0;
	*store_coverage = 1;
      }
      
    }
    else if ( old_coverage < 1 ) {

      /* store new snow fall over old distribution */
      *store_snow = TRUE;
      *store_swq  = swq - old_swq;

    }

  }
  else if ( melt > 0 ) {

    /***********************************************
      Snowpack begins to melt or continues melting 
    ***********************************************/

    if ( *store_swq > 0 && swq < old_swq ) {
      
      /* Melt thin snowfall off previous distribution */
      *store_swq += swq - old_swq;

      if ( *store_swq <= 0 ) {

	/* Snowpack cover has melted - clear storage */
	*store_swq = 0;
	// restore buried cover fraction
	old_coverage    = *store_coverage;
	*store_coverage = 1;

      }

    }
    
    if ( *store_swq == 0 ) {

      /* compute the current coverage fraction */
      
      if ( (*swq_slope) == 0 ) {

	/* compute new distribution function */
	if ( old_depth > depth_full_snow_cover / 2. ) 
	  *swq_slope = -depth_full_snow_cover;
	else 
	  *swq_slope = -2. * old_depth;
	*max_swq   = -(*swq_slope);
	*store_snow = TRUE;

      }

      /* if currently raining, new swq may be higher than previous 
	 maximum swq even if melt occurs.  reset maximum swq if this
	 occurs. */
      
      old_max_swq = *max_swq;
      *max_swq = 2. * depth;

      if ( *max_swq < old_max_swq || old_max_swq == 0 ) {

	/* melt has occured, reduce coverage fraction */
	coverage = -(*max_swq) / (*swq_slope);
	coverage = ( coverage > 1 ) ? 1 : coverage;

      }
      else 
	/* rain or sublimation has increased swq, 
	   coverage fraction does not change */
	coverage = old_coverage;

/** NOTE if (swq - old_swq) > 0, but there was no snowfall, there was rain so the coverage fraction should not change */

    }
    else 
      coverage = old_coverage;
  }
  else 
    /* no change to snowpack */
    coverage = old_coverage;

  return(coverage);
}
