#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "vicNl.h"

static char vcid[] = "$Id$";

double calc_atmos_energy_bal(double  InOverSensible,
			     double  InUnderSensible,
			     double  LatentHeatOver,
			     double  LatentHeatUnder,
			     double  LatentHeatSubOver,
			     double  LatentHeatSubUnder,
			     double  Lv,
			     double  NetLongOver,
			     double  NetLongUnder,
			     double  NetShortOver,
			     double  NetShortUnder,
			     double  Ra,
			     double  Tair, 
			     double  atmos_density,
			     double  vp,
			     double  vpd,
			     double *Error,
			     double *LatentHeat,
			     double *LatentHeatSub,
			     double *NetLongAtmos,
			     double *NetShortAtmos,
			     double *SensibleHeat, 
			     double *VPcanopy,
			     double *VPDcanopy) {
/************************************************************************
  calc_atmos_energy_bal.c        Keith Cherkauer       February 6, 2001

  This routine was written to iteratively solve the atmospheric energy
  balance equation to estimate the canopy air temperature.

  Basic concept for this was taken from:
  Sellers et al., J. Clim., v.9, April 1996, pp. 676-705.
  Dickinsen, BATS manual, NCAR Tech. Note (NCAR/TN-387+STR), August 1993.

  Modifications:
  04-Jun-04 Added more descriptive error message to beginning of screen
	    dump in error_print_atmos_moist_bal.		TJB
  21-Sep-04 Added ErrorString to store error messages from
	    root_brent.						TJB
2007-Apr-06 Modified to handle grid cell errors by returning to the
            main subroutine, rather than ending the simulation. GCT/KAC
 2007-Aug-31 Checked root_brent return value against -998 rather than -9998.    JCA
************************************************************************/

  double AtmosLatent;
  double F; // canopy closure fraction, not currently used by VIC
  double InLatent;
  double InSensible;
  double NetRadiation;
  double T_lower;
  double T_upper;
  double Tcanopy;
  double VP_lower;
  double VP_upper;
  double gamma;
  char ErrorString[MAXSTRING];
  
  F = 1;

  // compute incoming sensible heat
  InSensible = InOverSensible + InUnderSensible;

  // compute net radiation
  (*NetLongAtmos ) = (F * NetLongOver + (1. - F) * NetLongUnder);

  (*NetShortAtmos) = (NetShortOver + NetShortUnder);

  NetRadiation = (*NetShortAtmos + *NetLongAtmos); 

  // compute total latent heat flux
  (*LatentHeat)    = (LatentHeatOver + LatentHeatUnder);

  (*LatentHeatSub) = (LatentHeatSubOver + LatentHeatSubUnder);

  InLatent = (*LatentHeat) + (*LatentHeatSub);

  /******************************
    Find Canopy Air Temperature
  ******************************/

  /* set initial bounds for root brent **/
  T_lower = (Tair) - CANOPY_DT;
  T_upper = (Tair) + CANOPY_DT;

  // iterate for canopy air temperature
  Tcanopy = root_brent(T_lower, T_upper, ErrorString, func_atmos_energy_bal, 
		       (*LatentHeat) + (*LatentHeatSub), 
		       NetRadiation, Ra, Tair, atmos_density, InSensible, 
		       SensibleHeat);

  if ( Tcanopy <= -998 ) {
    // handle error flag from root brent
    (*Error) = error_calc_atmos_energy_bal(Tcanopy, (*LatentHeat) 
					   + (*LatentHeatSub), 
					   NetRadiation, Ra, Tair, 
					   atmos_density, InSensible, 
					   SensibleHeat, ErrorString);
    return ( ERROR );
  }
  // compute varaibles based on final temperature
  (*Error) = solve_atmos_energy_bal(Tcanopy, (*LatentHeat) + (*LatentHeatSub), 
				    NetRadiation, Ra, Tair, atmos_density, 
				    InSensible, SensibleHeat);

  /*****************************
    Find Canopy Vapor Pressure
  *****************************/

  /* set initial bounds for root brent **/
/*   VP_lower = (vp) - CANOPY_VP; */
/*   VP_upper = (vp) + CANOPY_VP; */

/*   gamma = svp_slope(Tair); */

  // iterate for canopy vapor pressure
/*   (*VPcanopy) = root_brent(VP_lower, VP_upper, ErrorString, func_atmos_moist_bal,  */
/* 			   InLatent, Lv, Ra, atmos_density, gamma, vp,  */
/* 			   &AtmosLatent); */

/*   if ( (*VPcanopy) <= -998 )  */
    // handle error flag from root brent
/*     (*Error) = error_calc_atmos_moist_bal((*VPcanopy), InLatent,  */
/* 					  Lv, Ra, atmos_density, gamma, vp,  */
/* 					  &AtmosLatent, ErrorString); */
  
  // compute varaibles based on final vapor pressure
/*   (*Error) = solve_atmos_moist_bal( (*VPcanopy), InLatent,  */
/* 				    Lv, Ra, atmos_density, gamma, vp,  */
/* 				    &AtmosLatent); */

  // compute vapor pressure deficit in canopy
/*   (*VPDcanopy) = vpd + (*VPcanopy - vp); */


  // Bi-pass above computations
  // (*VPcanopy) = vp;


  return(Tcanopy);

}

double solve_atmos_energy_bal(double Tcanopy, ...) {

  va_list ap;

  double error;

  va_start(ap, Tcanopy);
  error = func_atmos_energy_bal(Tcanopy, ap);
  va_end(ap);

  return error;

}

double error_calc_atmos_energy_bal(double Tcanopy, ...) {

  va_list ap;

  double error;

  va_start(ap, Tcanopy);
  error = error_print_atmos_energy_bal(Tcanopy, ap);
  va_end(ap);

  return error;

}

double error_print_atmos_energy_bal(double Tcanopy, va_list ap) {

  double  LatentHeat;
  double  NetRadiation;
  double  Ra;
  double  Tair;
  double  atmos_density;
  double  InSensible;

  double *SensibleHeat;
  char *ErrorString;
 
  // extract variables from va_arg
  LatentHeat    = (double)  va_arg(ap, double);
  NetRadiation  = (double)  va_arg(ap, double);
  Ra            = (double)  va_arg(ap, double);
  Tair          = (double)  va_arg(ap, double);
  atmos_density = (double)  va_arg(ap, double);
  InSensible    = (double)  va_arg(ap, double);

  SensibleHeat  = (double *)va_arg(ap, double *);
  ErrorString   = (char *)va_arg(ap, char *);

  // print variable values
  fprintf(stderr, "%s", ErrorString);
  fprintf(stderr, "ERROR: calc_atmos_energy_bal failed to converge to a solution in root_brent.  Variable values will be dumped to the screen, check for invalid values.\n");
  fprintf(stderr, "LatentHeat = %f\n",  LatentHeat);
  fprintf(stderr, "NetRadiation = %f\n",  NetRadiation);
  fprintf(stderr, "Ra = %f\n",  Ra);
  fprintf(stderr, "Tair = %f\n",  Tair);
  fprintf(stderr, "atmos_density = %f\n",  atmos_density);
  fprintf(stderr, "InSensible = %f\n",  InSensible);

  fprintf(stderr, "*SensibleHeat = %f\n", *SensibleHeat);
 
  fprintf(stderr, "Finished writing calc_atmos_energy_bal variables.\nTry increasing CANOPY_DT to get model to complete cell.\nThen check output for instabilities.\n");

  return( ERROR );
    
}

double solve_atmos_moist_bal(double VPcanopy, ...) {

  va_list ap;

  double error;

  va_start(ap, VPcanopy);
  error = func_atmos_moist_bal(VPcanopy, ap);
  va_end(ap);

  return error;

}

double error_calc_atmos_moist_bal(double VPcanopy, ...) {

  va_list ap;

  double error;

  va_start(ap, VPcanopy);
  error = error_print_atmos_moist_bal(VPcanopy, ap);
  va_end(ap);

  return error;

}

double error_print_atmos_moist_bal(double VPcanopy, va_list ap) {

 
  double  InLatent;
  double  Lv;
  double  Ra;
  double  atmos_density;
  double  gamma;
  double  vp;
  double *AtmosLatent;
  char *ErrorString;

  // extract variables from va_arg
  InLatent      = (double)  va_arg(ap, double);
  Lv            = (double)  va_arg(ap, double);
  Ra            = (double)  va_arg(ap, double);
  atmos_density = (double)  va_arg(ap, double);
  gamma         = (double)  va_arg(ap, double);
  vp            = (double)  va_arg(ap, double);
  AtmosLatent   = (double *)va_arg(ap, double *);
  ErrorString   = (char *)  va_arg(ap, char *);

  // print variable values
  fprintf(stderr, "%s", ErrorString);
  fprintf(stderr, "InLatent = %f\n",  InLatent);
  fprintf(stderr, "Lv = %f\n",  Lv);
  fprintf(stderr, "Ra = %f\n",  Ra);
  fprintf(stderr, "atmos_density = %f\n",  atmos_density);
  fprintf(stderr, "gamma = %f\n",  gamma);
  fprintf(stderr, "vp = %f\n", vp);
  fprintf(stderr, "AtmosLatent = %f\n", *AtmosLatent);
 
  vicerror("Finished writing calc_atmos_moist_bal variables.\nTry increasing CANOPY_VP to get model to complete cell.\nThen check output for instabilities.");

  return(0.0);
    
}

