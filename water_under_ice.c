/*
 * SUMMARY:      water_under_ice.c 
 * USAGE:        
 *
 * AUTHOR:       Laura Bowling
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:     8-Oct-1996 at 08:50:06
 * LAST-MOD: Mon Nov 12 16:38:32 2001 by Keith Cherkauer <cherkaue@u.washington.edu>
 * DESCRIPTION:  
 * DESCRIP-END.
 * FUNCTIONS:   
 * COMMENTS:     
 */

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

#define MAX_ITER 50

//static char vcid[] = "$Id$";

/*****************************************************************************
  Function name: water_under_ice()

  Purpose      : 

  Required     :


  Comments     :

  Modifications:
  15-Jun-04 Initialize mixmax to 0.						TJB
  04-Oct-04 Merged with Laura Bowling's updated lake model code.		TJB
  2006-Nov-07 Removed LAKE_MODEL option.					TJB
  2007-Apr-03 Modified to handle grid cell errors by returning to the
	      main subroutine, rather than ending the simulation.		KAC via GCT
  2007-Nov-06 Replaced lake.fraci with lake.areai.  Added workaround for
	      non-convergence of temperatures.					LCB via TJB
  2009-Dec-11 Replaced "assert" statements with "if" statements.		TJB
*****************************************************************************/
int water_under_ice(int     freezeflag, 
		    double  sw_ice,
		    double  wind,
		    double *Ti,
		    double *water_density,
		    double  lat,
		    int     numnod,
		    double  dz,
		    double  surfdz,
		    double  Tcutoff,
		    double *qw,
		    double *surface,
		    double *deltaH,
		    double *water_cp,
		    int     mixdepth,
		    double  hice,
		    double  sdepth,
		    double  dt,
		    double *energy_out_bottom)
	      
{
  double Tnew[MAX_LAKE_NODES];
  int k;
  int iterations;
  double jouleold;
  double joulenew;
  double error;
  double de[MAX_LAKE_NODES]; 
  double epsilon = 0.0001;
  double qw_init, qw_mean, qw_final;
  double sw_underice_visible, sw_underice_nir;
  
  iterations = 0;

  for(k=0; k<numnod; k++)
    Tnew[k] = Ti[k];

  // compute the eddy diffusivity 
  eddy(freezeflag, wind, Ti, water_density, de, lat, numnod, dz, surfdz);
  
  // estimate the flux out of the water
  qw_init =  0.57*(Ti[0]-Tcutoff)/(surfdz/2.);
  *qw = qw_init;
  qw_mean = -999.;

  energycalc(Ti, &jouleold, numnod,dz, surfdz, surface, water_cp, water_density);

  while((fabs(qw_mean - *qw) > epsilon) && iterations < MAX_ITER) {
    
    if(iterations == 0)
      *qw = qw_init;
    else
      *qw = qw_mean;

    // compute shortwave that transmitted through the lake ice 
    sw_underice_visible = a1*sw_ice*exp(-1.*(lamisw*hice+lamssw*sdepth));
    sw_underice_nir = a2*sw_ice*exp(-1.*(lamilw*hice+lamslw*sdepth));
	
    /* --------------------------------------------------------------------
     * Calculate the lake temperatures at different levels for the
     * new timestep.
     * -------------------------------------------------------------------- */
	
    temp_area (sw_underice_visible, sw_underice_nir, -1.*(*qw) , Ti, Tnew,
	       water_density, de, dt, surface, numnod, 
	       dz, surfdz, &joulenew, water_cp, energy_out_bottom);

    // recompute storage of heat in the lake
    *deltaH = (joulenew - jouleold)/(surface[0]*dt*SECPHOUR);
      
    /* --------------------------------------------------------------------
     * Do the convective mixing of the lake water.
     * -------------------------------------------------------------------- */

    tracer_mixer (Tnew, &mixdepth, freezeflag,
    		  surface,numnod, dz, surfdz, water_cp);
    
    qw_final = 0.57*(Tnew[0]-Tcutoff)/(surfdz/2.);
    
    qw_mean = (qw_final + *qw)/2.;  
    
    iterations += 1;

  }

  if(fabs(qw_mean - *qw) <= epsilon) {
    // Temperature reached convergence
    for(k=0; k<numnod; k++)
      Ti[k] = Tnew[k];
    *qw = qw_mean;
    return(0);
  }
  else {
//    fprintf(stderr, "Lake temps under ice failed to converge; temporary work around used.\n");
    *qw = 0.0;
    for(k=0; k<numnod; k++)
      Ti[k] = Tcutoff;
    energycalc(Ti, &joulenew, numnod,dz, surfdz, surface, water_cp, water_density);
    *deltaH = (joulenew - jouleold)/(surface[0]*dt*SECPHOUR);
    return(0);
    // return (ERROR);
  }

}
