/*
 * SUMMARY:      StabilityCorrection.c - Calculate the stability correction
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen and Pascal Storck
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu, pstorck@u.washington.edu
 * ORIG-DATE:    Apr-1996
 * LAST-MOD: Fri Mar  2 18:42:20 2001 by Keith Cherkauer <cherkaue@u.washington.edu>
 * DESCRIPTION:  Calculate the stability correction for exchange of sensible
 *               heat between the surface and the atmosphere 
 * DESCRIP-END.
 * FUNCTIONS:    StabilityCorrection()
 * COMMENTS:     
 */

#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

static char vcid[] = "$Id$";

/*****************************************************************************
  Function name: StabilityCorrection()

  Purpose      : Calculate atmospheric stability correction for non-neutral
                 conditions

  Required     :
    double Z          - Reference height (m)
    double d          - Displacement height (m)
    double TSurf      - Surface temperature (C)
    double Tair       - Air temperature (C)
    double Wind       - Wind speed (m/s)
    double Z0         - Roughness length (m)

  Returns      :
    double Correction - Multiplier for aerodynamic resistance

  Modifies     : None
    
  Comments     :
*****************************************************************************/
double StabilityCorrection(double Z, double d, double TSurf, double Tair, 
			  double Wind, double Z0)
{
  double Correction;		/* Correction to aerodynamic resistance */
  double Ri;                     /* Richardson's Number */
  double RiCr = 0.2;             /* Critical Richardson's Number */
  double RiLimit;                /* Upper limit for Richardson's Number */

  Correction = 1.0;

    /* Calculate the effect of the atmospheric stability using a Richardson 
       Number approach */ 

  if (TSurf != Tair) {

    /* Non-neutral conditions */
    
    Ri = G * (Tair - TSurf) * (Z - d)/
      (((Tair+273.15)+(TSurf+273.15))/2.0 * Wind * Wind);
    
    RiLimit = (Tair+273.15)/
      (((Tair+273.15)+(TSurf+273.15))/2.0 * (log((Z-d)/Z0) + 5));
    
    if (Ri > RiLimit)
      Ri = RiLimit;
    
    if (Ri > 0.0) 
      Correction = (1 - Ri/RiCr) * (1 - Ri/RiCr);
      
    else { 
      if (Ri < -0.5)
        Ri = -0.5;

      Correction = sqrt(1 - 16 * Ri);
    }
  }

  return Correction;

  /*  double Eta;*/			/* intermediate product */
    
  /* calculate the effect of atmospheric stability on aerodynamic resistance 
     using the method of Choudhury et al., Agric. For. Met., 37, 75-88, 
     1986 */
  
  /* Eq. A4, Choudhury et al [1986] */
  
  /*    Eta = 5.0 * (Z - d) * G * (TSurf - Tair)/
        ((Tair + 273.15) * Wind*Wind);
        
        if (TSurf < Tair) { */
  
  /* stable conditions */
  
  /* If Eta smaller or equal than -1, the correction increases again
     because of the quadratic form of the equation.  However, the
     correction should monnotonically decrease for increasing
     differences between the surface and air temperature.  Therefore a
     lower bound of -.8 is imposed on Eta (The lower bound of Eta is
     not put at -1, because this would result in no sensible heat
     exchange at all at the soil surface, which is unrealistic.  A lower
     bound of -.8 results in a correction of 25, i.e. an increase in the
     resistance of 25 times) */
  
  /*    if (Eta < -0.8)
        Eta = -0.8;
        Correction = (1 + Eta)*(1 + Eta);
        }
        
        else */
  
  /* unstable conditions */
  
  /*      Correction = pow((double) (1 + Eta), (double) 0.75);
          }*/
  
  /*   return Correction; */
}

