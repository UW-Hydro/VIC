/*
 * SUMMARY:      CalcAerodynamic.c - Calculate the aerodynamic resistances
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen and Pascal Storck
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu, pstorck@u.washington.edu
 * ORIG-DATE:    Thu Mar 27 18:00:10 1997
 * LAST-MOD: Thu Mar  8 13:24:10 2001 by Keith Cherkauer <cherkaue@u.washington.edu>
 * DESCRIPTION:  Calculate the aerodynamic resistances
 * DESCRIP-END.
 * FUNCTIONS:    CalcAerodynamic()
 * COMMENTS:     Modified for use with the vicNl model 3-12-98
 *		 by Keith Cherkauer
 */

#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
#include <math.h>

static char vcid[] = "$Id$";

  
/*****************************************************************************
  Function name: CalcAerodynamic()

  Purpose      : Calculate the aerodynamic resistance for each vegetation 
                 layer, and the wind 2m above the layer boundary.  In case of 
                 an overstory, also calculate the wind in the overstory.
                 The values are normalized based on a reference height wind 
                 speed, Uref, of 1 m/s.  To get wind speeds and aerodynamic 
                 resistances for other values of Uref, you need to multiply 
                 the here calculated wind speeds by Uref and divide the 
                 here calculated aerodynamic resistances by Uref
                 
  Required     :
    int NVegLayers - Number of vegetation layers
    char OverStory - flag for presence of overstory.  Only used if NVegLayers 
                     is equal to 1
    double Zref[0]     - Reference height for windspeed
    double n        - Attenuation coefficient for wind in the overstory
    double Height  - Height of the vegetation layers (top layer first)
    double Trunk    - Multiplier for Height[0] that indictaes the top of the 
                     trunk space
    double *U      - Vector of length 3, contains wind for vegetation 
                     conditions listed below:
                     [0] is always wind speed for the snow-free case,
                     [2] is always wind speed for the snow covered case,
                     [1] will contain the wind speed in the canopy if
                     OverStory is TRUE, otherwise it is unused.
    double *Ra     - Vector of length 3, contains aerodynamic resistance 
                     values for the conditions outlined for *U.  
    double *Zref   - Vector of length 3, contains reference height 
                     values for the conditions outlined for *U.  
    double *Z0     - Vector of length 3, contains roughness length  
                     values for the conditions outlined for *U.  
    double *d      - Vector of length 3, contains displacement height 
                     values for the conditions outlined for *U.  

  Returns      : int

  Modifies     :
    double *U
    double *Ra     
    double *Zref     
    double *Z0   
    double *d     
   
  Comments     :
*****************************************************************************/
int  CalcAerodynamic(char    OverStory,     /* overstory flag */
                     double  Height,        /* vegetation height */
                     double  Trunk,         /* trunk ratio parameter */
		     double  Z0_SNOW,       /* snow roughness */
		     double  Z0_SOIL,       /* soil roughness */
                     double  n,             /* wind attenuation parameter */
                     double *Ra,            /* aerodynamic resistances */
                     double *U,             /* adjusted wind speed */
                     double *displacement,  /* vegetation displacement */
                     double *ref_height,    /* vegetation reference height */
                     double *roughness)     /* vegetation roughness */
{
  /******************************************************************
  Modifications:
  2007-Apr-04 Modified to catch and return error flags from surface_fluxes
              subroutine.                                      GCT/KAC
  2009-Jun-09 Modified to use extension of veg_lib structure to contain
	      bare soil information.				TJB
  *******************************************************************/


  double d_Lower;
  double d_Upper;
  double K2;
  double Uh;
  double Ut;
  double Uw;
  double Z0_Lower;
  double Z0_Upper;
  double Zt;
  double Zw;
  double tmp_wind;

  tmp_wind = U[0];

  K2 = von_K * von_K;
  
  /* No OverStory, thus maximum one soil layer */
  
  if ( OverStory == FALSE ) {
    
    /* vegetation cover */
    Z0_Lower = roughness[0];
    d_Lower  = displacement[0];
    
    /* No snow */
    U[0]  = log((2. + Z0_Lower)/Z0_Lower)/log((ref_height[0] - d_Lower)/Z0_Lower);
    /****** DHSVM ******
    Ra[0] = log((2. + Z0_Lower)/Z0_Lower) * log((ref_height[0] - d_Lower)/Z0_Lower)
            /K2;
    ***** Old VIC *****/
    Ra[0] = log((2. + (1.0/0.63 - 1.0) * d_Lower) / Z0_Lower)
          * log((2. + (1.0/0.63 - 1.0) * d_Lower) / (0.1*Z0_Lower)) / K2;
    /******************/

    /* Copy bare parameters into canopy top parameters */
    U[1] = U[0];
    Ra[1] = Ra[0];

    /** Set aerodynamic resistance terms for canopy */
    /* not currently used */
    ref_height[1]   = ref_height[0];
    roughness[1]    = roughness[0];
    displacement[1] = displacement[0];

    /* Snow */
    U[2] = log((2. + Z0_SNOW)/Z0_SNOW)/log(ref_height[0]/Z0_SNOW);
    Ra[2] = log((2. + Z0_SNOW)/Z0_SNOW) * log(ref_height[0]/Z0_SNOW)/K2;
    /** Set aerodynamic resistance terms for snow */
    ref_height[2]   = 2. + Z0_SNOW;
    roughness[2]    = Z0_SNOW;
    displacement[2] = 0.;

  }
  
  /* Overstory present, one or two vegetation layers possible */
  else {
    Z0_Upper = roughness[0];
    d_Upper  = displacement[0];
    
    Z0_Lower = Z0_SOIL;
    d_Lower  = 0; 
    
    Zw = 1.5 * Height - 0.5 * d_Upper;
    Zt = Trunk * Height;

    if (Zt < (Z0_Lower+d_Lower)) {
      fprintf(stderr,"ERROR: CalcAerodynamic - Trunk space height below \"center\" of lower boundary");
      return( ERROR );
    }

    /* Resistance for overstory */
    Ra[1] = log((ref_height[0]-d_Upper)/Z0_Upper)/K2
          * (Height/(n*(Zw-d_Upper)) 
          * (exp(n*(1-(d_Upper+Z0_Upper)/Height))-1)
          + (Zw-Height)/(Zw-d_Upper)
          + log((ref_height[0]-d_Upper)/(Zw-d_Upper)));
    
    /* Wind at different levels in the profile */
    Uw = log((Zw-d_Upper)/Z0_Upper) / log((ref_height[0]-d_Upper)/Z0_Upper);
    Uh = Uw - (1-(Height-d_Upper)/(Zw-d_Upper))
       / log((ref_height[0]-d_Upper)/Z0_Upper);
    U[1] = Uh * exp(n * ((Z0_Upper+d_Upper)/Height - 1.));
    Ut = Uh * exp(n * (Zt/Height - 1.));
    
    /* resistance at the lower boundary */
    

    /***** Old VIC *****/
    U[0]  = log((2. + Z0_Upper)/Z0_Upper)/log((ref_height[0] - d_Upper)/Z0_Upper);
    Ra[0] = log((2. + (1.0/0.63 - 1.0) * d_Upper) / Z0_Upper)
          * log((2. + (1.0/0.63 - 1.0) * d_Upper) / (0.1*Z0_Upper)) / K2;
    /******************/



    /* Snow */
    /* case 1: the wind profile to a height of 2m above the lower boundary is 
       entirely logarithmic */
    if (Zt > (2. + Z0_SNOW)) {
      U[2] = Ut*log((2.+Z0_SNOW)/Z0_SNOW)/log(Zt/Z0_SNOW);
      Ra[2] = log((2.+Z0_SNOW)/Z0_SNOW) * log(Zt/Z0_SNOW)/(K2*Ut);  
    }
    
    /* case 2: the wind profile to a height of 2m above the lower boundary 
       is part logarithmic and part exponential, but the top of the overstory 
       is more than 2 m above the lower boundary */
    else if (Height > (2. + Z0_SNOW)) {
      U[2] = Uh * exp(n * ((2. + Z0_SNOW)/Height - 1.));
      Ra[2] = log(Zt/Z0_SNOW) * log(Zt/Z0_SNOW)/
        (K2*Ut) +
        Height * log((ref_height[0]-d_Upper)/Z0_Upper) / (n*K2*(Zw-d_Upper)) *
        (exp(n*(1-Zt/Height)) - exp(n*(1-(Z0_SNOW+2.)/Height)));
    }
    
    /* case 3: the top of the overstory is less than 2 m above the lower 
       boundary.  The wind profile above the lower boundary is part 
       logarithmic and part exponential, but only extends to the top of the 
       overstory */
    else {
      U[2] = Uh;
      Ra[2] = log(Zt/Z0_SNOW) * log(Zt/Z0_SNOW)/
        (K2*Ut) +
        Height * log((ref_height[0]-d_Upper)/Z0_Upper) / (n*K2*(Zw-d_Upper)) *
        (exp(n*(1-Zt/Height)) - 1);
      fprintf(stderr, "WARNING:  Top of overstory is less than 2 meters above the lower boundary\n");
    }

    /** Set aerodynamic resistance terms for canopy */
    /* not currently used */
    ref_height[1]   = ref_height[0];
    roughness[1]    = roughness[0];
    displacement[1] = displacement[0];
    ref_height[0]   = 2.;
    roughness[0]    = Z0_Lower;
    displacement[0] = d_Lower;

    /** Set aerodynamic resistance terms for snow */
    ref_height[2]   = 2. + Z0_SNOW;
    roughness[2]    = Z0_SNOW;
    displacement[2] = 0.;

  }

  if ( tmp_wind > 0. ) {
    U[0] *= tmp_wind;
    Ra[0] /= tmp_wind;
    if(U[1]!=-999) {
      U[1] *= tmp_wind;
      Ra[1] /= tmp_wind;
    }
    if(U[2]!=-999) {
      U[2] *= tmp_wind;
      Ra[2] /= tmp_wind;
    }
  }
  else {
    U[0] *= tmp_wind;
    Ra[0] = HUGE_RESIST;
    if(U[1]!=-999)
      U[1] *= tmp_wind;
    Ra[1] = HUGE_RESIST;
    if(U[2]!=-999)
      U[2] *= tmp_wind;
    Ra[2] = HUGE_RESIST;
  }
  return (0);

}
