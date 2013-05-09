#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

static char vcid[] = "$Id$";

double func_canopy_energy_bal(double Tfoliage, va_list ap)
/*********************************************************************
  func_canopy_energy_bal    Keith Cherkauer         January 27, 2001

  This routine iterates to determine the temperature of the canopy,
  and solve the resulting fluxes between the canopy and the atmosphere 
  and the canopy and the ground.

  Modifications:
  2004-Sep-28 Added Ra_used to store the aerodynamic resistance used in
	      flux calculations.					TJB
  2009-Jan-16 Modified aero_resist_used and Ra_used to become arrays of
	      two elements (surface and overstory); added
	      options.AERO_RESIST_CANSNOW.				TJB
  2009-May-17 Added AR_406_LS to options.AERO_RESIST_CANSNOW.		TJB
  2009-Sep-14 Replaced 0.622 with EPS in equation for vapor flux.	TJB
  2012-Jan-28 Removed AR_COMBO and GF_FULL.				TJB
 ********************************************************************/
{

  extern option_struct   options;

  /* General Model Parameters */
  int     band;
  int     month;
  int     rec;

  double  delta_t;
  double  elevation;

  double *Wcr;
  double *Wpwp;
  double *depth;
#if SPATIAL_FROST
  double *frost_fract;
#endif

  /* Atmopheric Condition and Forcings */
  double  AirDens;
  double  EactAir;
  double  Press;
  double  Le;
  double  Tcanopy;
  double  Vpd;
  double  mu;

  double *Evap;
  double *Ra;
  double *Ra_used;
  double *Rainfall;
  double *Wind;

  /* Vegetation Terms */
  int     UnderStory;
  int     iveg;
  int     veg_class;

  double *displacement;
  double *ref_height;
  double *roughness;

  float  *root;

  /* Water Flux Terms */
  double  IntRain;
  double  IntSnow;

  double *Wdew;

  layer_data_struct *layer_wet;
  layer_data_struct *layer_dry;
  veg_var_struct    *veg_var_wet;
  veg_var_struct    *veg_var_dry;

  /* Energy Flux Terms */
  double  LongOverIn;
  double  LongUnderOut;
  double  NetShortOver;

  double *AdvectedEnergy;
  double *LatentHeat;
  double *LatentHeatSub;
  double *LongOverOut;
  double *NetLongOver;
  double *NetRadiation;
  double *RefreezeEnergy;
  double *SensibleHeat;
  double *VaporMassFlux;

  /* Internal Variables */
  double  EsSnow;
  double  LongOut;
  double  Ls;
  double  RestTerm;
  double  Tmp;
  double  prec[2];

  /** Read variables from variable length argument list **/

  /* General Model Parameters */
  band    = (int) va_arg(ap, int);
  month   = (int) va_arg(ap, int);
  rec     = (int) va_arg(ap, int);

  delta_t   = (double) va_arg(ap, double);
  elevation = (double) va_arg(ap, double);

  Wcr         = (double *) va_arg(ap, double *);
  Wpwp        = (double *) va_arg(ap, double *);
  depth       = (double *) va_arg(ap, double *);
#if SPATIAL_FROST
  frost_fract = (double *) va_arg(ap, double *);
#endif

  /* Atmopheric Condition and Forcings */
  AirDens  = (double) va_arg(ap, double);
  EactAir  = (double) va_arg(ap, double);
  Press    = (double) va_arg(ap, double);
  Le       = (double) va_arg(ap, double);
  Tcanopy     = (double) va_arg(ap, double);
  Vpd      = (double) va_arg(ap, double);
  mu       = (double) va_arg(ap, double);

  Evap     = (double *) va_arg(ap, double *);
  Ra       = (double *) va_arg(ap, double *);
  Ra_used  = (double *) va_arg(ap, double *);
  Rainfall = (double *) va_arg(ap, double *);
  Wind     = (double *) va_arg(ap, double *);

  /* Vegetation Terms */
  UnderStory = (int) va_arg(ap, int);
  iveg       = (int) va_arg(ap, int);
  veg_class  = (int) va_arg(ap, int);

  displacement = (double *) va_arg(ap, double *);
  ref_height   = (double *) va_arg(ap, double *);
  roughness    = (double *) va_arg(ap, double *);

  root = (float *) va_arg(ap, float *);

  /* Water Flux Terms */
  IntRain = (double) va_arg(ap, double);
  IntSnow = (double) va_arg(ap, double);

  Wdew    = (double *) va_arg(ap, double *);

  layer_wet   = (layer_data_struct *) va_arg(ap, layer_data_struct *);
  layer_dry   = (layer_data_struct *) va_arg(ap, layer_data_struct *);
  veg_var_wet = (veg_var_struct *) va_arg(ap, veg_var_struct *);
  veg_var_dry = (veg_var_struct *) va_arg(ap, veg_var_struct *);

  /* Energy Flux Terms */
  LongOverIn         = (double) va_arg(ap, double);
  LongUnderOut       = (double) va_arg(ap, double);
  NetShortOver       = (double) va_arg(ap, double);

  AdvectedEnergy     = (double *) va_arg(ap, double *);
  LatentHeat         = (double *) va_arg(ap, double *);
  LatentHeatSub      = (double *) va_arg(ap, double *);
  LongOverOut        = (double *) va_arg(ap, double *);
  NetLongOver        = (double *) va_arg(ap, double *);
  NetRadiation       = (double *) va_arg(ap, double *);
  RefreezeEnergy     = (double *) va_arg(ap, double *);
  SensibleHeat       = (double *) va_arg(ap, double *);
  VaporMassFlux      = (double *) va_arg(ap, double *);

  /* Calculate the net radiation at the canopy surface, using the canopy 
     temperature.  The outgoing longwave is subtracted twice, because the 
     canopy radiates in two directions */

  Tmp = Tfoliage + KELVIN;
  *LongOverOut = STEFAN_B * (Tmp * Tmp * Tmp * Tmp);
  *NetRadiation = NetShortOver + LongOverIn + LongUnderOut 
    - 2 * (*LongOverOut);

  *NetLongOver  = LongOverIn - (*LongOverOut);

  if ( IntSnow > 0 ) {

    Ra_used[0] = Ra[0];
    Ra_used[1] = Ra[1];

    /** Added multiplication by 10 to incorporate change in canopy resistance due
	to smoothing by intercepted snow **/
    if (options.AERO_RESIST_CANSNOW == AR_406 || options.AERO_RESIST_CANSNOW == AR_406_LS || options.AERO_RESIST_CANSNOW == AR_406_FULL)
      Ra_used[1] *= 10.;

    /** Calculate the vapor mass flux between intercepted snow in 
	the canopy and the surrounding air mass **/
  
    EsSnow = svp(Tfoliage); 
    
    /* Apply stability correction to aerodynamic resistance */
    if (options.AERO_RESIST_CANSNOW == AR_410) {
      if (Wind[1] > 0.0) {
        Ra_used[1] /= StabilityCorrection(ref_height[1], displacement[1], Tfoliage, 
				          Tcanopy, Wind[1], roughness[1]);
      }
      else
        Ra_used[1] = HUGE_RESIST;
    }

    *VaporMassFlux = AirDens * ( EPS / Press ) * (EactAir - EsSnow) 
      / Ra_used[1] / RHO_W; 

    if (Vpd == 0.0 && *VaporMassFlux < 0.0)
      *VaporMassFlux = 0.0;
  
    /* Calculate the latent heat flux */
    
    Ls = (677. - 0.07 * Tfoliage) * JOULESPCAL * GRAMSPKG;
    *LatentHeatSub = Ls * *VaporMassFlux * RHO_W;
    *LatentHeat = 0;
    *Evap = 0;
    veg_var_wet->throughfall = 0;

    if (options.AERO_RESIST_CANSNOW == AR_406)
      Ra_used[1] /= 10;
  }
  else {

    if (options.AERO_RESIST_CANSNOW == AR_406_FULL || options.AERO_RESIST_CANSNOW == AR_410) {
      Ra_used[0] = Ra[0];
      Ra_used[1] = Ra[1];
    }
    else {
      Ra_used[0] = Ra[0];
      Ra_used[1] = Ra[0];
    }

    Wdew[WET] = IntRain * 1000.;
    prec[WET] = *Rainfall * 1000;
    prec[DRY] = 0;
    *Evap = canopy_evap(layer_wet, layer_dry, veg_var_wet, veg_var_dry, FALSE, 
			veg_class, month, mu, Wdew, delta_t, *NetRadiation, 
			Vpd, NetShortOver, Tcanopy, Ra_used[1], displacement[1], 
			roughness[1], ref_height[1], elevation, prec, 
			depth, Wcr, Wpwp, 
#if SPATIAL_FROST
			frost_fract,
#endif
			root);
    Wdew[WET] /= 1000.;

    *LatentHeat = Le * *Evap * RHO_W;
    *LatentHeatSub = 0;

  }

  /* Calculate the sensible heat flux */

  *SensibleHeat = AirDens * Cp * (Tcanopy - Tfoliage) / Ra_used[1];

  /* Calculate the advected energy */

  *AdvectedEnergy = (4186.8 * Tcanopy * Rainfall[0]) / (delta_t);

  /* Calculate the amount of energy available for refreezing */
  
  RestTerm = *SensibleHeat + *LatentHeat + *LatentHeatSub + *NetRadiation 
    + *AdvectedEnergy;
  
  if ( IntSnow > 0 ) {
    /* Intercepted snow present, check if excess energy can be used to 
       melt or refreeze it */

    *RefreezeEnergy = (IntRain * Lf * RHO_W) / (delta_t);

    if (Tfoliage == 0.0 && RestTerm > -(*RefreezeEnergy)) {
      *RefreezeEnergy = -RestTerm;  /* available energy input over cold content
				       used to melt, i.e. Qrf is negative value
				       (energy out of pack)*/ 
      RestTerm = 0.0;
    }
    else {
      RestTerm += *RefreezeEnergy; /* add this positive value to the pack */
    }

  }
  else *RefreezeEnergy = 0;

  return (RestTerm);

}
