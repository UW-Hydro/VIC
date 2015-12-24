#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

static char vcid[] = "$Id$";

/**********************************************************************
	canopy_evap.c	Dag Lohmann		September 1995

  This routine computes the evaporation, traspiration and throughfall
  of the vegetation types for multi-layered model.

  The value of x, the fraction of precipitation that exceeds the 
  canopy storage capacity, is returned by the subroutine.

  UNITS:	moist (mm)
		evap (mm)
		prec (mm)
		melt (mm)

  VARIABLE TYPE        NAME          UNITS DESCRIPTION
  atmos_data_struct    atmos         N/A   atmospheric forcing data structure
  layer_data_struct   *layer         N/A   soil layer variable structure
  veg_var_struct      *veg_var       N/A   vegetation variable structure
  char                 CALC_EVAP     N/A   TRUE = calculate evapotranspiration
  int                  veg_class     N/A   vegetation class index number
  int                  month         N/A   current month
  global_param_struct  global        N/A   global parameter structure
  double               mu            fract wet (or dry) fraction of grid cell
  double               ra            s/m   aerodynamic resistance
  double               prec          mm    precipitation
  double               displacement  m     displacement height of surface cover
  double               roughness     m     roughness height of surface cover
  double               ref_height    m     measurement reference height

  Modifications:
  9/1/97	Greg O'Donnell
  4-12-98  Code cleaned and final version prepared, KAC
  06-25-98 modified for new distributed precipitation data structure KAC
  01-19-00 modified to function with new simplified soil moisture 
           scheme                                                  KAC
  5-8-2001 Modified to close the canopy energy balance.       KAC

**********************************************************************/

double canopy_evap(layer_data_struct *layer,
                   veg_var_struct    *veg_var, 
		   char               CALC_EVAP,
                   int                veg_class, 
                   int                month, 
                   double            *Wdew,
                   double             delta_t,
                   double             rad,
		   double             vpd,
		   double             net_short,
		   double             air_temp,
                   double             ra,
                   double             displacement,
                   double             roughness,
                   double             ref_height,
		   double             elevation,
                   double             ppt,
		   double            *depth,
		   double            *Wmax,
		   double            *Wcr,
		   double            *Wpwp,
		   double            *frost_fract,
		   float             *root,
		   double            *dryFrac,
		   double             shortwave,
		   double             Catm,
		   double            *CanopLayerBnd)
/********************************************************************** 
  CANOPY EVAPORATION

  Calculation of evaporation from the canopy, including the
  possibility of potential evaporation exhausting ppt+canopy storage
  2.16 + 2.17
  Index [0] refers to current time step, index [1] to next one
  If f < 1.0 then veg_var->canopyevap = veg_var->Wdew + ppt
  and Wdew = 0.0

  DEFINITIONS:
  Wdmax - max monthly dew holding capacity
  Wdew - dew trapped on vegetation

  Modified 
     04-14-98 to work within calc_surf_energy_balance.c  KAC
     07-24-98 fixed problem that caused hourly precipitation
              to evaporate from the canopy during the same
	      time step that it falls (OK for daily time step, 
	      but causes oscilations in surface temperature
	      for hourly time step)                      KAC, Dag
	        modifications:
     6-8-2000 Modified to use spatially distributed soil frost if 
              present.                                           KAC
     5-8-2001 Modified to close the canopy energy balance.       KAC
  2009-Jun-09 Moved computation of canopy resistance rc out of penman()
	      and into separate function calc_rc().			TJB
  2012-Jan-16 Removed LINK_DEBUG code					BN
  2013-Jul-25 Save dryFrac for use elsewhere.				TJB
  2014-Mar-28 Removed DIST_PRCP option.					TJB
  2014-Apr-25 Switched LAI from veg_lib to veg_var.			TJB
  2014-May-05 Added logic to handle LAI = 0.				TJB
  2014-May-05 Replaced HUGE_RESIST with RSMAX.				TJB
**********************************************************************/ 

{

  /** declare global variables **/
  extern veg_lib_struct *veg_lib; 
  extern option_struct options;

  /** declare local variables **/
  int                i;
  double             f;		/* fraction of time step used to fill canopy */
  double             throughfall;
  double             Evap;
  double             tmp_Evap;
  double             canopyevap;
  double             tmp_Wdew;
  double             layerevap[MAX_LAYERS];
  double             rc;
  int flag_irr;

  Evap = 0;

  /* Initialize variables */
  for ( i = 0; i < options.Nlayer; i++ ) layerevap[i] = 0;
  canopyevap = 0;
  throughfall = 0;
  tmp_Wdew = *Wdew;

  /****************************************************
    Compute Evaporation from Canopy Intercepted Water
  ****************************************************/

  veg_var->Wdew = tmp_Wdew;
  if (tmp_Wdew > veg_var->Wdmax) {
    throughfall = tmp_Wdew - veg_var->Wdmax;
    tmp_Wdew    = veg_var->Wdmax;
  }

  flag_irr = veg_lib[veg_class].irr_active[month-1];

  rc = calc_rc((double)0.0, net_short, veg_lib[veg_class].RGL,
               air_temp, vpd, veg_var->LAI, (double)1.0, FALSE,flag_irr);
  if (veg_var->LAI > 0)
    canopyevap = pow((tmp_Wdew/veg_var->Wdmax),(2.0/3.0))
                 * penman(air_temp, elevation, rad, vpd, ra, rc, veg_lib[veg_class].rarc)
                 * delta_t / SEC_PER_DAY;
  else
    canopyevap = 0;

  if (canopyevap > 0.0 && delta_t == SEC_PER_DAY)
    /** If daily time step, evap can include current precipitation **/
    f = min(1.0,((tmp_Wdew + ppt) / canopyevap));
  else if (canopyevap > 0.0)
    /** If sub-daily time step, evap can not exceed current storage **/
    f = min(1.0,((tmp_Wdew) / canopyevap));
  else
    f = 1.0;
  canopyevap *= f;

  /* compute fraction of canopy that is dry */
  if (veg_var->Wdmax > 0)
    *dryFrac = 1.0-f*pow((tmp_Wdew/veg_var->Wdmax),(2.0/3.0));
  else
    *dryFrac = 0;

  tmp_Wdew += ppt - canopyevap;
  if (tmp_Wdew < 0.0) 
    tmp_Wdew = 0.0;
  if (tmp_Wdew <= veg_var->Wdmax) 
    throughfall += 0.0;
  else {
    throughfall += tmp_Wdew - veg_var->Wdmax;
    tmp_Wdew = veg_var->Wdmax;
  }

  /*******************************************
    Compute Evapotranspiration from Vegetation
  *******************************************/
  if(CALC_EVAP)
    transpiration(layer, veg_var, veg_class, month, rad, vpd, net_short, 
                  air_temp, ra, ppt, *dryFrac, delta_t, 
                  elevation, depth, Wmax, Wcr, Wpwp, 
                  layerevap, 
                  frost_fract, 
                  root,
                  shortwave,
                  Catm,
                  CanopLayerBnd);

  veg_var->canopyevap = canopyevap;
  veg_var->throughfall = throughfall;
  veg_var->Wdew = tmp_Wdew;
  tmp_Evap = canopyevap;
  for(i=0;i<options.Nlayer;i++) {
    layer[i].evap  = layerevap[i];
    tmp_Evap          += layerevap[i];
  }
    
  Evap += tmp_Evap / (1000. * delta_t);

  return (Evap);

}

/**********************************************************************
	EVAPOTRANSPIRATION ROUTINE
**********************************************************************/

void transpiration(layer_data_struct *layer,
		   veg_var_struct *veg_var,
		   int veg_class, 
		   int month, 
		   double rad,
		   double vpd,
		   double net_short,
		   double air_temp,
		   double ra,
		   double ppt,
		   double dryFrac,
		   double delta_t,
		   double elevation,
		   double *depth,
		   double *Wmax,
		   double *Wcr,
		   double *Wpwp,
		   double *layerevap,
		   double *frost_fract,
		   float  *root,
                   double  shortwave,
                   double  Catm,
                   double *CanopLayerBnd)
/**********************************************************************
  Computes evapotranspiration for unfrozen soils
  Allows for multiple layers.

  modifications:
  6-8-2000 Modified to use spatially distributed soil frost if 
           present.                                               KAC
  2006-Oct-16 Modified to initialize ice[] for all soil layers
	      before computing available moisture (to avoid using
	      uninitialized values later on).				TJB
  2009-Jun-09 Moved computation of canopy resistance rc out of penman()
	      and into separate function calc_rc().			TJB
  2013-Jul-25 Save dryFrac for use elsewhere.				TJB
  2013-Jul-25 Added photosynthesis terms.				TJB
  2014-Apr-25 Switched LAI from veg_lib to veg_var.			TJB
**********************************************************************/
{
  extern veg_lib_struct *veg_lib;
  extern option_struct options;

  int    i;
  int    frost_area;
  double gsm_inv;               	/* soil moisture stress factor */
  double moist1, moist2;                /* tmp holding of moisture */
  double evap;                          /* tmp holding for evap total */
  double Wcr1;                          /* tmp holding of critical water for upper layers */
  double root_sum;                      /* proportion of roots in moist>Wcr zones */
  double spare_evap;                    /* evap for 2nd distribution */
  double avail_moist[MAX_LAYERS];         /* moisture available for trans */
  double ice[MAX_LAYERS];
  double gc;
  double *gsLayer;
  int    cidx;
  int flag_irr;
  /********************************************************************** 
     EVAPOTRANSPIRATION

     Calculation of the evapotranspirations
     2.18

     First part: Soil moistures and root fractions of both layers
     influence each other

     Re-written to allow for multi-layers.
  **********************************************************************/


  flag_irr = veg_lib[veg_class].irr_active[month-1];

  /**************************************************
    Set ice content in all individual layers
    **************************************************/
  for(i=0;i<options.Nlayer;i++){
    ice[i] = 0;
    for ( frost_area = 0; frost_area < options.Nfrost; frost_area++ ) {
      ice[i] += layer[i].ice[frost_area] * frost_fract[frost_area];
    }
  }

  /**************************************************
    Compute moisture content in combined upper layers
    **************************************************/
  moist1 = 0.0;
  Wcr1 = 0.0;  
  for(i=0;i<options.Nlayer-1;i++){
    if(root[i] > 0.) {
      avail_moist[i] = 0;
      for ( frost_area = 0; frost_area < options.Nfrost; frost_area++ ) {
	avail_moist[i] += ((layer[i].moist - layer[i].ice[frost_area]) * frost_fract[frost_area]);
      }
      moist1+=avail_moist[i];
      Wcr1 += Wcr[i];
    }
    else avail_moist[i]=0.;
  }

  /*****************************************
    Compute moisture content in lowest layer
    *****************************************/
  i = options.Nlayer - 1;
  moist2 = 0;
  for ( frost_area = 0; frost_area < options.Nfrost; frost_area++ )
    moist2 += ((layer[i].moist - layer[i].ice[frost_area]) * frost_fract[frost_area]);
  avail_moist[i]=moist2;

  /** Set photosynthesis inhibition factor **/
  if (layer[0].moist > veg_lib[veg_class].Wnpp_inhib*Wmax[0])
    veg_var->NPPfactor = veg_lib[veg_class].NPPfactor_sat + (1 - veg_lib[veg_class].NPPfactor_sat) * (Wmax[0] - layer[0].moist) / (Wmax[0] - veg_lib[veg_class].Wnpp_inhib*Wmax[0]);
  else
    veg_var->NPPfactor = 1.0;

  /******************************************************************
    CASE 1: Moisture in both layers exceeds Wcr, or Moisture in
    layer with more than half of the roots exceeds Wcr.

    Potential evapotranspiration not hindered by soil dryness.  If
    layer with less than half the roots is dryer than Wcr, extra
    evaporation is taken from the wetter layer.  Otherwise layers
    contribute to evapotransipration based on root fraction.
  ******************************************************************/

  if( options.SHARE_LAYER_MOIST &&
      ( (moist1>=Wcr1 && moist2>=Wcr[options.Nlayer-1] && Wcr1>0.) ||
        (moist1>=Wcr1 && (1-root[options.Nlayer-1])>= 0.5) ||
        (moist2>=Wcr[options.Nlayer-1] && root[options.Nlayer-1]>=0.5) ) ) {

    gsm_inv=1.0;

    /* compute whole-canopy stomatal resistance */
    if (!options.CARBON || options.RC_MODE == RC_JARVIS) {
      /* Jarvis scheme, using resistance factors from Wigmosta et al., 1994 */
      veg_var->rc = calc_rc(veg_lib[veg_class].rmin, net_short,
		   veg_lib[veg_class].RGL, air_temp, vpd,
			    veg_var->LAI, gsm_inv, FALSE,flag_irr);
      if (options.CARBON) {
        for (cidx=0; cidx<options.Ncanopy; cidx++) {
          if (veg_var->LAI > 0)
            veg_var->rsLayer[cidx] = veg_var->rc / veg_var->LAI;
          else
            veg_var->rsLayer[cidx] = HUGE_RESIST;
          if (veg_var->rsLayer[cidx] > RSMAX) veg_var->rsLayer[cidx] = RSMAX;
        }
      }
    }
    else {
      /* Compute rc based on photosynthetic demand from Knorr 1997 */
      calc_rc_ps(veg_lib[veg_class].Ctype, veg_lib[veg_class].MaxCarboxRate,
		 veg_lib[veg_class].MaxETransport, veg_lib[veg_class].CO2Specificity,
		 veg_var->NscaleFactor, air_temp, shortwave, veg_var->aPARLayer, elevation, Catm,
		 CanopLayerBnd, veg_var->LAI, gsm_inv, vpd, veg_var->rsLayer, &(veg_var->rc));
    }

    /* compute transpiration */
    evap = penman(air_temp, elevation, rad, vpd, ra, veg_var->rc,
		  veg_lib[veg_class].rarc) * delta_t / SEC_PER_DAY * dryFrac;

    /** divide up evap based on root distribution **/
    /** Note the indexing of the roots **/
    root_sum=1.0;
    spare_evap=0.0;
    for(i=0;i<options.Nlayer;i++){
      if(avail_moist[i]>=Wcr[i]){
        layerevap[i]=evap*(double)root[i];
      }
      else {
          
        if (avail_moist[i] >= Wpwp[i]) 
          gsm_inv = (avail_moist[i] - Wpwp[i]) /
                    (Wcr[i] - Wpwp[i]);
        else 
          gsm_inv=0.0;
	    
        layerevap[i]  = evap*gsm_inv*(double)root[i];
        root_sum     -= root[i];
        spare_evap    = evap*(double)root[i]*(1.0-gsm_inv);
      }
    }

    /** Assign excess evaporation to wetter layer **/
    if(spare_evap>0.0){
      for(i=0;i<options.Nlayer;i++){
        if(avail_moist[i] >= Wcr[i]){
          layerevap[i] += (double)root[i]*spare_evap/root_sum;
        }
      }
    }
  }

  /*********************************************************************
    CASE 2: Independent evapotranspirations

    Evapotranspiration is restricted by low soil moisture. Evaporation
    is computed independantly from each soil layer.
  *********************************************************************/

  else {

    /* Initialize conductances for aggregation over soil layers */
    gc =  0;
    if (options.CARBON) {
      gsLayer = (double *)calloc(options.Ncanopy, sizeof(double));
      for (cidx=0; cidx<options.Ncanopy; cidx++) {
        gsLayer[cidx] = 0;
      }
    }

    for (i=0;i<options.Nlayer;i++) {

      /** Set evaporation restriction factor **/
      if(avail_moist[i] >= Wcr[i])
	gsm_inv = 1.0;
      else if(avail_moist[i] >= Wpwp[i] && avail_moist[i] < Wcr[i])
	gsm_inv = (avail_moist[i] - Wpwp[i]) / (Wcr[i] - Wpwp[i]);
      else 
	gsm_inv = 0.0;

      if(gsm_inv > 0.0){

        /* compute whole-canopy stomatal resistance */
        if (!options.CARBON || options.RC_MODE == RC_JARVIS) {
          /* Jarvis scheme, using resistance factors from Wigmosta et al., 1994 */
          veg_var->rc = calc_rc(veg_lib[veg_class].rmin, net_short,
		       veg_lib[veg_class].RGL, air_temp, vpd,
				veg_var->LAI, gsm_inv, FALSE,flag_irr);
          if (options.CARBON) {
            for (cidx=0; cidx<options.Ncanopy; cidx++) {
              if (veg_var->LAI > 0)
                veg_var->rsLayer[cidx] = veg_var->rc / veg_var->LAI;
              else
                veg_var->rsLayer[cidx] = HUGE_RESIST;
              if (veg_var->rsLayer[cidx] > RSMAX) veg_var->rsLayer[cidx] = RSMAX;
            }
          }
        }
        else {
          /* Compute rc based on photosynthetic demand from Knorr 1997 */
          calc_rc_ps(veg_lib[veg_class].Ctype, veg_lib[veg_class].MaxCarboxRate,
		     veg_lib[veg_class].MaxETransport, veg_lib[veg_class].CO2Specificity,
		     veg_var->NscaleFactor, air_temp, shortwave, veg_var->aPARLayer, elevation, Catm,
		     CanopLayerBnd, veg_var->LAI, gsm_inv, vpd, veg_var->rsLayer, &(veg_var->rc));
        }

        /* compute transpiration */
        layerevap[i] = penman(air_temp, elevation, rad, vpd, ra, veg_var->rc,
			      veg_lib[veg_class].rarc) * delta_t / SEC_PER_DAY * dryFrac * (double)root[i];

        if (veg_var->rc > 0)
          gc += 1/(veg_var->rc);
        else
          gc += HUGE_RESIST;

        if (options.CARBON) {
          for (cidx=0; cidx<options.Ncanopy; cidx++) {
            if (veg_var->rsLayer[cidx] > 0)
              gsLayer[cidx] += 1/(veg_var->rsLayer[cidx]);
            else
              gsLayer[cidx] += HUGE_RESIST;
          }
        }

      }
      else {
	layerevap[i] = 0.0;
        gc += 0;
        if (options.CARBON) {
          for (cidx=0; cidx<options.Ncanopy; cidx++) {
            gsLayer[cidx] += 0;
          }
	}
      }

    } // end loop over layers

    /* Now, take the inverse of the conductance */
    if (gc > 0)
      veg_var->rc = 1/gc;
    else
      veg_var->rc = HUGE_RESIST;
    if (veg_var->rc > RSMAX) veg_var->rc = RSMAX;

    if (options.CARBON) {
      for (cidx=0; cidx<options.Ncanopy; cidx++) {
        if (gsLayer[cidx] > 0)
          veg_var->rsLayer[cidx] = 1/gsLayer[cidx];
        else
          veg_var->rsLayer[cidx] = HUGE_RESIST;
        if (veg_var->rsLayer[cidx] > RSMAX) veg_var->rsLayer[cidx] = RSMAX;
      }
    }

    if (options.CARBON) free((char *)gsLayer);

  }

  /****************************************************************
    Check that evapotransipration does not cause soil moisture to 
    fall below wilting point.
  ****************************************************************/
  for ( i = 0; i < options.Nlayer; i++ ) {
    if ( ice[i] > 0 ) {
      if ( ice[i] >= Wpwp[i] ) {
	// ice content greater than wilting point can use all unfrozen moist
	if ( layerevap[i] > avail_moist[i] ) layerevap[i] = avail_moist[i];
      }
      else {
	// ice content less than wilting point restrict loss of unfrozen moist
	if ( layerevap[i] > layer[i].moist - Wpwp[i] ) 
	  layerevap[i] = layer[i].moist - Wpwp[i];
      }
    }
    else {
      // No ice restrict loss of unfrozen moist
      if ( layerevap[i] > layer[i].moist - Wpwp[i] ) 
	layerevap[i] = layer[i].moist - Wpwp[i];
    }
    if ( layerevap[i] < 0.0 ) layerevap[i] = 0.0;
  }

}
