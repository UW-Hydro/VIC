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

double canopy_evap(layer_data_struct *layer_wet,
                   layer_data_struct *layer_dry,
                   veg_var_struct    *veg_var_wet, 
                   veg_var_struct    *veg_var_dry, 
		   char               CALC_EVAP,
                   int                veg_class, 
                   int                month, 
                   double             mu,
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
                   double            *prec,
		   double            *depth,
		   double            *Wcr,
		   double            *Wpwp,
#if SPATIAL_FROST
		   double            *frost_fract,
#endif
		   float             *root)
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
**********************************************************************/ 

{

  /** declare global variables **/
  extern veg_lib_struct *veg_lib; 
  extern option_struct options;

  /** declare local variables **/
  int                Ndist;
  int                dist;
  int                i;
  double             ppt;		/* effective precipitation */
  double             f;		/* fraction of time step used to fill canopy */
  double             throughfall;
  double             Evap;
  double             tmp_Evap;
  double             canopyevap;
  double             tmp_Wdew;
  double             layerevap[MAX_LAYERS];
  double             rc;
  layer_data_struct *tmp_layer;
  veg_var_struct    *tmp_veg_var;

  if(options.DIST_PRCP) Ndist = 2;
  else Ndist = 1;

  Evap = 0;

  for ( dist = 0; dist < Ndist; dist++ ) {

    /* Initialize variables */
    for ( i = 0; i < options.Nlayer; i++ ) layerevap[i] = 0;
    canopyevap = 0;
    throughfall = 0;

    /* Set parameters for distributed precipitation */
    if(dist==0) {
      tmp_layer   = layer_wet;
      tmp_veg_var = veg_var_wet;
      ppt         = prec[WET];
      tmp_Wdew    = Wdew[WET];
    }
    else {
      tmp_layer   = layer_dry;
      tmp_veg_var = veg_var_dry;
      ppt         = prec[DRY];
      mu          = (1. - mu);
      tmp_Wdew    = Wdew[DRY];
    }      

    if(mu > 0) {

      /****************************************************
        Compute Evaporation from Canopy Intercepted Water
      ****************************************************/

      /** Due to month changes ..... Wdmax based on LAI **/
      tmp_veg_var->Wdew = tmp_Wdew;
      if (tmp_Wdew > veg_lib[veg_class].Wdmax[month-1]) {
	throughfall = tmp_Wdew - veg_lib[veg_class].Wdmax[month-1];
	tmp_Wdew    = veg_lib[veg_class].Wdmax[month-1];
      }
      
      rc = calc_rc((double)0.0, net_short, veg_lib[veg_class].RGL,
		   air_temp, vpd, veg_lib[veg_class].LAI[month-1],
		   (double)1.0, FALSE);
      canopyevap = pow((tmp_Wdew / veg_lib[veg_class].Wdmax[month-1]),(2.0/3.0))
		   * penman(air_temp, elevation, rad,
			    vpd, ra, rc, veg_lib[veg_class].rarc)
		   * delta_t / SEC_PER_DAY;

      if (canopyevap > 0.0 && delta_t == SEC_PER_DAY)
	/** If daily time step, evap can include current precipitation **/
	f = min(1.0,((tmp_Wdew + ppt) / canopyevap));
      else if (canopyevap > 0.0)
	/** If sub-daily time step, evap can not exceed current storage **/
	f = min(1.0,((tmp_Wdew) / canopyevap));
      else
	f = 1.0;
      canopyevap *= f;
    
      tmp_Wdew += ppt - canopyevap;
      if (tmp_Wdew < 0.0) 
	tmp_Wdew = 0.0;
      if (tmp_Wdew <= veg_lib[veg_class].Wdmax[month-1]) 
	throughfall += 0.0;
      else {
	throughfall += tmp_Wdew - veg_lib[veg_class].Wdmax[month-1];
	tmp_Wdew = veg_lib[veg_class].Wdmax[month-1];
      }

      /*******************************************
        Compute Evapotranspiration from Vegetation
      *******************************************/
      if(CALC_EVAP)
	transpiration(tmp_layer, veg_class, month, rad, vpd, net_short, 
		      air_temp, ra, ppt, f, delta_t, tmp_veg_var->Wdew, 
		      elevation, depth, Wcr, Wpwp, &tmp_Wdew, &canopyevap, 
		      layerevap, 
#if SPATIAL_FROST
		      frost_fract, 
#endif
		      root);

    }

    tmp_veg_var->canopyevap = canopyevap;
    tmp_veg_var->throughfall = throughfall;
    tmp_veg_var->Wdew = tmp_Wdew;
    tmp_Evap = canopyevap;
    for(i=0;i<options.Nlayer;i++) {
      tmp_layer[i].evap  = layerevap[i];
      tmp_Evap          += layerevap[i];
    }
    
    Evap += tmp_Evap * mu / (1000. * delta_t);

  }

  return (Evap);

}

/**********************************************************************
	EVAPOTRANSPIRATION ROUTINE
**********************************************************************/

void transpiration(layer_data_struct *layer,
		   int veg_class, 
		   int month, 
		   double rad,
		   double vpd,
		   double net_short,
		   double air_temp,
		   double ra,
		   double ppt,
		   double f,
		   double delta_t,
		   double Wdew,
		   double elevation,
		   double *depth,
		   double *Wcr,
		   double *Wpwp,
		   double *new_Wdew,
		   double *canopyevap,
		   double *layerevap,
#if SPATIAL_FROST
		   double *frost_fract,
#endif
		   float  *root)
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

**********************************************************************/
{
  extern veg_lib_struct *veg_lib;
  extern option_struct options;

  int    i;
#if SPATIAL_FROST
  int    frost_area;
#endif
  double gsm_inv;               	/* soil moisture stress factor */
  double moist1, moist2;                /* tmp holding of moisture */
  double evap;                          /* tmp holding for evap total */
  double Wcr1;                          /* tmp holding of critical water for upper layers */
  double root_sum;                      /* proportion of roots in moist>Wcr zones */
  double spare_evap;                    /* evap for 2nd distribution */
  double avail_moist[MAX_LAYERS];         /* moisture available for trans */
  double ice[MAX_LAYERS];
  double rc;

  /********************************************************************** 
     EVAPOTRANSPIRATION

     Calculation of the evapotranspirations
     2.18

     First part: Soil moistures and root fractions of both layers
     influence each other

     Re-written to allow for multi-layers.
  **********************************************************************/
 
  /**************************************************
    Set ice content in all individual layers
    **************************************************/
  for(i=0;i<options.Nlayer;i++){
#if SPATIAL_FROST
    ice[i] = 0;
    for ( frost_area = 0; frost_area < FROST_SUBAREAS; frost_area++ ) {
      ice[i] += layer[i].ice[frost_area] * frost_fract[frost_area];
    }
#else
    ice[i]         = layer[i].ice;
#endif
  }

  /**************************************************
    Compute moisture content in combined upper layers
    **************************************************/
  moist1 = 0.0;
  Wcr1 = 0.0;  
  for(i=0;i<options.Nlayer-1;i++){
    if(root[i] > 0.) {
#if SPATIAL_FROST
      avail_moist[i] = 0;
      for ( frost_area = 0; frost_area < FROST_SUBAREAS; frost_area++ ) {
	avail_moist[i] += ((layer[i].moist - layer[i].ice[frost_area]) 
			   * frost_fract[frost_area]);
      }
#else
      avail_moist[i] = layer[i].moist - layer[i].ice;
#endif
      moist1+=avail_moist[i];
      Wcr1 += Wcr[i];
    }
    else avail_moist[i]=0.;
  }

  /*****************************************
    Compute moisture content in lowest layer
    *****************************************/
  i = options.Nlayer - 1;
#if SPATIAL_FROST
  moist2 = 0;
  for ( frost_area = 0; frost_area < FROST_SUBAREAS; frost_area++ )
    moist2 += ((layer[i].moist - layer[i].ice[frost_area]) 
	       * frost_fract[frost_area]);
#else
  moist2 = layer[i].moist - layer[i].ice;
#endif

  avail_moist[i]=moist2;

  /******************************************************************
    CASE 1: Moisture in both layers exceeds Wcr, or Moisture in
    layer with more than half of the roots exceeds Wcr.

    Potential evapotranspiration not hindered by soil dryness.  If
    layer with less than half the roots is dryer than Wcr, extra
    evaporation is taken from the wetter layer.  Otherwise layers
    contribute to evapotransipration based on root fraction.
  ******************************************************************/

  if( (moist1>=Wcr1 && moist2>=Wcr[options.Nlayer-1] && Wcr1>0.) ||
      (moist1>=Wcr1 && (1-root[options.Nlayer-1])>= 0.5) ||
      (moist2>=Wcr[options.Nlayer-1] &&
      root[options.Nlayer-1]>=0.5) ){
    gsm_inv=1.0;
    rc = calc_rc(veg_lib[veg_class].rmin, net_short, veg_lib[veg_class].RGL,
		 air_temp, vpd, veg_lib[veg_class].LAI[month-1], gsm_inv, FALSE);
    evap = penman(air_temp, elevation, rad, vpd, ra, rc,
		  veg_lib[veg_class].rarc)
	   * delta_t / SEC_PER_DAY
	   * (1.0-f*pow((Wdew/veg_lib[veg_class].Wdmax[month-1]), (2.0/3.0)));

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

    for(i=0;i<options.Nlayer;i++){
      /** Set evaporation restriction factor **/
      if(avail_moist[i] >= Wcr[i])
	gsm_inv=1.0;
      else if(avail_moist[i] >= Wpwp[i])
	gsm_inv=(avail_moist[i] - Wpwp[i]) /
	  (Wcr[i] - Wpwp[i]);
      else 
	gsm_inv=0.0;

      if(gsm_inv > 0.0){
	/** Compute potential evapotranspiration **/
        rc = calc_rc(veg_lib[veg_class].rmin, net_short, veg_lib[veg_class].RGL,
		     air_temp, vpd, veg_lib[veg_class].LAI[month-1], gsm_inv, FALSE);
        layerevap[i] = penman(air_temp, elevation, rad, vpd, ra, rc,
			      veg_lib[veg_class].rarc)
		       * delta_t / SEC_PER_DAY * (double)root[i]
		       * (1.0-f*pow((Wdew/veg_lib[veg_class].Wdmax[month-1]),
					 (2.0/3.0)));
      }
      else layerevap[i] = 0.0;

    }
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
