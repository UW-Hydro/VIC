#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vicNl.h>

static char vcid[] = "$Id: $";

double func_lake_energy_bal(double Tlakebot, va_list ap)
/**********************************************************************
 * Solves both the T profile of the lake water column and of the soil
 * under the lake, for a given lake bottom temperature Tlb.  Outputs
 * the energy balance error at the lake bottom water-soil interface.
 **********************************************************************/
{

  extern option_struct options;

  int i,j,k;
  int ErrorFlag;
  double Le;
  double windi;
  double windw;
  int mixdepth;
  int freezeflag;
  double T[MAX_LAKE_NODES];
  double Ti[MAX_LAKE_NODES];
  double LWnetw,LWneti;
  double Qhw,Qhi;
  double Qew,Qei;
  double Qnet_ice;
  double qw, qf;
  float              tmp_float;
  double             tmp_double;
  layer_data_struct  tmp_layer;
  veg_var_struct     tmp_veg_var;
  double energy_out_bottom;
  double energy_out_bottom_ice;
  double heat_cond_to_soil;
  double heat_cond_to_soil_ice;
  double energy_ice_formation;
  double energy_ice_melt_bot;
  double temp_refreeze_energy;
  double tmp_deltaH;
  double tmp_fusion;
  double temphi;
  double lake_bottom_tempavg;
  double T2;
  double Ts_old;
  double T1_old;
  double delta_t;
  double D1;
  double D2;
  double max_moist;
  double error;

  /* Variables from variable list */

  /* general model terms */
  int                rec;
  int                nrecs;
  int                dt;

  /* major data structures */
  lake_var_struct   *lake;
  snow_data_struct  *lake_snow;
  soil_con_struct   *soil_con;
  dmy_struct         dmy;

  /* atm forcings */
  double             snowfall;
  double             rainfall;
  double             tair;
  double             wind;
  double             sw_water;
  double             sw_ice;
  double             longin;
  double             vp;
  double             vpd;
  double             pressure;
  double             air_density;

  double             wind_h;

  /* lake water/ice state */
  double             fracprv;
  double            *water_density;
  double            *water_cp;
  double             Tcutoff;
  double             sumjoulb;
  double            *Tavg;
  double            *ice_water_eq;
  double            *new_ice_area;
  double            *new_ice_height;
  double            *new_ice_water_eq;
  double            *lake_tempavg;
  double            *tempi;
  double            *areai;
  double            *hice;

  /* lake soil state */
  int                Nnodes;
  int                NOFLUX;
  int               *FIRST_SOLN;
  double             ice0;
  double            *moist_node;
  double            *ice_node;
  double            *T_node;
  double            *Tnew_node;
  char              *Tnew_fbflag;
  int               *Tnew_fbcount;
  double            *T1;

  /************************************
    Read variables from variable list 
  ************************************/

  /* general model terms */
  rec                  = (int) va_arg(ap, int);
  nrecs                = (int) va_arg(ap, int);
  dt                   = (int) va_arg(ap, int);

  /* major data structures */
  lake                 = (lake_var_struct *) va_arg(ap, lake_var_struct *);
  lake_snow            = (snow_data_struct *) va_arg(ap, snow_data_struct *);
  soil_con             = (soil_con_struct *) va_arg(ap, soil_con_struct *);
  dmy                  = (dmy_struct) va_arg(ap, dmy_struct);

  /* atm forcings */
  snowfall             = (double) va_arg(ap, double);
  rainfall             = (double) va_arg(ap, double);
  tair                 = (double) va_arg(ap, double);
  wind                 = (double) va_arg(ap, double);
  sw_water             = (double) va_arg(ap, double);
  sw_ice               = (double) va_arg(ap, double);
  longin               = (double) va_arg(ap, double);
  vp                   = (double) va_arg(ap, double);
  vpd                  = (double) va_arg(ap, double);
  pressure             = (double) va_arg(ap, double);
  air_density          = (double) va_arg(ap, double);

  wind_h               = (double) va_arg(ap, double);

  /* lake water/ice state */
  fracprv              = (double) va_arg(ap, double);
  water_density        = (double *) va_arg(ap, double *);
  water_cp             = (double *) va_arg(ap, double *);
  Tcutoff              = (double) va_arg(ap, double);
  sumjoulb             = (double) va_arg(ap, double);
  Tavg                 = (double *) va_arg(ap, double *);
  ice_water_eq         = (double *) va_arg(ap, double *);
  new_ice_area         = (double *) va_arg(ap, double *);
  new_ice_height       = (double *) va_arg(ap, double *);
  new_ice_water_eq     = (double *) va_arg(ap, double *);
  lake_tempavg         = (double *) va_arg(ap, double *);
  tempi                = (double *) va_arg(ap, double *);
  areai                = (double *) va_arg(ap, double *);
  hice                 = (double *) va_arg(ap, double *);

  /* lake soil state */
  Nnodes               = (int) va_arg(ap, int);
  NOFLUX               = (int) va_arg(ap, int);
  FIRST_SOLN           = (int *) va_arg(ap, int *);
  ice0                 = (double) va_arg(ap, double);
  moist_node           = (double *) va_arg(ap, double *);
  ice_node             = (double *) va_arg(ap, double *);
  T_node               = (double *) va_arg(ap, double *);
  Tnew_node            = (double *) va_arg(ap, double *);
  Tnew_fbflag          = (char *) va_arg(ap, char *);
  Tnew_fbcount         = (int *) va_arg(ap, int *);
  T1                   = (double *) va_arg(ap, double *);

  /* Initialize iteration variables */
  for ( k = 0; k < lake->activenod; k++ ) {
    T[k] = lake->temp[k];
    Ti[k] = lake->temp[k];
  }
  lake_snow->albedo = lake->snow.albedo;
  lake_snow->coldcontent = lake->snow.coldcontent;
  lake_snow->coverage = lake->snow.coverage;
  lake_snow->density = lake->snow.density;
  lake_snow->depth = lake->snow.depth;
  lake_snow->last_snow = lake->snow.last_snow;
  lake_snow->max_snow_depth = lake->snow.max_snow_depth;
  lake_snow->MELTING = lake->snow.MELTING;
  lake_snow->pack_temp = lake->snow.pack_temp;
  lake_snow->pack_water = lake->snow.pack_water;
  lake_snow->snow = lake->snow.snow;
  lake_snow->store_coverage = lake->snow.store_coverage;
  lake_snow->store_snow = lake->snow.store_snow;
  lake_snow->store_swq = lake->snow.store_swq;
  lake_snow->surf_temp = lake->snow.surf_temp;
  lake_snow->surf_temp_fbflag = lake->snow.surf_temp_fbflag;
  lake_snow->surf_temp_fbcount = lake->snow.surf_temp_fbcount;
  lake_snow->surf_water = lake->snow.surf_water;
  lake_snow->swq = lake->snow.swq;
  lake_snow->snow_distrib_slope = lake->snow.snow_distrib_slope;
  lake->energy.advection=0.0;
  lake->energy.deltaCC = 0.0;
  lake->energy.deltaH = 0.0;
  lake->energy.fusion = 0.0;
  lake->energy.grnd_flux = 0.0;
  lake->energy.snow_flux = 0.0;
  lake->energy.lake_soil_heat_flux = 0.0;
  lake->energy.lake_soil_net_short = 0.0;
  lake->energy.Tsurf = lake->temp[0];
  for (k=0; k<Nnodes; k++) {
    T_node[k] = lake->energy.T[k];
    moist_node[k] = soil_con->max_moist_node[k];
    ice_node[k] = lake->energy.ice[k];
  }
  tmp_deltaH = 0.0;
  tmp_fusion = 0.0;
  lake->snowmlt = 0.0;
  lake->evapw=0.0;
  *tempi = lake->tempi;
  *areai = lake->areai;
  *hice = lake->hice;
  *ice_water_eq = lake->ice_water_eq;
  qw=0.0;
  *new_ice_height = *new_ice_area = *new_ice_water_eq = 0.0;
  energy_ice_formation = 0.0;
  energy_out_bottom = energy_out_bottom_ice = 0.0;
  heat_cond_to_soil = heat_cond_to_soil_ice = 0.0;
  temp_refreeze_energy = 0.0;
  freezeflag = 0;


  /**********************************************************************
   *  Calculate initial energy balance over ice-free water.
   **********************************************************************/
   
  if( (1.-fracprv) > SMALL && lake->activenod > 0) {
    freezeflag=1;         /* Calculation for water, not ice. */
    windw = wind*log((2. + ZWATER)/ZWATER)/log(wind_h/ZWATER);

    ErrorFlag = water_energy_balance(lake->activenod, lake->surface, &lake->evapw, 
				     dt, freezeflag, lake->dz, lake->surfdz,
				     (double)soil_con->lat, Tcutoff, tair, windw, 
				     pressure, vp, air_density, longin, sw_water, 
				     sumjoulb, wind_h, &Qhw, &Qew, &LWnetw, T, 
				     water_density, &lake->energy.deltaH, 
				     &energy_ice_formation, fracprv, 
				     new_ice_area, water_cp, new_ice_height,
				     &energy_out_bottom, new_ice_water_eq,
				     lake->volume-lake->ice_water_eq, soil_con->Zsum_node[1],
				     0.5*(lake->energy.kappa_node[0]+lake->energy.kappa_node[1]),
				     Tlakebot, T_node[1], &heat_cond_to_soil);
    if ( ErrorFlag == ERROR ) return (ERROR);

    /* --------------------------------------------------------------------
     * Do the convective mixing of the lake water.
     * -------------------------------------------------------------------- */

    mixdepth = 0;        /* Set to zero for this time step. */
    tracer_mixer( T, &mixdepth, freezeflag, lake->surface, 
		  lake->activenod, lake->dz, lake->surfdz, water_cp );

    lake->energy.AtmosLatent      = ( 1. - fracprv ) * Qew;
    lake->energy.AtmosSensible    = ( 1. - fracprv ) * Qhw;
    lake->energy.NetLongAtmos     = ( 1. - fracprv ) * LWnetw;
    lake->energy.NetShortAtmos    = ( 1. - fracprv ) * sw_water;
    lake->energy.refreeze_energy  = energy_ice_formation*(1.-fracprv);
    lake->energy.deltaH          *= ( 1. - fracprv );
    lake->energy.Tsurf            = (1. - fracprv)*T[0];
    lake->energy.lake_soil_heat_flux = (heat_cond_to_soil)*(1.-fracprv);
    lake->energy.lake_soil_net_short = (energy_out_bottom)*(1.-fracprv);

  }          /* End of water fraction calculations. */
  else {
    // ice covers 100% of lake, reset open water fluxes
    mixdepth = 0;
    LWnetw = 0;
    Qew    = 0;
    Qhw    = 0;

    lake->energy.AtmosLatent      = 0.0;
    lake->energy.AtmosSensible    = 0.0;
    lake->energy.NetLongAtmos     = 0.0;
    lake->energy.NetShortAtmos    = 0.0;
    lake->energy.refreeze_energy  = 0.0;
    lake->energy.deltaH           = 0.0;
    lake->energy.Tsurf            = 0.0;
    lake->energy.lake_soil_heat_flux = 0.0;
    lake->energy.lake_soil_net_short = 0.0;

  }
 
  /**********************************************************************
   *  Calculate initial energy balance over ice.
   **********************************************************************/

  if ( fracprv >= FRACLIM ) {

    freezeflag = 0;         /* Calculation for ice. */
    Le = (677. - 0.07 * tair) * JOULESPCAL * GRAMSPKG; /* ice*/
    windi = ( wind * log((2. + soil_con->snow_rough) / soil_con->snow_rough) 
		/ log(wind_h/soil_con->snow_rough) );
    if ( windi < 1.0 ) windi = 1.0;
    lake->aero_resist = (log((2. + soil_con->snow_rough) / soil_con->snow_rough)
			   * log(wind_h/soil_con->snow_rough) / (von_K*von_K)) / windi;

    /* Calculate snow/ice temperature and change in ice thickness from surface melting. */
    ErrorFlag = ice_melt( wind_h+soil_con->snow_rough, lake->aero_resist, &(lake->aero_resist),
			  Le, lake_snow, lake, dt,  0.0, soil_con->snow_rough, 1.0, 
			  rainfall, snowfall,  windi, Tcutoff, tair, sw_ice, 
			  longin, air_density, pressure,  vpd,  vp, &lake->snowmlt, 
			  &lake->energy.advection, &lake->energy.deltaCC, 
			  &lake->energy.snow_flux, &Qei, &Qhi, &Qnet_ice, 
			  &temp_refreeze_energy, &LWneti, fracprv);
    if ( ErrorFlag == ERROR ) return (ERROR);

    lake->energy.refreeze_energy += temp_refreeze_energy*fracprv;
    *tempi = lake_snow->surf_temp;  

    /**********************************************************************
     *  Adjust temperatures of water column in ice fraction.
     **********************************************************************/

    if (lake->activenod > 0) {
      ErrorFlag = water_under_ice( freezeflag, sw_ice, wind, Ti, water_density, 
				   (double)soil_con->lat, lake->activenod, lake->dz, lake->surfdz,
				   Tcutoff, &qw, lake->surface, &temphi, water_cp, 
				   mixdepth, lake->hice, lake_snow->swq*RHO_W/RHOSNOW,
				   (double)dt, &energy_out_bottom_ice, soil_con->Zsum_node[1],
				   0.5*(lake->energy.kappa_node[0]+lake->energy.kappa_node[1]),
				   Tlakebot, T_node[1], &heat_cond_to_soil_ice);
      if ( ErrorFlag == ERROR ) return (ERROR);
    }
    else {
      temphi = -sumjoulb;
      heat_cond_to_soil_ice = 0;
      energy_out_bottom_ice = 0;
    }
	
    /**********************************************************************
     *  Calculate change in ice thickness and fraction
     *  within fraction that already has ice.
     **********************************************************************/
    /* Check to see if ice has already melted (from the top) in this time step. */
    if(lake->ice_water_eq > 0.0) {
      ErrorFlag = lakeice(tempi, Tcutoff, sw_ice, Ti[0], 
			  fracprv, dt, lake->energy.snow_flux, qw, 
			  &energy_ice_melt_bot, lake_snow->swq*RHO_W/RHOSNOW, 
			  lake->energy.deltaCC, rec, dmy, &qf,
			  ice_water_eq, lake->volume-*new_ice_water_eq, lake->surface[0]);
      if ( ErrorFlag == ERROR ) return (ERROR);
    }
    else {
      energy_ice_melt_bot = 0;
    }

    lake->energy.AtmosLatent += fracprv * Qei;
    lake->energy.advection       *= fracprv;
    lake->energy.AtmosSensible   += fracprv * Qhi;
    lake->energy.NetLongAtmos    += fracprv * LWneti;
    lake->energy.NetShortAtmos   += fracprv * sw_ice;
    lake->energy.deltaH          += fracprv * temphi;
    lake->energy.lake_soil_heat_flux  += (heat_cond_to_soil_ice)*fracprv;
    lake->energy.lake_soil_net_short  += (energy_out_bottom_ice)*fracprv;
//    lake->energy.refreeze_energy += energy_ice_melt_bot;
    lake->energy.refreeze_energy += energy_ice_melt_bot*fracprv;
    lake->energy.Tsurf += fracprv*lake_snow->surf_temp;

  }
  else {
    /* No Lake Ice Fraction */
    LWneti = 0;
    Qei    = 0.;
    Qhi    = 0;
    qf=0.0;
    temphi =0.0;
    lake->energy.refreeze_energy=0.0;
    if(fracprv > 0.0) {
      energy_ice_melt_bot = (lake->hice*RHOICE + (snowfall/1000.)*RHO_W)*Lf/(dt*SECPHOUR);
      *areai = 0.0;
      *hice = 0.0;
      *ice_water_eq = 0.0;

    }
    else {
      energy_ice_melt_bot = 0.0;
      *areai = 0.0;
      *hice = 0.0;
      *ice_water_eq = 0.0;
    }
  }

  /**********************************************************************
   * Average water temperature between ice-covered and ice-free columns.
   **********************************************************************/

  if(lake->activenod > 0) {
    // Average ice-covered and non-ice water columns.
    colavg( Tavg, T, Ti, fracprv, lake->density, lake->activenod,
            lake->dz, lake->surfdz);

    // Calculate depth average temperature of the lake
    *lake_tempavg = 0.0;
    for(i=0; i< lake->activenod; i++) {
      *lake_tempavg += Tavg[i]/lake->activenod;
    }
  }
  else
    *lake_tempavg = *tempi;

  /**********************************************************************
   * Compute T profile of soil under lake
   **********************************************************************/

  delta_t             = (double)dt * SECPHOUR;
  T2                  = T_node[Nnodes-1]; // soil column bottom temp
  Ts_old              = T_node[0]; // previous surface temperature
  T1_old              = T_node[1]; // previous first node temperature
  D1                  = soil_con->Zsum_node[1]-soil_con->Zsum_node[0];
  D2                  = soil_con->Zsum_node[2]-soil_con->Zsum_node[1];
  max_moist           = soil_con->max_moist[0] / (soil_con->depth[0]*1000.);

  // Compute T profile of soil under lake
  error = solve_surf_energy_bal(Tlakebot, rec, nrecs, (int)0, (int)2,
				 (int)0, (int)0, delta_t,
				 lake->energy.Cs[0], lake->energy.Cs[1],
				 D1, D2, T1_old, T2, Ts_old, lake->energy.T,
				 soil_con->bubble[0], soil_con->dp,
				 soil_con->expt[0], ice0,
				 lake->energy.kappa[0], lake->energy.kappa[1], 
				 max_moist, max_moist, &tmp_float, &tmp_double, (int)0, (int)0,
				 lake->energy.lake_soil_net_short, 
				 (double)0, (double)0, (double)0, (double)0, (double)0, 
				 (double)0, lake->energy.lake_soil_heat_flux, (double)0, (double)1, (double)1,
				 (double)0, (double)0, (double)0, (double)0, &tmp_double,
				 &tmp_double, &tmp_double, 
				 &tmp_double, &tmp_double, 
				 &tmp_double, &tmp_double, 
				 &tmp_double, &tmp_double, 
				 (double)0, (double)0, (double)0, (double)0, (double)0, (double)0, 
				 (double)0, (double)0, (double)0, (double)0, (double)0, &tmp_double, 
				 &tmp_double, &tmp_double, 
				 &tmp_double, &tmp_double, 
				 Nnodes, lake->energy.Cs_node, 
				 T_node, Tnew_node, 
				 Tnew_fbflag, Tnew_fbcount, 
				 soil_con->alpha, soil_con->beta, 
				 soil_con->bubble_node, soil_con->Zsum_node, 
				 soil_con->expt_node, soil_con->gamma, 
				 ice_node, lake->energy.kappa_node, 
				 soil_con->max_moist_node, moist_node, 
				 soil_con, &tmp_layer, &tmp_layer, 
				 &tmp_veg_var, &tmp_veg_var, 
				 (int)0, NOFLUX, options.EXP_TRANS, (int)0, FIRST_SOLN,
				 &tmp_double, &tmp_double, T1,
				 &tmp_deltaH, &tmp_fusion, 
				 &(lake->energy.grnd_flux), 
				 &(lake->energy.latent), 
				 &(lake->energy.latent_sub), 
				 &(lake->energy.sensible), 
				 &(lake->energy.snow_flux),
				 &(lake->energy.error));

  /* add the lake soil's deltaH and fusion to those of the lake water */
  lake->energy.deltaH += tmp_deltaH;
  lake->energy.fusion += tmp_fusion;

  Tnew_node[0] = Tlakebot;

  return(error);

}
