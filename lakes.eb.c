#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vicNl.h>

static char vcid[] = "$Id$";

int lakemain(atmos_data_struct  *atmos, 
	     lake_con_struct    lake_con,
	     double             snowprec,
	     double             rainprec,
	     soil_con_struct    *soil_con, 
	     veg_con_struct     *veg_con, 
#if EXCESS_ICE
	     int                SubsidenceUpdate,
	     double             total_meltwater,
#endif
	     int                dt, 
	     dist_prcp_struct    *prcp,
	     int                hour, 
	     int                rec, 
	     double             wind_h,
	     global_param_struct *gp,
	     dmy_struct          *dmy,
	     int                iveg,
	     int                band )
{
  /**********************************************************************
      lakes           Valentijn Pauwels             September 1998

  This subroutine sets up the boundary conditions for simulations with 
  the lake module and serves as the primary driver for the lake model.
  (Based on Hostetler and Bartlein, WRR, 26 (10), 2603-2612, 10/1990)
 
  Modifications:
  12/2000 Model was converted to C and incorporated into VIC.  LCB
  11-18-02 Modifications were made to improve handling of snow and ice
           and to make the lake algorithm interact with the wetland 
           algorithm.                                          LCB
  27-Aug-04 Added logic to handle blowing_flux and surface_flux.		TJB
  28-Sep-04 Added aero_resist_used to store the aerodynamic resistance
	    used in flux calculations.						TJB
  04-Oct-04 Merged with Laura Bowling's updated lake model code.		TJB
  23-Feb-05 Merged with Laura Bowling's second update to lake model code.	TJB
  2006-Jul-18 Changed sin(lat) to sin(fabs(lat)) in ks computation so that
	      southern hemisphere locations can be handled correctly.		TJB
  2006-Nov-07 Removed LAKE_MODEL option.					TJB
  2007-Apr-03 Modified to catch and return error flags from surface_fluxes
              subroutine.							KAC
  2007-Apr-24 Passes soil_con.Zsum_node to distribute_node_moisture_properties. JCA
  2007-Aug-09 Added features for EXCESS_ICE option.				JCA
  2007-Aug-16 Added ErrorFlag return value from initialize_prcp.		JCA
  2007-Oct-24 Changed get_sarea, get_volume, and get_depth to return exit
	      status so that errors can be trapped and communicated up the
	      chain of function calls.						KAC via TJB
  2007-Oct-24 Changed lakeice() to return exit status.				KAC via TJB
  2007-Nov-06 Major changes to lake model, including:				LCB via TJB
	        * New drainage parameterization: flow over broad-crested wier
	        * Snow on lake ice now has physics consistent with snow pack
	          on land
	        * Fixes for treatment of lake ice (previously, could grow in
	          excess of available lake water)
	        * ice_water_eq is now the state variable for ice volume; hice
	          represents the maximum ice thickness, not the average
	        * Lake surface area is now the larger of liquid water surface
	          area and lake ice surface area
	        * Fixes for some crashes in extreme cases and water balance
	          errors

  Parameters :

  soil_con.lat  	Latitude of the lake (degrees).
  soil_con.lng 	        Longitude of the lake (degrees).
  atmos                 Met forcings
  lake_energy.latent	Latent heat flux (W/m2).
  lake_energy.sensible	Sensible heat flux (W/m2).
  lake_con.eta_a	Decline of solar radiation input with depth (m-1).  
  lake_con.surface[numnod]	Area of the lake (m2).
  lake.temp[numnod]	Temperature of the lake water (K).
  lake.tempi	        Temperature of the lake ice (K).
  lake.hice           	Maximum height of the lake ice (m over ice area).
  lake.fraci	        Fractional coverage of ice (-).
  lake_snow.sdepth	Height of the snow on top of the lake (m over ice area).
  lake_snow.swe         SWE of snow on top of the lake (m).
  lake.numnod        	Number of nodes in the lake (-).
  global_param.dt	Time step size (hrs).

Note: lake->sarea and lake->areai are updated at the begining of each time
step and do not change during the time step.  They represent the lake surface
area at the beginning of the time step, and are the values that are used
DURING the time step, to determine the area of wetland and lake flux
calculations.  Therefore they must be used in put_data() when determining 
the grid cell average fluxes.  
**********************************************************************/

  extern option_struct   options;
  energy_bal_struct      lake_energy;
  snow_data_struct       lake_snow;
  snow_data_struct     **snow;
  lake_var_struct       *lake;
  cell_data_struct    ***cell;
  veg_var_struct      ***veg_var;

  double oldvolume;
  double oldsnow;
  int i, ErrorFlag;
  double lakefrac;
  double fraci;

  // DEBUG variables
  double V1, V2;
  double SWQ1, SWQ2, EVAP, PRECIP, DELTA, ROUT, RIN;
  double MOIST1, MOIST2, DELTAMOIST;
  double WDEW1, WDEW2, DDEW;
  double VAPOR;
  double DVOL, DSWE;

  veg_var = prcp->veg_var;
  cell    = prcp->cell;
  snow    = prcp->snow;
  lake    = &prcp->lake_var;

  /* Update sarea to equal new surface area from previous time step. */
  lake->sarea = lake->surface[0];	
  lake->areai = lake->new_ice_area;

  // Initialize DEBUG variables
  V1 = lake->volume*1000./lake_con.basin[0];
  SWQ1 = snow[iveg][band].swq*1000.;
  MOIST1= 0.0;
  WDEW1 = veg_var[WET][iveg][band].Wdew;
  for(i=0; i<options.Nlayer; i++) {
        MOIST1 += cell[WET][iveg][band].layer[i].moist;
  }

  // Initialize lake data structures 
  ErrorFlag = initialize_prcp(prcp, &lake_energy, &lake_snow, iveg, band, 
			      *soil_con, lake, rec, lake_con, &lakefrac, &fraci);
  if ( ErrorFlag == ERROR )
    // Return failure flag to main routine
    return ( ErrorFlag );

  /**********************************************************************
   * Solve the energy budget for the exposed land.
   **********************************************************************/	 
  ErrorFlag = wetland_energy(rec, atmos, prcp, dmy, gp, soil_con, 
#if EXCESS_ICE
			     SubsidenceUpdate,
#endif
			     iveg, band, lakefrac, lake_con, veg_con);

  if ( ErrorFlag == ERROR )
    // Return failure flag to main routine
    return ( ErrorFlag );

  /**********************************************************************
   * Solve the energy budget for the lake.
   **********************************************************************/
 
  oldvolume = lake->volume;
  oldsnow = lake_snow.swq;
 
  atmos->out_prec += (snowprec + rainprec) * lake_con.Cl[0] * lakefrac;
  ErrorFlag = solve_lake(snowprec, rainprec, atmos->air_temp[hour],
			 atmos->wind[hour], atmos->vp[hour] / 1000.,
			 atmos->shortwave[hour], atmos->longwave[hour],
			 atmos->vpd[hour] / 1000., atmos->pressure[hour] / 1000.,
			 atmos->density[hour], lake, lake_con, *soil_con,
			 dt, rec, &lake_energy, &lake_snow, wind_h, dmy[rec], fraci);
  if ( ErrorFlag == ERROR ) return (ERROR);

  /**********************************************************************
   * Solve the water budget for the lake.
   **********************************************************************/

  ErrorFlag = water_balance(lake, lake_con, dt, prcp, rec, iveg, band, lakefrac, *soil_con, 
#if EXCESS_ICE
		 SubsidenceUpdate, total_meltwater,
#endif		 
		 snowprec+rainprec, oldvolume, oldsnow-lake_snow.swq, lake_snow.vapor_flux, fraci);
  if ( ErrorFlag == ERROR ) return (ERROR);

  /**********************************************************************
   * Reallocate fluxes between lake and wetland fractions.
   **********************************************************************/

  update_prcp(prcp, &lake_energy, &lake_snow, lakefrac, lake,
	      lake_con, iveg, band, *soil_con);

  /**********************************************************************
   * DEBUG water balance for lake and wetland combined.
   **********************************************************************/

  V2 = lake->volume*1000./lake_con.basin[0];
  SWQ2= snow[iveg][band].swq*1000.;
  PRECIP = (rainprec+snowprec);
  ROUT = ((lake->runoff_out+lake->baseflow_out))*soil_con->cell_area/lake_con.basin[0];
  RIN = (( lake->runoff_in  + lake->baseflow_in ))*soil_con->cell_area/lake_con.basin[0];
  EVAP = (lakefrac*lake->evapw);
  MOIST2= 0.0;
  for(i=0; i<options.Nlayer; i++) {
    MOIST2 += cell[WET][iveg][band].layer[i].moist;
    EVAP += cell[WET][iveg][band].layer[i].evap;
  }
  EVAP += veg_var[WET][iveg][band].canopyevap;
  VAPOR = (snow[iveg][band].vapor_flux + snow[iveg][band].canopy_vapor_flux)*1000.;
  WDEW2 = veg_var[WET][iveg][band].Wdew;
  DDEW = WDEW2 - WDEW1;
  DSWE = SWQ2-SWQ1;
  DVOL = V2-V1;
  DELTAMOIST = (MOIST2-MOIST1)+DDEW + DSWE + DVOL;
  DELTA = DELTAMOIST - (RIN + PRECIP) + EVAP + VAPOR + ROUT;

//  if(fabs(DELTA) > SMALL) {
//    fprintf(stderr, "rec %d Error = %.5e (DeltaMoist=%f Dv=%f DSWE=%f) (Rin=%f, R=%f S=%f) (Rout= %f, vapor=%f evap=%f) ice=%.2f lake=%.2f\n",rec,DELTA,DELTAMOIST,DVOL,DSWE, RIN, (rainprec), (snowprec),ROUT, VAPOR, EVAP, fraci, lakefrac);
//    exit(0);
//  }

  return (0);

}

 /**********************************************************************
  * End of lake main.
  **********************************************************************/

int solve_lake(double             snow,
	       double             rain,
	       double             tair, 
	       double             wind, 
	       double             vp, 
	       double             shortin, 
	       double             longin, 
	       double             vpd, 
	       double             pressure, 
	       double             air_density, 
	       lake_var_struct   *lake, 
	       lake_con_struct    lake_con, 
	       soil_con_struct    soil_con, 
	       int                dt, 
	       int                rec, 
	       energy_bal_struct *lake_energy, 
	       snow_data_struct  *lake_snow, 
	       double             wind_h,
	       dmy_struct         dmy,
	       double             fracprv)
{

/**********************************************************************
* This subroutine solves the energy budget for open water bodies.
*
* Parameter description :
*
  lat  	Latitiude of the lake (degrees).
  lng 	Longitude of the lake (degrees).
  tair	Air temperature (C).
  vp	Air vapor pressure (kPa).
  pressure      	Air pressure (kPa).
  wind		Wind speed (m/s).
  longwave		Downwelling long wave radiation (W/m2).
  shortwave		Incoming short wave radiation (W/m2).
  prec	         Precipitation (mm).
  lake_energy.latent	 Latent heat flux (W/m2).
  lake_energy.sensible	 Sensible heat flux (W/m2).
  lake_con.eta_a	 Decline of solar radiation input with depth (m-1).  
  lake_con.surface[numnod]	Area of the lake (m2).
  lake.temp[numnod]	Temperature of the lake water (C).
  lake.tempi                Temperature of the lake ice (C).
  lake.hice      	         Height of the lake ice (m).
  lake.fraci	         Fractional coverage of ice (-).
  lake_snow.depth	         Height of the snow on top of the lake (m).
  lake.numnod        	Number of nodes in the lake (-).
  dt		         Time step size (hrs).

  Modifications:
  2006-Oct-16 Now set mixdepth=0 for case of complete ice cover; this
	      guarantees that it is initialized for all cases.  TJB
  2006-Nov-15 Convert swq and surf_water from mm over lake to mm over ice
	      fraction at beginning of function; this was needed to avoid
	      a water budget error since swq and surf_water were being
	      converted to mm over lake at end of the function.  TJB
  2007-Apr-21 Added initialization of energy_ice_melt_bot and qf for case
	      in which fracprv >= FRACLIM but hice == 0.0.		TJB
  2007-Apr-23 Added initialization of lake_energy->Tsurf.		TJB
  2007-Nov-06 Lake snow physics now consistent with land snow pack.
	      Lake ice now limited by available lake water.		LCB via TJB
  2008-Apr-21 Corrected omissions in previous mods.  Also added code to
	      handle previously-unhandled case in deltaH computation.	TJB
**********************************************************************/

  double LWnetw,LWneti;
  double sw_water, sw_ice;
  double T[MAX_LAKE_NODES];   /* temp of the water column, open fraction. */
  double Tnew[MAX_LAKE_NODES];  
  double Ti[MAX_LAKE_NODES];   /* temp of the water column, ice fraction. */
  double water_density[MAX_LAKE_NODES], water_cp[MAX_LAKE_NODES];
  float albi, albw;
  double Ts, tempalbs;
  double Tcutoff, Tcutk;  /* Lake freezing temperature (K). */
  double Qhw, Qhi;
  double Qgw, Qgi;
  double Qew, Qei;
  double eflux, eadd;
  double Tki, Tkw;     /* Surface temp. of ice and water in Kelvin. */
  double qbot, qw;
  int freezeflag;
  int mixdepth;
  double new_ice_area; /* Ice area formed by freezing in the open water portion. */
  int k, i, ErrorFlag;
  double tw1, tw2, r1, r2, rtu, rnet;
  double Ls, Le;
  double sumjoula, sumjoulb, sumjouli;
  double temp;
  double water_energy;
  double temphw, temphi;
  double SW_under_ice;
  int last_snow;
  double cc;
  double ice_energy;
  double error_over_water;
  double sw_underice_visible;
  double sw_underice_nir;
  double T_lower, T_upper;
  double Tsurf;
  double new_ice_height;
  double rainfall, snowfall;
  double windw, windi;
  double energy_ice_formation;
  double energy_ice_melt_bot;
  double melt, deltaCC_ice, Qnet_ice;
  double joule_intermediate;
  double energy_out_bottom;
  double energy_out_bottom_ice;    
  double sw_out;
  int index;
  double tempdepth, in;
  double qf;
  double inputs, outputs, internal, phasechange;
  double new_ice_water_eq;
  
  /**********************************************************************
   * 1. Initialize variables.
   **********************************************************************/
  
  lake_energy->advection=0.0;
  lake_energy->deltaCC = 0.0;
  lake_energy->grnd_flux = 0.0;
  lake_energy->snow_flux = 0.0;
  lake->snowmlt = 0.0;
  qbot=qw=0.0;
  new_ice_height = new_ice_area = new_ice_water_eq = 0.0;
  lake->evapw=0.0;
  energy_ice_formation = 0.0;
  energy_out_bottom = energy_out_bottom_ice = 0.0;
  lake_snow->vapor_flux=0.0;
  lake_energy->Tsurf = lake->temp[0];

  if(lake->activenod > 0 || lake->areai > 0.0) {

    /* --------------------------------------------------------------------
     * Calculate the water freezing point.
     * -------------------------------------------------------------------- */

    rhoinit(&Tcutoff, pressure);          /* degrees C */

    /* --------------------------------------------------------------------
     * Initialize liquid water and ice temperature profiles to be used
     * later in the module.
     * -------------------------------------------------------------------- */

    for ( k = 0; k < lake->activenod; k++ ) {
      T[k] = lake->temp[k];
      Ti[k] = lake->temp[k];
      water_density[k] = calc_density(T[k]);
      water_cp[k] = specheat(T[k]);
    }

    energycalc( lake->temp, &sumjoulb, lake->activenod, lake->dz, 
		lake->surfdz,lake->surface, water_cp, water_density);

    /**********************************************************************
     * 2. Calculate added precipitation and total snow height.
     **********************************************************************/

    if ( fracprv >= 1.0) {  /* areai is relevant */
      rainfall = rain;
      snowfall = snow;

      // If there is no snow, add the rain directly to the lake.
      if(lake_snow->swq <=0.0 && rainfall > 0.0) {
        lake->volume += (rainfall/1000.)*lake->areai;
        rainfall = 0.0;
      }
    }
    else if ( fracprv > FRACLIM && fracprv < 1.0 ) { /*sarea is relevant */
      rainfall = rain;
      snowfall = snow;

      /* Precip over open water directly increases lake volume. */
      lake->volume += ( (snow/1000. + rain/1000.)*(1-fracprv)
                        * lake->sarea );

      if(lake_snow->swq <= 0.0 && rainfall > 0.0) {
        lake->volume += (rain/1000.)*fracprv*lake->sarea;
        rainfall = 0.0; /* Because do not want it added to snow->surf_water */
      }
    }
    else {
      lake->volume += ((rain+snow)/1000.)*lake->sarea;
      rainfall = 0.0;
      snowfall = 0.0;
    }

    /**********************************************************************
     * 3. Calculate incoming solar radiation over water and ice.
     **********************************************************************/

    /* --------------------------------------------------------------------
     * Calculate the albedo of the lake for ice, snow and liquid water.
     * To be consistent with VIC snow model, tempalbs = NEW_SNOW_ALB if snow
     * is falling, but aging routing does not get reset.
     * -------------------------------------------------------------------- */

    alblake(Tcutoff, tair, &lake->SAlbedo, &tempalbs, &albi, &albw, snowfall,
	    lake_snow->coldcontent, dt, &lake_snow->last_snow, 
	    lake_snow->swq, lake_snow->depth, &lake_snow->MELTING, dmy.day_in_year );

    /* --------------------------------------------------------------------
     * Calculate the incoming solar radiaton for both the ice fraction
     * and the liquid water fraction of the lake.
     * -------------------------------------------------------------------- */

    if (lake_snow->swq > SNOWCRIT *RHOSNOW/RHO_W) {
      sw_ice = shortin * (1.-tempalbs);
      lake_energy->AlbedoLake = (fracprv) * tempalbs + (1. - fracprv) * albw;
    }
    else if (lake_snow->swq > 0. && lake_snow->swq <= SNOWCRIT*RHOSNOW/RHO_W) {
      sw_ice = shortin * (1.- (albi+tempalbs)/2.);
      lake_energy->AlbedoLake = (fracprv) * (albi+tempalbs)/2. + (1. - fracprv) * albw;
    }
    else if (fracprv > 0. && lake_snow->swq <= 0.) {
      sw_ice = shortin * (1.-albi);
      lake_energy->AlbedoLake = (fracprv) * albi + (1. - fracprv) * albw;
    }
    else /* lake->fraci = 0 */ {
      sw_ice=0.0;
      lake_energy->AlbedoLake = albw;
    }

    sw_water = shortin * (1.-albw);

    /**********************************************************************
     * 4. Calculate initial energy balance over water.
     **********************************************************************/

    if( (1.-fracprv) > SMALL && lake->activenod > 0) {
      freezeflag=1;         /* Calculation for water, not ice. */
      windw = wind*log((2. + ZWATER)/ZWATER)/log(wind_h/ZWATER);

      ErrorFlag = water_energy_balance( lake->activenod, lake->surface, &lake->evapw, 
					dt, freezeflag, lake->dz, lake->surfdz,
					(double)soil_con.lat, Tcutoff, tair, windw, 
					pressure, vp, air_density, longin, sw_water, 
					sumjoulb, wind_h, &Qhw, &Qew, &LWnetw, T, 
					water_density, &lake_energy->deltaH, 
					&energy_ice_formation, fracprv, 
					&new_ice_area, water_cp, &new_ice_height,
					&energy_out_bottom, &new_ice_water_eq,
					lake->volume-lake->ice_water_eq);
      if ( ErrorFlag == ERROR ) return (ERROR);
 
      /* --------------------------------------------------------------------
       * Do the convective mixing of the lake water.
       * -------------------------------------------------------------------- */

      mixdepth = 0;        /* Set to zero for this time step. */
      tracer_mixer( T, &mixdepth, freezeflag, lake->surface, 
		    lake->activenod, lake->dz, lake->surfdz, water_cp );

    }          /* End of water fraction calculations. */
    else {
      // ice covers 100% of lake, reset open water fluxes
      mixdepth = 0;
      LWnetw = 0;
      Qew    = 0;
      Qhw    = 0;
    }
 
    /**********************************************************************
     *  6. Calculate initial energy balance over ice.
     **********************************************************************/

    if ( fracprv >= FRACLIM ) {

      freezeflag = 0;         /* Calculation for ice. */
      Le = (677. - 0.07 * tair) * JOULESPCAL * GRAMSPKG; /* ice*/
      windi = ( wind * log((2. + soil_con.snow_rough) / soil_con.snow_rough) 
		/ log(wind_h/soil_con.snow_rough) );
      if ( windi < 1.0 ) windi = 1.0;
      lake->aero_resist = (log((2. + soil_con.snow_rough) / soil_con.snow_rough)
			   * log(wind_h/soil_con.snow_rough) / (von_K*von_K)) / windi;

      /* Calculate snow/ice temperature and change in ice thickness from 
         surface melting. */
      ErrorFlag = ice_melt( wind_h+soil_con.snow_rough, lake->aero_resist, &(lake->aero_resist_used),
			    Le, lake_snow, lake, dt,  0.0, soil_con.snow_rough, 1.0, 
			    rainfall, snowfall,  windi, Tcutoff, tair, sw_ice, 
			    longin, air_density, pressure,  vpd,  vp, &lake->snowmlt, 
			    &lake_energy->advection, &lake_energy->deltaCC, 
			    &lake_energy->snow_flux, &Qei, &Qhi, &Qnet_ice, 
			    &lake_energy->refreeze_energy, &LWneti, fracprv);
      if ( ErrorFlag == ERROR ) return (ERROR);

      lake->tempi = lake_snow->surf_temp;  

      /**********************************************************************
       *  7. Adjust temperatures of water column in ice fraction.
       **********************************************************************/

      /* --------------------------------------------------------------------
       * Calculate inputs to temp_area..
       * -------------------------------------------------------------------- */

      if (lake->activenod > 0) {
        ErrorFlag = water_under_ice( freezeflag, sw_ice, wind, Ti, water_density, 
				     (double)soil_con.lat, lake->activenod, lake->dz, lake->surfdz,
				     Tcutoff, &qw, lake->surface, &temphi, water_cp, 
				     mixdepth, lake->hice, lake_snow->swq*RHO_W/RHOSNOW,
				     (double)dt, &energy_out_bottom_ice);     
        if ( ErrorFlag == ERROR ) return (ERROR);
      }
      else
        temphi = -sumjoulb;
	
      /**********************************************************************
       *   8.  Calculate change in ice thickness and fraction
       *    within fraction that already has ice.
       **********************************************************************/
      /* Check to see if ice has already melted (from the top) in this time step. */
      if(lake->ice_water_eq > 0.0) {
        ErrorFlag = lakeice(&lake->tempi, Tcutoff, sw_ice, Ti[0], 
			    fracprv, dt, lake_energy->snow_flux, qw, 
			    &energy_ice_melt_bot, lake_snow->swq*RHO_W/RHOSNOW, 
			    lake_energy->deltaCC, rec, dmy, &qf,
			    &lake->ice_water_eq, lake->volume-new_ice_water_eq, lake->surface[0]);
        if ( ErrorFlag == ERROR ) return (ERROR);
      }
    }
    else {
      /* No Lake Ice Fraction */
      LWneti = 0;
      Qei    = 0.;
      Qhi    = 0;
      qf=0.0;
      temphi =0.0;
      lake_energy->refreeze_energy=0.0;
      if(fracprv > 0.0) {
        energy_ice_melt_bot = (lake->hice*RHOICE + (snowfall/1000.)*RHO_W)*Lf/(dt*SECPHOUR);
        lake->areai = 0.0;
        lake->hice = 0.0;
        lake->ice_water_eq = 0.0;

      }
      else {
        energy_ice_melt_bot = 0.0;
        lake->areai = 0.0;
        lake->hice = 0.0;
        lake->ice_water_eq = 0.0;
      }
    }

    /**********************************************************************
     * 9. Average water temperature.
     **********************************************************************/

    if(lake->activenod > 0) {
      // Average ice-covered and non-ice water columns.
      colavg( lake->temp, T, Ti, fracprv,lake->density, lake->activenod,
              lake->dz, lake->surfdz);

      // Calculate depth average temperature of the lake
      lake->tempavg = 0.0;
      for(i=0; i< lake->activenod; i++) {
        lake->tempavg += lake->temp[i]/lake->activenod;
      }
    }
    else
      lake->tempavg = -99;

    /**********************************************************************
     * 10. Calculate the final water heat content and energy balance.
     **********************************************************************/

    /* Incoming energy. */ 
    inputs = (sw_ice + LWneti + lake_energy->advection + Qhi + Qei);
    outputs = energy_out_bottom_ice;
    internal = temphi;
    phasechange = -1*(lake_energy->refreeze_energy)-1.*energy_ice_melt_bot;

    lake_energy->error = inputs - outputs - internal - phasechange;
    

    lake_energy->AtmosLatent      = ( 1. - fracprv ) * Qew + fracprv * Qei; 
    lake_energy->advection       *= fracprv;     
    lake_energy->AtmosSensible    = ( 1. - fracprv ) * Qhw + fracprv * Qhi;
    lake_energy->NetLongAtmos     = ( ( 1. - fracprv ) * LWnetw + fracprv 
					* LWneti );
    lake_energy->NetShortAtmos    = ( ( 1. - fracprv ) * sw_water 
					+ fracprv * sw_ice );

    lake_energy->refreeze_energy += energy_ice_melt_bot;
    lake_energy->refreeze_energy *= fracprv;
    lake_energy->refreeze_energy += energy_ice_formation*(1.-fracprv);

    lake_energy->snow_flux       = 0.0;
    lake_energy->deltaH          *= ( 1. - fracprv );
    lake_energy->deltaH          += fracprv*temphi;
    lake_energy->deltaH          *= -1.;
    lake_energy->deltaCC         = 0.0;
    lake_energy->grnd_flux        = -1.*(energy_out_bottom*(1.-fracprv) +
					   energy_out_bottom_ice*fracprv);

    lake_energy->Tsurf = fracprv*lake_snow->surf_temp + (1. - fracprv)*T[0];

    lake_energy->error            = ( lake_energy->NetShortAtmos 
				      + lake_energy->NetLongAtmos 
				      + lake_energy->AtmosSensible 
				      + lake_energy->AtmosLatent 		   
				      + lake_energy->deltaH 
				      + lake_energy->grnd_flux
				      + lake_energy->refreeze_energy 
				      + lake_energy->advection );

    temphw=temphi=0.0;

    /**********************************************************************
     * 11. Final accounting for passing variables back to VIC.
     **********************************************************************/

    /* Adjust water and ice evaporation to be representative of 
       entire lake. */
    lake->evapw           *= ( (1. - fracprv ) * dt * SECPHOUR ); // in mm
    lake_snow->vapor_flux *= fracprv; // in meters 
    lake_snow->blowing_flux *= fracprv; // in meters 
    lake_snow->surface_flux *= fracprv; // in meters 	

    // Adjust snow variables to represent average depth over lakefraction, based on
    // ice fraction at the beginning of the time step. Converted to average over
    // wetland land tile in update_prcp.
    lake_snow->swq        *= fracprv;
    lake_snow->surf_water *= fracprv;
    lake_snow->pack_water *= fracprv;
    lake_snow->depth = lake_snow->swq * RHO_W / RHOSNOW;

    /* Update ice area to include new ice growth in water fraction. */

    lake->new_ice_area  = lake->areai;

    if(new_ice_area > 0.0 ) {
      lake->new_ice_area += new_ice_area;
      lake->ice_water_eq += new_ice_water_eq;
    }

    if(lake->ice_water_eq > 0.0) {
      if(fabs(fracprv-1.0) > FRACLIM ) {
        lake->hice = (lake->ice_water_eq/lake->new_ice_area)*RHO_W/RHOICE;
      }
      else {
        ErrorFlag = ice_depth(lake_con, lake->volume, lake->ice_water_eq, &(lake->hice));
        if (ErrorFlag == ERROR) return(ERROR);
      }
    }
    else
      lake->hice = 0.0;

    /* Change area of ice-covered fraction if ice has thinned. */
    if(lake->hice <= 0.0) {
      lake->new_ice_area = 0.0;
      lake->hice = 0.0;
    }
    else if(lake->hice < FRACMIN) {
      lake->new_ice_area = (lake->new_ice_area * lake->hice) / FRACMIN;
      lake->hice = FRACMIN;
    }
  }

  return (0);

}   /* End of solve_lake function. */

void latsens (double Tsurf, double Tcutk, double hice, double tair, double wind,
	      double pressure, double vp, double air_density, 
	      double *evap, double *qsen, double wind_h)
{
/**********************************************************************
 * Calculate the partitioning of the energy balance into latent and
 * sensible heat fluxes.
 *
 * Parameters :
 *
 * Tsurf	Lake surface temperature (K).
 * Tcutk        Melting temperature (K)
 * hice		Ice height (m).
 * tair		Air temperature (C).
 * wind		Wind speed (m/s).
 * pressure	Air pressure (kPa).
 * density	Air density (kg/m3).
 * vp           Vapor pressure (kPa). 
 * evap		Evaporation rate (m/s).
 * qsen		Sensible heat flux (W/m2).
 * eog		Vapor pressure at lake level (kPa).
 **********************************************************************/

    float dragcoeff;
    double eog,delT;
    double qair, qlake;  /* Specific humidity of the atmosphere and lake, respectively. */
    double delq;		/* Difference in absolute humidity between the lake */
			/* surface and higher up (-). */

/**********************************************************************
 * Calculate the drag coefficient.
 **********************************************************************/

    if(hice > 0.)
      dragcoeff = lkdrag (Tsurf, tair + KELVIN, wind, ZSNOW, wind_h);
    else
      dragcoeff = lkdrag (Tsurf, tair + KELVIN, wind, ZWATER, wind_h);
       
/**********************************************************************
 * Determine the coefficients to be used in the calculation of
 * the vapor pressure at lake level depending on whether the lake is 
 * covered with ice or not.
 **********************************************************************/

	if ( (hice <= 0.) && (Tsurf > Tcutk) ) {

/* --------------------------------------------------------------------
 * Lake is not covered with ice (Handbook of Hydrology eq. 4.2.2).
 * eog in kPa.
 * -------------------------------------------------------------------- */
	  eog=.611*exp(17.269*(Tsurf-KELVIN)/(Tsurf+237.3-KELVIN));
	}
       
	else {

/* --------------------------------------------------------------------
 * Lake is covered with ice.
 * -------------------------------------------------------------------- */
       
	  eog=.611*exp(21.874*(Tsurf-KELVIN)/(Tsurf-7.66));
	}


/**********************************************************************
 * Calculate the specific humidity at lake level and the difference
 * between this humidity and the air humidity at measurement height.
 **********************************************************************/

      qlake=0.622*(eog/(pressure-0.378*eog));
      qair = 0.622*(vp/(pressure - 0.378*vp));
      delq = qair - qlake;
 
      /**********************************************************************
       * Calculate the evaporation rate.  Eq. 4 in Hostetler (1991), 
       * after Brutsaert (1982).  evap in mm/s
       **********************************************************************/

      *evap=-1*dragcoeff*wind*air_density*delq;

/**********************************************************************
 * Calculate the difference in lake surface temperature and the
 * air temperature at measurement height and consequently the
 * sensible heat flux.  Hostetler (1991) eq 5, in W/m^2.
 **********************************************************************/

      delT=tair+KELVIN-Tsurf;
      *qsen=dragcoeff*wind*air_density*Cp;
      *qsen=*qsen*delT;
}


void alblake (double  Tcutoff, 
	      double  Tair, 
	      double *snowalbedo, 
	      double *albs, 
	      float  *albi, 
	      float  *albw, 
	      double  newsnow, 
	      double  coldcontent, 
	      int     dt, 
	      int    *last_snow, 
	      double  swq,
	      double  depth,
	      char   *MELTING,
	      int     day_in_year)
{
/**********************************************************************
 * Calculate the albedo of snow, ice and water of the lake.
 *
 * Parameters :
 *
 * tair		Air temperature (C).
 * albs		Snow albedo (-).
 * albi		Ice albedo (-).
 * albw		Water albedo (-).
 *
 * Modifications:
 * 2007-Nov-06 Lake snow physics now consistent with land snow pack.	LCB via TJB
 * 2008-Apr-21 Corrected some omissions from the previous mods.		LCB via TJB
 **********************************************************************/
  double albgl, albgs;

  if( (Tair-Tcutoff) > 0.0) {
    if( (Tair-Tcutoff) < 20.) {
      albgl=0.4 - 0.011*(Tair-Tcutoff);
      albgs=0.6 - 0.0245*(Tair-Tcutoff);
    }
    else {	
      albgl=0.4 - 0.011*20.;
      albgs = 0.6 - 0.0245*20.;
    }
  }
  else {
    albgl=0.4;
    albgs=0.6;
  }

  *albi=0.5*albgs + 0.5*albgl; 

  // update number of days since last significant snowfall
  if ( newsnow > TraceSnow )
    *last_snow = 1;
  else if ( swq == 0 )
    *last_snow = 0;
  else
    *last_snow +=1;

  /** Record if snowpack is melting this time step **/
  if (swq > 0.0) {
    if ( coldcontent >= 0 && day_in_year > 60 // ~ March 1
       && day_in_year < 273 // ~ October 1
       ) {
      *MELTING = TRUE;
    }
    else {
      *MELTING = FALSE;
    }
  }
  else {
    *MELTING = FALSE;
  }

  if ( *MELTING && newsnow > TraceSnow )
    *MELTING = FALSE;

  // compute snow surface albedo
  if(swq > 0.0)
// Uncomment this if we take the update to snow_albedo() for Sun et al. 1999
//    *snowalbedo = snow_albedo(newsnow, swq, depth, *snowalbedo, coldcontent, dt, *last_snow, *MELTING);
    *snowalbedo = snow_albedo(newsnow, swq, coldcontent, dt, *last_snow, *MELTING);
  else if(swq == 0.0 && newsnow > 0.0)
    *snowalbedo = NEW_SNOW_ALB;
  else
    *snowalbedo = 0.0;

  if(newsnow > 0.0)
    *albs = NEW_SNOW_ALB;
  else
    *albs = *snowalbedo;

  *albw = 0.15;

}


void colavg (double *finaltemp, double *T,double *Ti,float lakeprv,double *density, int numnod, double dz, double surfdz)
  {

/**********************************************************************
 * Calculate the temperature of the lake water at the different node
 * levels and calculate the water density.
 *
 * Parameters :
 * 
 * finaltemp    Averaged water column temperature (C).
 * T		Water temperature for water fraction(C).
 * Ti		Water temperature for ice fraction(C).
 * lakeprv	Fraction of lake covered with ice (-).
 * density  	Water density at each node (kg/m3).
 * numnod	Number of nodes in the lake (-).
 * dz           Thickness of each water layer (m). 
 * surfdz       Thickness of the surface water layer.
 **********************************************************************/
    int j;
    double water_densityw, water_densityi;
    double temp;
    float z;
  
      for(j=0; j<numnod; j++)
	{

/* --------------------------------------------------------------------
 * Calculate the densities of the ice and water fractions.
 * -------------------------------------------------------------------- */

         water_densityw = calc_density(T[j]);
         water_densityi = calc_density(Ti[j]);
         water_densityw = water_densityw+1000.;
	 water_densityi = water_densityi+1000.; 

/* --------------------------------------------------------------------
 * Calculate the depth differences.
 * -------------------------------------------------------------------- */

         z=dz;
         if (j == 0)
	   z=surfdz;

/* --------------------------------------------------------------------
 * Calculate the lake (water) temperature as a weight of ice and water
 * temperatures.
 * -------------------------------------------------------------------- */

	 /*       temp=((1.-lakeprv)*T[j]*z*water_densityw*sh_water+
              lakeprv*Ti[j]*z*water_densityi*sh_ice)/
              ((z*water_densityw*sh_water+z*water_densityi*sh_ice)*0.5);
	 */
	 temp=((1.-lakeprv)*T[j]*z*water_densityw+
	       lakeprv*Ti[j]*z*water_densityi)/
	   ((1.-lakeprv)*z*water_densityw+lakeprv*z*water_densityi);
	 
   finaltemp[j]=temp;
	
/* --------------------------------------------------------------------
 * Recalculate the water density at the different nodes.
 * -------------------------------------------------------------------- */

	 density[j] = calc_density(finaltemp[j]);
	}
  }

double calc_density(double ts)
  {

/**********************************************************************
 * Calculate the water density.
 *
 * Parameters :
 *
 * ts		Temperature at the node (C).
 *
 * Returns:
 *
 * rhostps	Water density at the node (kg/m3).
 **********************************************************************/

      double rhow,t,rhostps;

      t=ts;

/* --------------------------------------------------------------------
 * Calculate the water density as a function of temperature.
 * -------------------------------------------------------------------- */

      rhow=999.842594 + .06793952*t - .009095290*t*t +
	.0001001685*t*t*t - .000001120083*t*t*t*t + .000000006536332*t*t*t*t*t;

     

/* --------------------------------------------------------------------
 * Return the difference between 1000 kg/m3 and the calculated
 * density, not the actual value.
 * -------------------------------------------------------------------- */

      rhostps=rhow-1000.;
      return (rhostps);
  }

void eddy (int freezeflag, double wind, double *T, double *water_density, 
	   double *de, double lat, int numnod, double dz, double surfdz)
   {

/**********************************************************************
 * Calculate the eddy diffusivity.
 *
 * Parameters :
 *
 * freezeflag	= 1 if liquid water, 0 if ice.
 * wind		Wind speed at 2 meter (m/s).
 * T		Water temperature at the different nodes (C).
 * density		Water density at the different nodes (C).
 * de		Eddy diffusivity at each node (m2/d).
 * lat	         Latitude of the pixel (degrees).
 * numnod	         Number of nodes in the lake (-).
 **********************************************************************/

      double ks, N2, ws, radmax, Po;
      double dpdz, rad;
      double zhalf[MAX_LAKE_NODES];  
      int k;
      double z;
      double Ri; /* Richardson's number. */      

/**********************************************************************
 * Calculate the density at all nodes and the distance between all nodes.
 **********************************************************************/

      for(k=0; k< numnod; k++)
	{
         zhalf[k] = dz;
	}

/**********************************************************************
 * Calculate the distance between the first node and the surface node.
 **********************************************************************/

      zhalf[0] = (surfdz + dz) * 0.5;

/**********************************************************************
 * If there is ice only molecular diffusivity is taken in account,
 * no eddy diffusivity.
 **********************************************************************/

      if (freezeflag != 1) {
	for(k=0; k< numnod; k++)
	  {
            de[k]=DM;
	  }
      }

      /**********************************************************************
       * Avoid too low wind speeds for computational stability.
       **********************************************************************/
     
      else {
	
	if(wind < 1.0)
	  wind = 1.0;

      /**********************************************************************
       *Determine the latitudinaly dependent parameter of the Ekman
       * profile, the surface value of the friction velocity and the neutral
       * value of the Prandtl number (Hostetler and Bartlein eq. 6 and 7).
       **********************************************************************/
      
	ks=6.6*pow(sin((double)fabs(lat)*PI/180.),0.5)*pow(wind,-1.84);
	ws=0.0012*wind;
	Po=1.0;

      /**********************************************************************
       * Determine the maximum value of the second term in the calculation
       * of the Richardson number for computational stability.
       **********************************************************************/

	radmax=4.e4;
      
	for(k= 0; k<numnod-1; k++) {

	/**********************************************************************
         * Calculate the eddy diffusivity for each node.
	 **********************************************************************/

	/* --------------------------------------------------------------------
	 * Calculate the derivative of density with depth, the Brunt-Vaisala
	 * frequency and the depth of the node (Hostetler and Bartlein eq. 9).
	 * -------------------------------------------------------------------- */
	
	  dpdz=(water_density[k+1]-water_density[k])/zhalf[k];
	  N2=(dpdz/(1.e3+water_density[k]))*9.8;   
	  z=surfdz+((float)k)*dz;	 
	
	/* --------------------------------------------------------------------
	 * Calculate the second term in the calculation of the Richardson
	 * number, make sure this number does not get too large for
	 * computational stability.  Also make sure this term does not go
	 * below 1 to avoid negative Richardson numbers (Hostetler and Bartlein eq. 8).
	 * -------------------------------------------------------------------- */
 
	  if ((z*exp(ks*z)/ws) > 1.e8) 
	    rad = radmax;
	
	  else 
	    {
	      rad=1.+40.*N2*(von_K*z)*(von_K*z)/(ws*ws*exp(-2.*ks*z));
	      
	      if (rad > radmax) 
		rad=radmax;
	    
	      if (rad < 1.0) 
		rad=1.0;
	    }

	/* --------------------------------------------------------------------
	 * Calculate the Richardson number and the eddy diffusivity.
	 * (Hostetler and Bartlein eq. 8 and 5.
	 * -------------------------------------------------------------------- */

	  Ri=(-1.0+sqrt(rad))/20.0;
	  de[k]=DM+(von_K*ws*z/Po)*exp(-ks*z)
	    /(1.0+37.0*Ri*Ri);
	}

      /* --------------------------------------------------------------------
       * The eddy diffusivity of the last node is assumed to equal the
       * eddy diffusivity of the second last node.
       * -------------------------------------------------------------------- */
      
	de[numnod-1]=de[numnod-2];
      }
      
   }

void iceform (double *qfusion, 
	      double *T, 
	      double  Tcutoff, 
	      double  fracprv,
	      double *areaadd, 
	      int     numnod, 
	      int     dt,  
	      double  dz, 
	      double  surfdz,	
	      double *cp, 
	      double *surface, 
	      double *new_ice_height, 
	      double *water_density,
	      double *new_ice_water_eq,
	      double  lvolume)
   {

/**********************************************************************
 * Calculate the form of new ice in the lake as long as the fractional
 * coverage of the ice is not 1 yet.
 *
 * Parameters :
 *
 * T		Water temperatures for all the nodes (C).
 * Tcutoff	Temperature at which water freezes (C).
 * fracprv	Fractional coverage of ice before this calculation (-).
 * fracadd	Added fractional ice coverage (-).
 * lake.fraci	New fractional ice coverage (-).
 * qfusion      Heat flux absorbed into the ice (W/m2).
 * hice		Ice height (m).
 * dt		Time step (s).
 * numnod	Number of nodes in the lake (-).
 * cp           Specific heat (J/Kg K) 
 *
 * Modifications:
 * 2007-Nov-06 Ice formation is now limited by available liquid water.
 *             ice_water_eq is now the state variable for lake ice water
 *             storage (used to be hice).				LCB via TJB
 **********************************************************************/
     double sum, extra;
     int j;
     double xfrac, di;

/**********************************************************************
 * Calculate the newly added ice to the ice layer.
 **********************************************************************/
     *qfusion = 0.0;
     sum = 0.0;
     *new_ice_water_eq = 0.0;

      for(j=0; j<numnod;j++) {
	  if (T[j]  < Tcutoff) {

/* --------------------------------------------------------------------
 * Calculate the ice growth (extra in J).
 * -------------------------------------------------------------------- */
	    
            if (j == 0) 
	      extra = ( (Tcutoff-T[j]) * surfdz * RHO_W
			* cp[j]  * (1.0-fracprv) * (surface[j]+surface[j+1])/2.);
	    else if (j == numnod-1)
	      extra = ( (Tcutoff-T[j]) * dz * RHO_W
			* cp[j]  * (1.0-fracprv) * surface[j]);
	    else
	      extra = ( (Tcutoff-T[j]) * dz * RHO_W
		        * cp[j]  * (1.0-fracprv) * (surface[j]+surface[j+1])/2.);
/* --------------------------------------------------------------------
 * If ice, the temperature of the water is the freezing temperature.
 * -------------------------------------------------------------------- */

            T[j]=Tcutoff;
	 
            sum+=extra;
	    }
	}

/**********************************************************************
 * Calculate the heat flux absorbed into the ice (in W/m2) and the 
 * thickness of the new ice.
 **********************************************************************/

      *new_ice_water_eq = sum/(RHO_W*Lf);

      if(lvolume > *new_ice_water_eq) {

        *qfusion=(sum/(dt*SECPHOUR*surface[0]*(1.0-fracprv))); /*W/m2*/

        di=sum/(Lf*RHOICE);    /* m^3 of ice formed */
      }
      else if (lvolume > 0.0 ){
        *new_ice_water_eq = lvolume;

        di = *new_ice_water_eq*RHO_W/RHOICE;

        // NEED TO CHANGE ICE TEMPERATURE TO ACCOUNT FOR EXTRA qfusion
        *qfusion=(*new_ice_water_eq*RHO_W/RHOICE)/(dt*SECPHOUR*surface[0]*(1.0-fracprv)); /*W/m2*/
        //      *deltaCC=(sum - *new_ice_water_eq*RHO_W/RHOICE)/(dt*SECPHOUR*surface[0]*(1.0-fracprv)); /*W/m2*/
      }
      else {
        *new_ice_water_eq = 0.0;

        di = 0.0;

        *qfusion=0.0; /*W/m2*/
      }

      /**********************************************************************
       * Calculate the added fractional coverage of ice, make sure the
       * total fractional coverage does not exceed 1.
       **********************************************************************/

      *new_ice_height = FRACMIN;
      *areaadd=di/(FRACMIN);

      if ( *areaadd > (1.0 - fracprv)*surface[0]) {
        *new_ice_height = di/((1.-fracprv)*surface[0]);
        *areaadd=(1.-fracprv)*surface[0];
      }

   }

void icerad (double  sw,
	     double  hi,
	     double  hs,
	     double *avgcond,
	     double *SWnet,
	     double *SW_under_ice)
{

/**********************************************************************
 * Calculate the radiation balance over ice.
 *
 * Paramterers :
 *
 * sw		Net solar radiation at top of snowpack (W/m2).
 * hi		Ice depth (m).
 * hs		Snow depth (m).
 * avgcond	Thermal conductivity of the ice and snow pack
 *		combined (W/mK). 
 * SWnet	Net short wave radiation at the top of the lake ice.
 * SW_under_ice	Incoming short wave radiation at the bottom of
 *		the snow-ice layer.
 **********************************************************************/
 
  double a, b, c, d;

/**********************************************************************
 * Calculate the thermal conductivity of the combined snow and ice
 * layer.
 **********************************************************************/

  *avgcond=(hs*CONDI+hi*CONDS)/(CONDI*CONDS);

/**********************************************************************
 * Calculate the incoming radiation at different levels in the
 * ice-snow layer.
 **********************************************************************/

/* --------------------------------------------------------------------
 * Calculation of constants. Patterson and Hamblin eq, 7
 * -------------------------------------------------------------------- */

  a=-1.*(1.-exp(-lamssw*hs))/(CONDS*lamssw);
  b=-1.*exp(-lamssw*hs)*(1-exp(-lamisw*hi))/(CONDI*lamisw);
  c=-1.*(1.-exp(-lamslw*hs))/(CONDS*lamslw);
  d=-1.*exp(-lamslw*hs)*(1-exp(-lamilw*hi))/(CONDI*lamilw);

/* --------------------------------------------------------------------
 * Solar radiation at bottom of snow pack. RHS of Patterson and Hamblin, eq. 7
 * -------------------------------------------------------------------- */

  *SWnet=sw*a1*(a+b)+sw*a2*(c+d);

/* --------------------------------------------------------------------
 * Solar radiation at bottom of snow/ice layer.
 * Patterson and Hamblin eq. 8 (qf-qo)
 * -------------------------------------------------------------------- */

  *SW_under_ice = ( a1*sw*(1-exp(-(lamssw*hs+lamisw*hi)))
		    +a2*sw*(1-exp(-(lamslw*hs+lamilw*hi))) );
}
  
int lakeice (double *tempi, double Tcutoff, double sw_ice,
	     double twater, double fracice, int dt, double snowflux,
	     double qw, double *energy_ice_melt_bot, double sdepth,
	     double SWabsorbed, int rec, dmy_struct dmy, double *qf,
	     double *ice_water_eq, double volume, double sarea)
/**********************************************************************
 * Calculate the growth and decrease in the lake ice cover.
 * Changed from original model since ice_melt() (based on VIC/DHSVM snow_melt())
 * now handles melt of snow and ice from the surface, for consistency with
 * the rest of VIC.  This routine now handles melt/freeze at the
 * bottom of the ice pack.
 *
 * Parameters :
 *
 * tempi      Temperature of the ice (K).
 * Tcutoff    Temperature at which water freezes (K).
 * sw_ice     Net short wave radiation over ice (W/m2).
 * hice               Ice height (m).
 * twater       Temperature of the top water layer (K).
 * *fracice   Revised fraction of the lake covered by ice due to growth/decrease of exisitng cover (-).
 * dt         Time step (s).
 * *qbot      Incoming short wave radiation at bottom of the ice (W/m2).
 * *qw                Heat storage in the lake ice (J/m3).
 * ds           Amount of snowmelt (m)
 * qnetice;   Heat flux absorbed into the ice (W/m2).
 *
 * Modifications:
 * 2007-Nov-06 Changed from original model since ice_melt() (based on
 *             VIC/DHSVM snow_melt()) now handles melt of snow and ice
 *             from the surface, for consistency with the rest of VIC.
 *             This routine now handles melt/freeze at the bottom of
 *             the ice pack.						LCB via TJB
 **********************************************************************/
{
  
  double condqw, evapl; 
  double qmelts;           /* Energy available to melt ice. ? */
  double dibot;            /* change in ice surface at the bottom of the pack. */
  double delta_ice_frac;   /* Change in ice covered fraction. */
  double delta_ice_depth;
  double xfrac;
  double tprev;
  double RefrozenWater;
  double new_water_eq;
     
  /**********************************************************************
   * Calculate fluxes at the base of the ice.
   **********************************************************************/
  
  /* --------------------------------------------------------------------
   * Flux of heat in the ice at the ice/water interface (P & H, eq. 8)
   * -------------------------------------------------------------------- */

  *qf = snowflux + sw_ice - SWabsorbed;  

  /* --------------------------------------------------------------------
   * Amount of heat used to melt the ice (positive means freezing).
   * -------------------------------------------------------------------- */
  
  *energy_ice_melt_bot = *qf - qw; // from Hostetler 1991

  /* --------------------------------------------------------------------
   * Calculate the growth of the ice pack at the bottom (in meters).
   * -------------------------------------------------------------------- */

  dibot=(*energy_ice_melt_bot/(RHOICE*Lf))*dt*SECPHOUR;

  /* --------------------------------------------------------------------
   * Calculate the water equivalent of the ice volume (in cubic meters).
   * Freezing occurs over surface area of water, not ice, if ice area exceeds
   * water area.
   * -------------------------------------------------------------------- */

  new_water_eq = dibot * sarea * (fracice) * RHOICE/RHO_W;

  /**********************************************************************
   * Calculate the new height of the ice pack and the fractional ice
   * cover.
   **********************************************************************/
  if(dibot > 0.0) { /*Freezing; check if enough unfrozen water is available. */
    if(volume - *ice_water_eq >= new_water_eq ) {
      *ice_water_eq += new_water_eq;
    }
    else { /* Freezing is restricted by available water. */
      if( (volume - *ice_water_eq) > 0.0) {
        dibot = (volume - *ice_water_eq) / (sarea*(fracice)*RHOICE/RHO_W);
        *ice_water_eq = volume;
      }
      else {
        dibot = 0.0;
      }
    }
  }
  else { /* Melt */

    *ice_water_eq += new_water_eq;

    // check that ice not completely melted in the current time step
    if ( *ice_water_eq <= 0.0 ) {
      *ice_water_eq = 0.0;
    }

  }

  // *energy_ice_melt_bottom is not currently adjusted if there is not enough water to freeze or ice to melt.
  // This energy should go to ? and warming the water?, respectively.

  return (0);

}

float lkdrag (float Tsurf, double Tair, double wind, double roughness, double Z1)
{

  /**********************************************************************
   * Calculate the lake drag coefficient.
   *
   * Parameter :
   *
   * Tsurf   	Lake surface temperature (K).
   * Tair		Air temperature (K).
   * wind		Wind speed (m/s).
   * dragcoeff	Drag coefficient (-).
   **********************************************************************/

     double cdrn, ribn, ribd, rib;
     double cdr, cdrmin;
      
/**********************************************************************
 * Calculate the Richardson number.
 **********************************************************************/

     cdrn=(von_K/log(Z1/roughness))*(von_K/log(Z1/roughness)); /* dimensionless */

     ribn=Z1*G*(1.-Tsurf/Tair);           /* m2/s2 */

      if ((Tsurf/Tair) <= 1.0) {
	ribd=wind*wind+0.1*0.1;
      }
      else {
	ribd=wind*wind+1.0*1.0;
	  }

      rib=ribn/ribd;                /* dimensionless */

/**********************************************************************
 * Calculate the drag coefficient using the Richardson number.
 **********************************************************************/

      if (rib < 0.) {

         cdr=cdrn*(1.0+24.5*sqrt(-cdrn*rib));
      }

      else {
         cdr=cdrn/(1.0+11.5*rib);
	   }

      if((.25*cdrn) > 6.e-4)
	cdrmin=.25*cdrn;
      else
	cdrmin=6.e-4;

      if (cdr < cdrmin) 
	cdr=cdrmin;

      return(cdr);

   }


void rhoinit(double *tfsp, 
	     double pressure)
{

/**********************************************************************
 * Calculate the temperature at which water 
 * freezes depending on salinity and air pressure.
 *
 * Paramters :
 *
 * tfsp		Lake freezing point (C).
 * pressure     Air pressure (Pa)
 **********************************************************************/

  double salinity;

/**********************************************************************
 * Salinity is assumed to be zero.
 **********************************************************************/

  salinity = 0.;

/**********************************************************************
 * Calculate the lake freezing temperature (C). Pressure in bars.
 **********************************************************************/

  *tfsp = ( -0.0575 * salinity + 1.710523e-3 * pow(salinity,1.5) 
	    - 2.154996e-4 * salinity * salinity - 7.53e-3 * (pressure) / 100. );
      

}
 
double specheat (double t)
    {

/**********************************************************************
 * Calculate the specific heat of the water depending on water
 * temperature.  Salinity is assumed to be zero.
 *
 * Paramterers :
 *
 * t	Water temperature (C).
 **********************************************************************/

      double cpt;               /* Specific heat (J/Kg K) */

      cpt=4217.4 - 3.720283*t + 0.1412855*t*t - 2.654387e-3*t*t*t +
	2.093236e-5*t*t*t*t;

      return cpt;
    }

void temp_area(double sw_visible, double sw_nir, double surface_force, 
	       double *T, double *Tnew, double *water_density, double *de,
	       int dt, double *surface, int numnod, double dz, double surfdz, double *temph, 
	       double *cp, double *energy_out_bottom)
    {
/********************************************************************** 				       
  Calculate the water temperature for different levels in the lake.
 
  Parameters :
 
  sw_visible   Shortwave rad in visible band entering top of water column
  sw_nir       Shortwave rad in near infrared band entering top of water column
  surface_force The remaining rerms i nthe top layer energy balance 
  T		Lake water temperature at different levels (K).
  water_density		Water density at different levels (kg/m3).
  de		Diffusivity of water (or ice) (m2/d).
  dt		Time step size (s).
  surface	Area of the lake at different levels (m2).
  numnod	Number of nodes in the lake (-).
  dz        Thickness of the lake layers. 

  Modifications:
  2007-Apr-23 Added initialization of temph.				TJB
  2007-Oct-24 Modified by moving closing bracket for if ( numnod==1 ) up
	      so that the code actually calls energycalc() even if the
	      lake is represented by only one node.			KAC via TJB

 **********************************************************************/

      double z[MAX_LAKE_NODES], zhalf[MAX_LAKE_NODES];
      double a[MAX_LAKE_NODES], b[MAX_LAKE_NODES], c[MAX_LAKE_NODES];
      double d[MAX_LAKE_NODES];
      double told[MAX_LAKE_NODES];
     
      double dist12;
      int k;
      double surface_1, surface_2, surface_avg, T1;
      double cnextra;
      double swtop; /* The solar radiation at the top of the water column. */
      double top, bot; /* The depth of the top and the bottom of the current water layer. */
      double water_density_new;
      double jouleold, joulenew;
      double energyinput;
      double energymixed;
      double term1, term2;
      double totalenergy;

/**********************************************************************
 * Calculate the distance between the centers of the surface and first
 * lake layers.
 **********************************************************************/

/**********************************************************************
 * Initialize the water density at all the nodes and the depth of all
 * and distance between all nodes.
 **********************************************************************/

      for(k=0; k<numnod; k++) {
	 if(k==0)
	   z[k] = surfdz;
	 else
	   z[k]=dz;
         zhalf[k]=dz;
      }
      zhalf[0]=0.5*(z[0]+z[1]);
      energyinput=0.0;
      energymixed = 0.0;

/**********************************************************************
 * Calculate the right hand side vector in the tridiagonal matrix system
 * of equations.
 **********************************************************************/

	
      surface_1 = surface[0];
      surface_2 = surface[1];
      surface_avg = (surface_1 + surface_2)/2.;

       T1 = (sw_visible*(1*surface_1-surface_2*exp(-lamwsw*surfdz)) + 
      	    sw_nir*(1*surface_1-surface_2*exp(-lamwlw*surfdz)))/surface_avg 
      	+ (surface_force*surface_1)/surface_avg;          /* W/m2 */

       totalenergy = sw_visible  + sw_nir + surface_force;
       energyinput +=T1*surface_avg;
       cnextra = 0.5*(surface_2/surface_avg)*(de[0]/zhalf[0])*((T[1]-T[0])/z[0]);
      
       energymixed += cnextra;

       *temph=0.0;
       
       if(numnod==1)
	 Tnew[0] = T[0]+(T1*dt*SECPHOUR)/((1.e3+water_density[0])*cp[0]*z[0]);
       else {	
	 
	 /* --------------------------------------------------------------------
	  * First calculate d for the surface layer of the lake.
	  * -------------------------------------------------------------------- */
	 
	 d[0]= T[0]+(T1*dt*SECPHOUR)/((1.e3+water_density[0])*cp[0]*z[0])+cnextra*dt*SECPHOUR;
	 
	 *energy_out_bottom = (surface_1 - surface_2)*(sw_visible*exp(-lamwsw*surfdz) + 
						       sw_nir*exp(-lamwlw*surfdz));
	 
	 
	 
/* --------------------------------------------------------------------
 * Calculate d for the remainder of the column.
 * --------------------------------------------------------------------*/

/* ....................................................................
 * All nodes but the deepest node.
 * ....................................................................*/

      for(k=1; k<numnod-1; k++) {
	top = (surfdz+(k-1)*dz);
        bot = (surfdz+(k)*dz);

         surface_1 = surface[k]; 
	 surface_2 = surface[k+1];
         surface_avg =( surface[k]  + surface[k+1]) / 2.;


	  T1 = (sw_visible*(surface_1*exp(-lamwsw*top)-surface_2*exp(-lamwsw*bot)) + 
	       sw_nir*(surface_1*exp(-lamwlw*top)-surface_2*exp(-lamwlw*bot)))/surface_avg; 

	  energyinput +=T1*surface_avg;
	  term1 = 0.5 *(1./surface_avg)*((de[k]/zhalf[k])*((T[k+1]-T[k])/z[k]))*surface_2;
	  term2 = 0.5 *(-1./surface_avg)*((de[k-1]/zhalf[k-1])*((T[k]-T[k-1])/z[k]))*surface_1;
	 
	  cnextra = term1 + term2;
	
	 energymixed += (term1+term2);
         d[k]= T[k]+(T1*dt*SECPHOUR)/((1.e3+water_density[k])*cp[k]*z[k])+cnextra*dt*SECPHOUR;

	 	 *energy_out_bottom += (surface_1 - surface_2)*(sw_visible*exp(-lamwsw*bot) + 
	 						sw_nir*exp(-lamwlw*bot));

      }
/* ....................................................................
 * Calculation for the deepest node.
 * ....................................................................*/
       k=numnod-1;
      surface_1 = surface[k];
      surface_2 = surface[k];
      surface_avg = surface[k];
     
      top = (surfdz+(k-1)*dz);
      bot = (surfdz+(k)*dz);

      //       T1 = (sw_visible*(exp(-lamwsw*top)) + 
      //	    sw_nir*(exp(-lamwlw*top)))*surface_1/surface_avg; 

      T1 = (sw_visible*(surface_1*exp(-lamwsw*top)-surface_2*exp(-lamwsw*bot)) + 
	    sw_nir*(surface_1*exp(-lamwlw*top)-surface_2*exp(-lamwlw*bot)))/surface_avg; 

      energyinput +=T1*surface_avg;
      energyinput /= surface[0];

      cnextra = 0.5 * (-1.*surface_1/surface_avg)*((de[k-1]/zhalf[k-1])*((T[k]-T[k-1])/z[k]));

            *energy_out_bottom = 0.;
      *energy_out_bottom += surface_2*(sw_visible*exp(-lamwsw*bot) + sw_nir*exp(-lamwlw*bot));
      *energy_out_bottom /= surface[0];

      energymixed += cnextra;
    
      d[k] = T[k]+(T1*dt*SECPHOUR)/((1.e3+water_density[k])*cp[k]*z[k])+cnextra*dt*SECPHOUR;

/**********************************************************************
 * Calculate arrays for tridiagonal matrix.
 **********************************************************************/

/* --------------------------------------------------------------------
 * Top node of the column.
 * --------------------------------------------------------------------*/

     
      surface_2 = surface[1];
      surface_avg = (surface[0] + surface[1] ) / 2.;

      b[0] = -0.5 * ( de[0] / zhalf[0] ) *
	(  dt*SECPHOUR / z[0] ) * surface_2/surface_avg;
      a[0] = 1. - b[0];

/* --------------------------------------------------------------------
 * Second to second last node of the column.
 * --------------------------------------------------------------------*/

      for(k=1;k<numnod-1;k++) {
         surface_1 = surface[k];
         surface_2 = surface[k+1];
	 surface_avg = ( surface[k]  + surface[k+1]) / 2.;

         b[k] = -0.5 * ( de[k] / zhalf[k] ) *
	   (  dt*SECPHOUR / z[k] )*surface_2/surface_avg;
         c[k] = -0.5 * ( de[k-1] / zhalf[k-1] ) *
	   (  dt*SECPHOUR / z[k] )*surface_1/surface_avg;
         a[k] = 1. - b[k] - c[k];

	   }
/* --------------------------------------------------------------------
 * Deepest node of the column.
 * --------------------------------------------------------------------*/

      surface_1 = surface[numnod-1];
      surface_avg = surface[numnod-1];
      c[numnod-1] = -0.5 * ( de[numnod-1] / zhalf[numnod-1] ) *
	(  dt*SECPHOUR / z[numnod-1] ) * surface_1/surface_avg;
      a[numnod-1] = 1. - c[numnod-1];

/**********************************************************************
 * Solve the tridiagonal matrix.
 **********************************************************************/
     
      tridia(numnod,c,a,b,d,Tnew);

    }

/**********************************************************************
 * Adjust energy fluxes for change in density -> should be fixed by
 * moving to lagrangian scheme
 **********************************************************************/
    
    energycalc(Tnew, &joulenew, numnod,dz, surfdz, surface, cp, water_density);
     
    *temph=0.0;

    *temph = joulenew;

}

void tracer_mixer (double *T, 
                   int    *mixdepth,
		   int     freezeflag,
                   double *surface,
                   int     numnod, 
                   double  dz, 
		   double surfdz, 
                   double *cp) {
/**********************************************************************
 * Simulate the convective mixing in the lake.
 *
 * Paramters :
 *
 * T		Water temperatures (K).
 * water_density		Water densities (kg/m3).
 * mixdepth	         Top depth of local instability (node number, -).
 * freezeflag	         0 for ice, 1 for liquid water.
 * surface	Area of the lake per node number (m2).
 * numnod	         Number of nodes in the lake (-).
 **********************************************************************/

  int    k,j,m;             /* Counter variables. */
  int    mixprev;
  double avet, avev;
  double vol; /* Volume of the surface layer (formerly vol_tr). */
  double heatcon; /*( Heat content of the surface layer (formerly vol). */ 
  double Tav, densnew;
  double rho_max;
  double water_density[MAX_LAKE_NODES];
  
  for ( k = 0; k < numnod; k++ )
    water_density[k] = calc_density(T[k]);

/**********************************************************************
 * Initialize the top depth of local instability.
 **********************************************************************/

  mixprev = 0;
      
  for ( k = 0; k < numnod-1; k++ ) {

/**********************************************************************
 * Check for instability at each slice in water column.
 **********************************************************************/
    avet=0.0;
    avev=0.0;
	
    if ( water_density[k] > water_density[k+1] ) {

/* --------------------------------------------------------------------
 * If there is instability apply the mixing scheme.
 * --------------------------------------------------------------------*/

      if (mixprev == 0 && (k+1) > *mixdepth) {

/* ....................................................................
 * Correct the top depth of local instability.
 * ....................................................................*/
	*mixdepth = k+1;

      }

/*----------------------------------------------------------------------
 * Mix from mixprev to k+1  
 *----------------------------------------------------------------------*/
	   
      for ( m = mixprev; m <= k+1; m++ ) {

/* --------------------------------------------------------------------
 * Apply the mixing scheme from the previous depth to the instability
 * up to this node.
 * --------------------------------------------------------------------*/

	if (m == 0) {
	  /* Calculate the heat content and volume of the surface layer. */
	  heatcon = surfdz*(1.e3+water_density[m])*cp[m]*surface[m];
	  vol = surfdz * surface[m];
	}
	else {
	  /* Calculate the heat content and volume of all layers but 
	     the surface layer. */
	  heatcon = dz*(1.e3+water_density[m])*cp[m]*surface[m];
	  vol = dz * surface[m];
	}
	
	/* Calculate the volumetric weighed average lake temperature. */
	avet = avet+T[m]*heatcon;
	avev = avev+heatcon;
	
      }
      
      Tav = avet / avev;

/* --------------------------------------------------------------------
 * Calculate the density of the surface layer.
 * --------------------------------------------------------------------*/

      densnew = calc_density(Tav);

/* --------------------------------------------------------------------
 * Calculate the maximum density up to the local level.
 * --------------------------------------------------------------------*/

      rho_max = 0.0;
	
      for ( j = 0; j < mixprev; j++ ) {
	if ( (1000.+water_density[j]) > rho_max ) {
	  rho_max=1000.+water_density[j];
	}
      }

/* --------------------------------------------------------------------
 * Adjust temperatures and density in the mixed part of column.
 * --------------------------------------------------------------------*/

      for ( j = mixprev; j <= k+1; j++ ) {
	T[j]=Tav;
	water_density[j] = densnew;
      }

/* --------------------------------------------------------------------
 * Check to make sure that the mixing has not generated new instabilities
 * above the previous depth to the instabilities.
 * --------------------------------------------------------------------*/

      if (rho_max > (1000.+densnew)) {
	
	/* If there are still instabilities iterate again..*/
	mixprev = 0;
	k=-1;
      }
    }

    else {
/**********************************************************************
 * If there are no instabilities up to now then the depth to the
 * instability has to be increased by 1 node.
 **********************************************************************/

      mixprev=k+1;
    }
  }

/**********************************************************************
 * Recalculate the water density.
 **********************************************************************/

  for ( k = 0; k < numnod; k++ ) {
    water_density[k] = calc_density(T[k]);
  }
}      

void tridia (int ne, 
	     double *a, 
	     double *b, 
	     double *c, 
	     double *y, 
	     double *x) {
/**********************************************************************
 * Solve a tridiagonal system of equations.
 * Parameters :
 * ns		The number of systems to be solved.
 * nd		First dimension of arrays (larger than or equal to ns).
 * ne		The number of unknowns in each system. This must 
 *		be larger than or equal to 2.
 * a		The sub diagonals of the matrices are stored in locations
 *		a(j,2) through a(j,ne).
 * b		The main diagonals of the matrices are stored in
 *		locations b(j,1) through b(j,ne).
 * c		The super diagonals of the matrices are stored in
 *		locations c(j,1) through c(j,ne-1).
 * y		The right hand side of the equations is stored in
 *		y(j,1) through y(j,ne).
 * x		The solutions of the systems are returned in
 *		locations x(j,1) through x(j,ne).
 *
 * History : Based on a streamlined version of the old NCAR ULIB subroutine
 * TRDI used in the PHOENIX climate model of Schneider and Thompson
 * (J.G.R., 1981). Revised by Starley Thompson to solve multiple systems
 * and vectorize well on the CRAY-1. Later revised to include a PARAMETER
 * statement to define loop limits and thus enable Cray short vector
 * loops.
 *
 * Algorithm:  LU decomposition followed by solution.  NOTE: This
 * subroutine executes satisfactorily if the input matrix is diagonally
 * dominant and non-singular.  The diagonal elements are used to pivot, and
 * no tests are made to determine singularity. If a singular or numerically
 * singular matrix is used as input a divide by zero or floating point
 * overflow will result.
 **********************************************************************/

     double alpha[MAX_LAKE_NODES], gamma[MAX_LAKE_NODES]; /* Work arrays dimensioned (nd,ne).*/

     int nm1, i, ib;
  
     nm1 = ne-1;

/**********************************************************************
 * Obtain the LU decompositions.
 **********************************************************************/

	  alpha[0] = 1./b[0];
	  gamma[0] = c[0]*alpha[0];


      for(i=1; i<nm1;i++) {
	  alpha[i] = 1./(b[i]-a[i]*gamma[i-1]);
	  gamma[i] = c[i]*alpha[i];
      }

/**********************************************************************
 * Solve the system.
 **********************************************************************/

      x[0] = y[0]*alpha[0];
	
      for(i=1; i<nm1; i++) {
	  x[i] = (y[i]-a[i]*x[i-1])*alpha[i];
      }

    
      x[nm1] = (y[nm1]-a[nm1]*x[nm1-1])/
	   (b[nm1]-a[nm1]*gamma[nm1-1]);

      for(i=nm1-1; i>=0; i--) {
            x[i] = x[i]-gamma[i]*x[i+1];
      }

}  


void energycalc (double *finaltemp, double *sumjoule, int numnod, double dz, double surfdz, double *surface, double *cp, double *density)
{
  
  /**********************************************************************
   * Calculate the thermal energy in the column.
   * finaltemp[numnodes]   Average temp for each node of the water column.
   * sumjoule              Thermal energy of water column, in joules 
   * numnod   
   **********************************************************************/
  
  double energy;
  int k;
  
  *sumjoule=0.0;
  
  for(k=0; k<numnod;k++)
    {
      // density = calc_density(finaltemp[k]);
      
      if(k==0)
	energy=(finaltemp[k]+KELVIN)*surfdz*(1000.+density[k])*cp[k]*(surface[k]+surface[k+1])/2.;
      else if (k==numnod-1)
	energy=(finaltemp[k]+KELVIN)*dz*(1000.+density[k])*cp[k]*surface[k]/2.;
      else
	energy=(finaltemp[k]+KELVIN)*dz*(1000.+density[k])*cp[k]*(surface[k]+surface[k+1])/2.;
      
      *sumjoule+=energy;
    }
}

int water_balance (lake_var_struct *lake, lake_con_struct lake_con, int dt, dist_prcp_struct *prcp,
		    int rec, int iveg,int band, double lakefrac, soil_con_struct soil_con, 
#if EXCESS_ICE
		    int SubsidenceUpdate, double total_meltwater,
#endif
		    double prec,  double oldvolume, double deltasnow, double vapor_flux, double fracprv)
/**********************************************************************
 * This routine calculates the water balance of the lake
 
  Modifications:
  2007-Aug-09 Added features for EXCESS_ICE option.				JCA
  2007-Oct-24 Added rec to call for water_balance so that error and warning
	      messages can report the model record number.			KAC via TJB
  2007-Oct-24 Modified loop for get_sarea to include the active node.		KAC via TJB
  2007-Nov-06 Ice water content now tracked via ice_water_eq.  Lake area is
	      now the larger of lake->areai and lake->sarea.  Drainage is now
	      modeled as flow over a broad-crested wier.			LCB via TJB
  2008-Apr-21 Corrected some omissions from the previous mods.			LCB via TJB
**********************************************************************/
{
  extern option_struct   options;
  int isave_n;
  double d_area, d_level, d_volume;
  double runoff_volume;
  double tempvolume,surfacearea,ldepth;
  double remain;
  double i;
  double m;
  float index;
  double newdepth;
  int j,k;
  double in;
  double Tnew[MAX_LAKE_NODES];
  double tempdepth;
  double wetland_runoff;
  double inflow, outflow, storage;
  double bpercent;
  double inflec, midrate;
  double circum;
  double snowmltvolume, evapvolume;
  double newfraction, Recharge;
  int ErrorFlag;
  cell_data_struct    ***cell;
  int lindex;
  double frac;

  cell    = prcp->cell;
  
  /**********************************************************************
   * 1. convert runoff input/output volume to rate
   **********************************************************************/
  
  isave_n = lake->activenod;   /* save initial no. of nodes for later */
  
  /* convert from mm over grid cell tomm over lake. */
  
  wetland_runoff = 0.0;
  
  if(fabs(lakefrac - 1.0) > SMALL) {
    wetland_runoff += cell[WET][iveg][band].runoff + cell[WET][iveg][band].baseflow;
    wetland_runoff *= lake_con.Cl[0]*(1.-lakefrac);
  }
  
  /* Assume baseflow enters lake in same proportion as runoff. */ 
  bpercent = lake_con.rpercent;
  
  /* Convert from mm to m3. */
#if EXCESS_ICE
  if(SubsidenceUpdate > 0 ) 
    runoff_volume = ( ( ( lake->runoff_in * lake_con.rpercent + lake->baseflow_in*bpercent + wetland_runoff + total_meltwater*lake_con.Cl[0]*lakefrac) / 1000. ) * soil_con.cell_area );  
  else
#endif
    runoff_volume = ( ( ( lake->runoff_in * lake_con.rpercent + lake->baseflow_in*bpercent + wetland_runoff) / 1000. ) * soil_con.cell_area );  
  
  
  /**********************************************************************
   * 2. calculate change in lake level
   *     grid cell runoff (m3/TS for non-lake area)
   *     snow meltwater (mm)
   *     open water evaporation (mm)
   *     precip (added in solvelake)  
   **********************************************************************/
  
  oldvolume = lake->volume;
  evapvolume = (lake->evapw/1000.) * lake->sarea;
  snowmltvolume = (lake->snowmlt/1000.) * lake->areai;

  // Add runoff from rest of grid cell and wetland to lake, remove evaporation
  if( evapvolume > (lake->volume + runoff_volume-lake->ice_water_eq) + snowmltvolume) {
    lake->evapw = 1000.* ((lake->volume + runoff_volume-lake->ice_water_eq) + snowmltvolume)/lake->sarea;
     evapvolume = (lake->evapw/1000.) * lake->sarea;
    lake->volume = lake->ice_water_eq;
  }
  else {
    lake->volume += ( runoff_volume + snowmltvolume-evapvolume);
  }
  inflow = ( runoff_volume + snowmltvolume-evapvolume);

  // Find new lake depth and surface area of the unfrozen portion for recharge calculations
  if(1.- fracprv < FRACLIM)
    ErrorFlag = get_depth(lake_con, lake->volume-lake->ice_water_eq, &ldepth);
  else
    ErrorFlag = get_depth(lake_con, lake->volume, &ldepth);
  if ( ErrorFlag == ERROR ) {
    fprintf(stderr, "Something went wrong in get_depth; record = %d, volume = %f, depth = %e\n",rec,lake->volume,ldepth);
    return ( ErrorFlag );
  }
  ErrorFlag = get_sarea(lake_con, ldepth, &surfacearea);
  if ( ErrorFlag == ERROR ) {
    fprintf(stderr, "Something went wrong in get_sarea; record = %d, depth = %f, sarea = %e\n",rec,ldepth,surfacearea);
    return ( ErrorFlag );
  }

  // Extract moisture for recharge first
  if(lake->areai > surfacearea)
    newfraction = lake->areai/lake_con.basin[0];
  else
    newfraction = surfacearea/lake_con.basin[0];

  if(newfraction > lakefrac) {

  // Check if sufficient liquid volume exists to 'fill' newly flooded soil
  // If not enough liquid water level goes to zero.

    Recharge = 0.0;
    for(j=0; j<options.Nlayer; j++) {
      Recharge += ((soil_con.max_moist[j]-cell[WET][iveg][band].layer[j].moist ) / 1000. )
        * (newfraction-lakefrac) * lake_con.basin[0];
    }

    if(lake->volume-lake->ice_water_eq > Recharge) {
      lake->volume -= Recharge;

      /* Update soil moisture storage terms. */
      for(j=0; j<options.Nlayer; j++) {
        cell[DRY][iveg][band].layer[j].moist += (soil_con.max_moist[j] -cell[DRY][iveg][band].layer[j].moist )*(newfraction-lakefrac)/(1.-lakefrac);
        cell[WET][iveg][band].layer[j].moist += (soil_con.max_moist[j] -cell[WET][iveg][band].layer[j].moist )*(newfraction-lakefrac)/(1.-lakefrac);
      }
    }
    else {

      Recharge = (lake->volume-lake->ice_water_eq)* 1000./((1.-lakefrac)*lake_con.basin[0]);
      lake->volume = lake->ice_water_eq;

      for(j=0; j<options.Nlayer; j++) {

        if(Recharge > (soil_con.max_moist[j]-cell[WET][iveg][band].layer[j].moist)) {
          Recharge -= (soil_con.max_moist[j]-cell[WET][iveg][band].layer[j].moist);
          cell[WET][iveg][band].layer[j].moist = soil_con.max_moist[j];
        }
        else {
          cell[WET][iveg][band].layer[j].moist += Recharge;
          Recharge = 0.0;
        }
      }

    }
  }

  /**********************************************************************
   * 3. Calculate runoff from lake.  Runoff out is in meters over lake surface.
   * Based on the equation for flow over a broad crested weir.
   **********************************************************************/

  // If wetland baseflow is negative, it has already been extracted from the lake. */
  if(cell[WET][iveg][band].baseflow < 0.0)
    lake->baseflow_out = 0.0;
  else {
    lindex = options.Nlayer-1;
    frac = soil_con.Ds * soil_con.Dsmax
      / (soil_con.Ws * soil_con.max_moist[lindex]);
    lake->baseflow_out = frac * (soil_con.max_moist[lindex] - soil_con.resid_moist[lindex] );

    lake->baseflow_out += (soil_con.Dsmax - soil_con.Ds * soil_con.Dsmax / soil_con.Ws);
    lake->baseflow_out *= lake_con.Cl[0]*(lakefrac)* soil_con.cell_area/1000.;

    // extract baseflow volume from the lake m^3

    if(lake->volume -lake->ice_water_eq >= lake->baseflow_out)
      lake->volume -= ( lake->baseflow_out );
    else {
      lake->baseflow_out = lake->volume - lake->ice_water_eq;
      lake->volume -= ( lake->baseflow_out );
    }
  }

  // Find new lake depth and surface area of the unfrozen portion for runoff calculations
  ErrorFlag = get_depth(lake_con, lake->volume-lake->ice_water_eq, &ldepth);
  if ( ErrorFlag == ERROR ) {
    fprintf(stderr, "Something went wrong in get_depth; record = %d, volume = %f, depth = %e\n",rec,lake->volume,ldepth);
    return ( ErrorFlag );
  }

  // Compute runoff volume in m^3 and extract runoff volume from lake
  if(ldepth <= lake_con.mindepth)
    lake->runoff_out = 0.0;
  else {
    circum=2*PI*pow(surfacearea/PI,0.5);
    lake->runoff_out = lake_con.wfrac*circum*SECPHOUR*((double)dt)*1.6*pow(ldepth-lake_con.mindepth, 1.5);
    if(lake->volume - lake->ice_water_eq-lake_con.minvolume >= lake->runoff_out)
      lake->volume -= ( lake->runoff_out );
    else {
      lake->runoff_out = lake->volume - lake->ice_water_eq-lake_con.minvolume ;
      lake->volume -= ( lake->runoff_out );
    }
  }
  outflow = lake->runoff_out + lake->baseflow_out;

  // check that lake volume does not exceed its maximum
  if (lake->volume - lake_con.maxvolume > SMALL) {
    if(lake->ice_water_eq > lake_con.maxvolume) {
      lake->runoff_out += (lake->volume - lake->ice_water_eq);
      lake->volume = lake->ice_water_eq;
    }
    else {
      lake->runoff_out += (lake->volume - lake_con.maxvolume);
      lake->volume = lake_con.maxvolume;
    }
  }
  else if (lake->volume < SMALL)
    lake->volume = 0.0;

  /**********************************************************************/
  /* End of runoff calculation */
  /**********************************************************************/

  // adjust from m^3 to mm over entire gridcell
  lake->runoff_out = ( lake->runoff_out ) * 1000. / soil_con.cell_area;

  if(cell[WET][iveg][band].baseflow < 0.0)
     lake->baseflow_out = cell[WET][iveg][band].baseflow* lake_con.Cl[0]*(1.-lakefrac);
  else
    lake->baseflow_out *= 1000. / soil_con.cell_area;

  
  /* Add in runoff and baseflow that did not pass through the lake. */
  lake->runoff_out   += (1. - lake_con.rpercent) * lake->runoff_in;
  lake->baseflow_out  += (1. - bpercent) * lake->baseflow_in;
  
  // Recalculate lake depth
  if(1.-fracprv < FRACLIM)
    ErrorFlag = get_depth(lake_con, lake->volume-lake->ice_water_eq, &(lake->ldepth));
  else
    ErrorFlag = get_depth(lake_con, lake->volume, &(lake->ldepth));
  if ( ErrorFlag == ERROR ) {
    fprintf(stderr, "Something went wrong in get_depth; record = %d, volume = %f, depth = %e\n",rec,lake->volume,lake->ldepth);
    return ( ErrorFlag );
  }

  /**********************************************************************
   *  4. Adjust the activenodes and lake area array. 
   **********************************************************************/
 
  if(lake->ldepth > MAX_SURFACE_LAKE && lake->ldepth < 2*MAX_SURFACE_LAKE) {
    /* Not quite enough for two full layers. */
    lake->surfdz = lake->ldepth/2.;
    lake->dz = lake->ldepth/2.;
    lake->activenod = 2;
  }
  else if(lake->ldepth >= 2* MAX_SURFACE_LAKE) {
    /* More than two layers. */	
    lake->surfdz = MAX_SURFACE_LAKE;
    lake->activenod = (int) (lake->ldepth/MAX_SURFACE_LAKE);
    if(lake->activenod > MAX_LAKE_NODES)
      lake->activenod = MAX_LAKE_NODES;
    lake->dz = (lake->ldepth-lake->surfdz)/((float)(lake->activenod-1));
  }	
  else if(lake->ldepth > SMALL) {
    lake->surfdz = lake->ldepth;
    lake->dz = 0.0;
    lake->activenod = 1;
  }
  else {
    lake->runoff_out += ( lake->volume -lake->ice_water_eq) * 1000. / soil_con.cell_area;
    lake->volume = lake->ice_water_eq;
    lake->surfdz = 0.0;
    lake->dz = 0.0;
    lake->activenod = 0;
    lake->ldepth = 0.0;
  } 	

  // lake_con.basin equals the surface area at specific depths as input by
  // the user in the lake parameter file or calculated in read_lakeparam(), 
  // lake->surface equals the area at the top of each dynamic solution layer 
  
  /* Re-calculate lake surface area  */
  for(k=0; k<=lake->activenod; k++) {
    if(k==0)
      ldepth = lake->ldepth;
    else
      ldepth = lake->dz*(lake->activenod - k);
    
    ErrorFlag = get_sarea(lake_con, ldepth, &(lake->surface[k]));
    if ( ErrorFlag == ERROR ) {
      fprintf(stderr, "Something went wrong in get_sarea; record = %d, depth = %f, sarea = %e\n",rec,ldepth,lake->surface[k]);
      return ( ErrorFlag );
    }
  }
  
  
  /*******************************************************************/  
  /* Adjust temperature distribution if number of nodes has changed. 
     Note:  This approach (and the lake model in general) does not preserve 
     the thermal energy of the water column. */
 
  index = 0.;
  if ( lake->activenod != isave_n ) {
    for(k=0; k< lake->activenod; k++) {
      Tnew[k] = 0.0;
      for(i=0; i< isave_n; i++) {
	index += (1./lake->activenod);
	Tnew[k] += lake->temp[(int)floor(index)];
      }
    }
    for(k=0; k< lake->activenod; k++) {
      if(isave_n > 0)
	lake->temp[k] = Tnew[k]/isave_n;
      else
	lake->temp[k] = prcp->energy[iveg][band].Tsurf;
    }	
  }
 
  if(lake->activenod == isave_n && isave_n == 0)
    lake->temp[k] = prcp->energy[iveg][band].Tsurf;	
 
  // DEBUG water balance
  storage = oldvolume - lake->volume;
 
  return(0);

}


void update_prcp(dist_prcp_struct *prcp, 
		 energy_bal_struct *lake_energy, 
		 snow_data_struct *lake_snow, 
		 double lakefraction,	
		 lake_var_struct *lake,
		 lake_con_struct lake_con,
		 int iveg,
		 int band,
		 soil_con_struct soil_con) 
/**************************************************************************
  This function combines the energy balance and snow data structures
  such that the fluxes leaving the lake module represent a weighted average
  from the open water and wetland portions.

  Modifications:
  2007-Oct-24 Modified to include redistibution of snow density, without it
	      snow left on the shore as the lake retreats can cause nans to
	      appear as snow is updated without an initial density.		KAC via TJB
  2007-Oct-24 modified to set lake ice fraction to 0 if new lake fraction is
	      0, previously lake ice was always adjusted by dividing old by new
	      fractions, which causes nans when the new fraction is 0.		KAC via TJB
  2007-Oct-24 Modified redistribution of snow variables between wetland and
	      lake ice.  Surface water is now stored as a lake variable and
	      used to reestablish surf_water at the start of a new time step.	KAC via TJB
  2007-Nov-06 Lake ice is now explicitly counted in the lake mass balance, so
	      the lake cannot disappear if there is still ice cover, previously
	      the volume of frozen water was not tracked or extracted from the
	      liquid water portion of the lake, so ice could grow indefinitely.	LCB via TJB
  2008-Jan-23 Added updates of wland_snow->pack_water and
	      pack_temp in conjunction with 2-layer snow pack over
	      lake ice.								LCB via TJB
  2008-Apr-21 Added snow surf_temp, pack_temp, and coldcontent.			LCB via TJB
**************************************************************************/
{
  extern option_struct   options;
  energy_bal_struct    **wland_energy;
  snow_data_struct     **wland_snow;
  cell_data_struct    ***cell;
  veg_var_struct      ***veg_var;
  int i;
  int index;	
#if SPATIAL_FROST
  int frost_area;
#endif // SPATIAL_FROST
  // DEBUG variables
  double istor, finstor;
  
  wland_energy  = prcp->energy;
  wland_snow    = prcp->snow;
  cell    = prcp->cell;
  veg_var = prcp->veg_var;
  
  // DEBUG wetland water balance
  istor = wland_snow[iveg][band].swq*(1.-lakefraction) + lake_snow->swq*(lakefraction);
  istor += veg_var[WET][iveg][band].Wdew*(1.-lakefraction);
  istor += 1000.*lake->volume/lake_con.basin[0];
  for(i=0; i<options.Nlayer; i++) {
#if EXCESS_ICE
    istor += cell[WET][iveg][band].layer[i].moist*(1.-lakefraction) + (lakefraction)*soil_con.effective_porosity[i]*soil_con.depth[i]*1000.;
#else
    istor += cell[WET][iveg][band].layer[i].moist*(1.-lakefraction) + (lakefraction)*soil_con.porosity[i]*soil_con.depth[i]*1000.;
#endif // EXCESS_ICE
  }

  /* Update energy balance variables. */
  wland_energy[iveg][band].AlbedoLake = lake_energy->AlbedoLake;

  wland_energy[iveg][band].NetLongAtmos *= (1.-lakefraction);
  wland_energy[iveg][band].NetLongAtmos += lakefraction * lake_energy->NetLongAtmos;
  
  wland_energy[iveg][band].NetShortAtmos *= (1.-lakefraction);
  wland_energy[iveg][band].NetShortAtmos += lakefraction * lake_energy->NetShortAtmos;
  
  wland_energy[iveg][band].snow_flux *= (1.-lakefraction);
  wland_energy[iveg][band].snow_flux += lakefraction * lake_energy->snow_flux;
  
  wland_energy[iveg][band].deltaH *= (1.-lakefraction);
  wland_energy[iveg][band].deltaH += lakefraction * lake_energy->deltaH;
  
  wland_energy[iveg][band].deltaCC *= (1.-lakefraction);
  wland_energy[iveg][band].deltaCC += lakefraction * lake_energy->deltaCC;
  
  wland_energy[iveg][band].grnd_flux *= (1.-lakefraction);
  wland_energy[iveg][band].grnd_flux += lakefraction * lake_energy->grnd_flux;
  
  wland_energy[iveg][band].refreeze_energy *= (1.-lakefraction);
  wland_energy[iveg][band].refreeze_energy += lakefraction * lake_energy->refreeze_energy;
  
  wland_energy[iveg][band].advection *= (1.-lakefraction);
  wland_energy[iveg][band].advection += lakefraction * lake_energy->advection;

  wland_energy[iveg][band].AtmosLatent *= (1.-lakefraction);
  wland_energy[iveg][band].AtmosLatent += lakefraction * lake_energy->AtmosLatent;

  wland_energy[iveg][band].AtmosSensible *= (1.-lakefraction);
  wland_energy[iveg][band].AtmosSensible += lakefraction * lake_energy->AtmosSensible;

  wland_energy[iveg][band].Tsurf *= (1.-lakefraction);
  wland_energy[iveg][band].Tsurf += lakefraction * lake_energy->Tsurf;

  wland_energy[iveg][band].error *= (1.-lakefraction);
  wland_energy[iveg][band].error += lakefraction * lake_energy->error;

  /* Update snow variables. */
  wland_snow[iveg][band].albedo *= (1.-lakefraction);
  wland_snow[iveg][band].albedo += lakefraction * lake->SAlbedo;

  wland_snow[iveg][band].coverage *= (1.-lakefraction);
  wland_snow[iveg][band].coverage += lakefraction * lake_snow->coverage;
 
  wland_snow[iveg][band].vapor_flux *= (1.-lakefraction);
  wland_snow[iveg][band].vapor_flux += lakefraction * lake_snow->vapor_flux;

  wland_snow[iveg][band].canopy_vapor_flux *= (1.-lakefraction);

  wland_snow[iveg][band].blowing_flux *= (1.-lakefraction);	
  wland_snow[iveg][band].blowing_flux += lakefraction * lake_snow->blowing_flux;
  
  wland_snow[iveg][band].surface_flux *= (1.-lakefraction);
  wland_snow[iveg][band].surface_flux += lakefraction * lake_snow->surface_flux;
  
  wland_snow[iveg][band].swq *= (1.-lakefraction);
  wland_snow[iveg][band].swq += lakefraction * lake_snow->swq;
  lake->swe = lake_snow->swq;
  
  wland_snow[iveg][band].last_snow = lake_snow->last_snow;
  
  wland_snow[iveg][band].surf_temp *= (1.-lakefraction);
  wland_snow[iveg][band].surf_temp += lakefraction * lake_snow->surf_temp;
  lake->surf_temp = lake_snow->surf_temp;

  wland_snow[iveg][band].pack_temp *= (1.-lakefraction);
  wland_snow[iveg][band].pack_temp += lakefraction * lake_snow->pack_temp;
  lake->pack_temp = lake_snow->pack_temp;
  
  wland_snow[iveg][band].surf_water *= (1.-lakefraction);
  wland_snow[iveg][band].surf_water += lakefraction * lake_snow->surf_water;
  lake->surf_water = lake_snow->surf_water;
  
  wland_snow[iveg][band].pack_water *= (1.-lakefraction);
  wland_snow[iveg][band].pack_water += lakefraction * lake_snow->pack_water;
  lake->pack_water = lake_snow->pack_water;
  
  wland_snow[iveg][band].coldcontent *= (1.-lakefraction);
  wland_snow[iveg][band].coldcontent += lakefraction * lake_snow->coldcontent;
  lake->coldcontent = lake_snow->coldcontent;

  wland_snow[iveg][band].depth *= (1.-lakefraction);
  wland_snow[iveg][band].depth += lakefraction * lake_snow->depth; 
  lake->sdepth = lake_snow->depth;
  
  if ( wland_snow[iveg][band].density == 0 && wland_snow[iveg][band].depth > 0 )
     wland_snow[iveg][band].density = 1000. *  wland_snow[iveg][band].swq
       / wland_snow[iveg][band].depth;

  /* Update canopy storage terms. */
  veg_var[DRY][iveg][band].Wdew *= (1.-lakefraction);
  veg_var[WET][iveg][band].Wdew *= (1.-lakefraction);

  /* Update evap terms */
  veg_var[DRY][iveg][band].canopyevap *= (1.-lakefraction);
  veg_var[WET][iveg][band].canopyevap *= (1.-lakefraction);

  for(i=0; i<options.Nlayer; i++) {
    cell[DRY][iveg][band].layer[i].evap *= (1.-lakefraction);
    cell[WET][iveg][band].layer[i].evap *= (1.-lakefraction);
  }

  /* Recharge due to lake expansion is now addressed in the water_balance() */

  for(i=0; i<options.Nlayer; i++) {
    cell[DRY][iveg][band].layer[i].moist *= (1.-lakefraction);
    cell[WET][iveg][band].layer[i].moist *= (1.-lakefraction);
    cell[DRY][iveg][band].layer[i].moist += soil_con.max_moist[i]*lakefraction;
    cell[WET][iveg][band].layer[i].moist += soil_con.max_moist[i]*lakefraction;

    // Do not rescale ice content since it is not currently computed sub-lake.
    /*
      #if SPATIAL_FROST
      for ( frost_area = 0; frost_area < FROST_SUBAREAS; frost_area++ ) {
      cell[DRY][iveg][band].layer[i].ice[frost_area] *= (1.-lakefraction);
      cell[WET][iveg][band].layer[i].ice[frost_area] *= (1.-lakefraction);
      }
      #else
      cell[DRY][iveg][band].layer[i].ice *= (1.-lakefraction);
      cell[WET][iveg][band].layer[i].ice *= (1.-lakefraction);
      #endif */
  }

  // Debug water balance
  finstor = wland_snow[iveg][band].swq;
  finstor += veg_var[WET][iveg][band].Wdew;
  finstor += 1000.*lake->volume/lake_con.basin[0];
  for(i=0; i<options.Nlayer; i++) {
    finstor += cell[WET][iveg][band].layer[i].moist;
  }

//  // DEBUG water balance
//  if(fabs(finstor-istor) > SMALL) {
//    fprintf(stderr, "error = %f\n", finstor-istor);
//    exit(0);
//  }

} 


/**************************************************************************
 This function assigns values to the local lake variables from the
 structure of weighted averages from the previous time step.

  Modifications:
  2006-Nov-07 Assigned value to MELTING.					TJB
  2007-Aug-16 Added ErrorFlag to return value of
               distribute_node_moisture_properties.				JCA
  2007-Aug-21 Added features for EXCESS_ICE option.
              Including moving distribute_node_moisture_properties call to
              before surface_fluxes call in wetland_energy.			JCA
  2007-Oct-24 Modified redistribution of snow variables between wetland and
	      lake ice.  Surface water is now stored as a lake variable and
	      used to reestablish surf_water at the start of a new time step.	KAC via TJB
  2007-Nov-06 Redistribution now takes into account evaporative fluxes from
	      wetland.  Recharge due to lake expansion is now addressed in
	      water_balance().  Replaced lake.fraci with lake.areai.		LCB via TJB
  2008-Jan-23 Added initialization of lake_snow->pack_water and
	      pack_temp in conjunction with 2-layer snow pack over
	      lake ice.								LCB via TJB
  2008-Apr-21 Added snow surf_temp, pack_temp, and coldcontent.			LCB via TJB
**************************************************************************/
int initialize_prcp(dist_prcp_struct *prcp, 
		    energy_bal_struct *lake_energy, 
		    snow_data_struct *lake_snow, 
		    int iveg,
		    int band,
		    soil_con_struct soil_con, 
		    lake_var_struct *lake,
		    int rec,
		    lake_con_struct lake_con,
		    double *lakefrac,
		    double *fraci)
{
  extern option_struct   options;
  energy_bal_struct    **wland_energy;
  snow_data_struct     **wland_snow;
  cell_data_struct    ***cell;
  veg_var_struct      ***veg_var;
  int i, index;
#if SPATIAL_FROST
  int    frost_area;
#endif // SPATIAL_FROST
  double tmp_volume;
  double totalmoist;
  double newmoist;
  double moist[MAX_LAYERS];
  double error;
  double sum_ice;
  double istor, finstor;
  int ErrorFlag;

  wland_energy  = prcp->energy;
  wland_snow    = prcp->snow;
  cell    = prcp->cell;
  veg_var = prcp->veg_var;

  istor = wland_snow[iveg][band].swq + veg_var[WET][iveg][band].Wdew;
  for(i=0; i<options.Nlayer; i++) {
        istor += cell[WET][iveg][band].layer[i].moist;
  }

  /* Update energy balance variables. This is probably unnecessary. */
  lake_energy->NetLongAtmos = wland_energy[iveg][band].NetLongAtmos;
  
  lake_energy->NetShortAtmos =  wland_energy[iveg][band].NetShortAtmos;
  
  lake_energy->snow_flux =  wland_energy[iveg][band].snow_flux;

  lake_energy->deltaH =  wland_energy[iveg][band].deltaH;

  lake_energy->deltaCC =  wland_energy[iveg][band].deltaCC;

  lake_energy->grnd_flux =  wland_energy[iveg][band].grnd_flux;

  lake_energy->refreeze_energy =  wland_energy[iveg][band].refreeze_energy;

  lake_energy->advection =  wland_energy[iveg][band].advection;

  lake_energy->AtmosLatent =  wland_energy[iveg][band].AtmosLatent;

  lake_energy->AtmosSensible =  wland_energy[iveg][band].AtmosSensible;

  lake_energy->error = wland_energy[iveg][band].error;

  /* Update snow variables. */
  lake_snow->last_snow = wland_snow[iveg][band].last_snow;
  lake_snow->coverage = wland_snow[iveg][band].coverage;

  if(lake->areai > 0.0) {   
    lake_snow->vapor_flux = wland_snow[iveg][band].vapor_flux;
    lake_snow->blowing_flux = wland_snow[iveg][band].blowing_flux;
    lake_snow->surface_flux = wland_snow[iveg][band].surface_flux;
    lake_snow->surf_temp = wland_snow[iveg][band].surf_temp;
    lake_snow->pack_temp = wland_snow[iveg][band].pack_temp;
  }
  else {
    lake_snow->vapor_flux = 0.;  
    lake_snow->blowing_flux = 0.;  
    lake_snow->surface_flux = 0.;  
    lake_snow->surf_temp = 0.;
    lake_snow->pack_temp = 0.;
  }

  if((lake->areai >= lake->sarea) && lake->areai > 0.0) {
    *lakefrac = lake->areai/lake_con.basin[0];
    *fraci = 1.0;
  }
  else if (lake->areai > 0.0) {
    *fraci = lake->areai/lake->sarea;
    *lakefrac = lake->sarea/lake_con.basin[0];
  }
  else { /* ice area = 0.0 */
    *fraci = 0.0;
    *lakefrac = lake->sarea/lake_con.basin[0];
  }

  //if(*fraci > 0.99) *fraci = 1.0;
  //if(*lakefrac  > .99) *lakefrac = 1.0;

  /* Snow must be distributed between lake ice and wetland */
	
  if(*lakefrac < 1.0) {
    wland_snow[iveg][band].swq = (wland_snow[iveg][band].swq - lake->swe*(*lakefrac))/(1.-(*lakefrac));
    wland_snow[iveg][band].albedo = (wland_snow[iveg][band].albedo - lake->SAlbedo*(*lakefrac))/(1.-(*lakefrac));
    wland_snow[iveg][band].depth = (wland_snow[iveg][band].depth - lake->sdepth*(*lakefrac))/(1.-(*lakefrac));
    wland_snow[iveg][band].surf_water = (wland_snow[iveg][band].surf_water - lake->surf_water*(*lakefrac))/(1.-(*lakefrac));
    wland_snow[iveg][band].pack_water = (wland_snow[iveg][band].pack_water - lake->pack_water*(*lakefrac))/(1.-(*lakefrac));
    wland_snow[iveg][band].coldcontent = (wland_snow[iveg][band].coldcontent - lake->coldcontent*(*lakefrac))/(1.-(*lakefrac));
    wland_snow[iveg][band].surf_temp = (wland_snow[iveg][band].surf_temp - lake->surf_temp*(*lakefrac))/(1.-(*lakefrac));
    wland_snow[iveg][band].pack_temp = (wland_snow[iveg][band].pack_temp - lake->pack_temp*(*lakefrac))/(1.-(*lakefrac));

    if ( wland_snow[iveg][band].swq < SMALL ) {
      wland_snow[iveg][band].swq = 0.0;
      wland_snow[iveg][band].depth = 0.0;
      wland_snow[iveg][band].surf_water = 0.0;
      wland_snow[iveg][band].pack_water = 0.0;
      wland_snow[iveg][band].albedo = 0.0;
      wland_snow[iveg][band].coldcontent = 0.0;
      wland_snow[iveg][band].surf_temp = 0.0;
      wland_snow[iveg][band].pack_temp = 0.0;
    }
    if(*fraci > 0.0) {
      lake_snow->swq = lake->swe/(*fraci);
      lake_snow->surf_water = lake->surf_water/(*fraci);
      lake_snow->pack_water = lake->pack_water/(*fraci);
      lake_snow->depth = lake->sdepth/(*fraci);
      lake_snow->albedo = lake->SAlbedo;
      lake_snow->coldcontent = lake->coldcontent;
      lake_snow->surf_temp = lake->surf_temp;
      lake_snow->pack_temp = lake->pack_temp;
    }
    else {
      lake_snow->swq = 0.0;
      lake_snow->surf_water = 0.0;
      lake_snow->pack_water = 0.0;
      lake_snow->depth = 0.0;
      lake_snow->albedo = 0.0;
      lake_snow->coldcontent = 0.0;
      lake_snow->surf_temp = 0.0;
      lake_snow->pack_temp = 0.0;
    }   
  }
  else {
    if(*fraci > 0.0) {
      lake_snow->swq = lake->swe/(*fraci);
      lake_snow->surf_water = lake->surf_water/(*fraci);
      lake_snow->pack_water = lake->pack_water/(*fraci);
      lake_snow->depth = lake->sdepth/(*fraci);
      lake_snow->albedo = lake->SAlbedo; 
      lake_snow->coldcontent = lake->coldcontent;
      lake_snow->surf_temp = lake->surf_temp;
      lake_snow->pack_temp = lake->pack_temp;
    }
    else {
      lake_snow->swq = 0.0;
      lake_snow->surf_water = 0.0;
      lake_snow->pack_water = 0.0;
      lake_snow->depth = 0.0; 
      lake_snow->albedo = 0.0; 
      lake_snow->coldcontent = 0.0;
      lake_snow->surf_temp = 0.0;
      lake_snow->pack_temp = 0.0;
    }
    wland_snow[iveg][band].swq =0.0;
    wland_snow[iveg][band].depth =0.0;
    wland_snow[iveg][band].surf_water =0.0;
    wland_snow[iveg][band].pack_water =0.0;
    wland_snow[iveg][band].albedo =0.0;
    wland_snow[iveg][band].coldcontent = 0.0;
    wland_snow[iveg][band].surf_temp = 0.0;
    wland_snow[iveg][band].pack_temp = 0.0;
  }
  if (lake_snow->swq == 0.0) {
    lake_snow->MELTING = FALSE;
  }

  if(*lakefrac < 1.0) {
    /* Update canopy storage terms. */
    veg_var[DRY][iveg][band].Wdew /= (1.-(*lakefrac));
    veg_var[WET][iveg][band].Wdew /= (1.-(*lakefrac));
  }

  /* Update Soil Moisture and Ice Content. */
  if (*lakefrac < 1.0) {
    for(i=0; i<options.Nlayer; i++) {
#if EXCESS_ICE
      cell[DRY][iveg][band].layer[i].moist -= soil_con.effective_porosity[i]*soil_con.depth[i]*1000.*(*lakefrac);
      cell[WET][iveg][band].layer[i].moist -= soil_con.effective_porosity[i]*soil_con.depth[i]*1000.*(*lakefrac);
#else
      cell[DRY][iveg][band].layer[i].moist -= soil_con.porosity[i]*soil_con.depth[i]*1000.*(*lakefrac);
      cell[WET][iveg][band].layer[i].moist -= soil_con.porosity[i]*soil_con.depth[i]*1000.*(*lakefrac);
#endif
      cell[DRY][iveg][band].layer[i].moist /= (1.-(*lakefrac));
      cell[WET][iveg][band].layer[i].moist /= (1.-(*lakefrac));
//#if SPATIAL_FROST	
//      for ( frost_area = 0; frost_area < FROST_SUBAREAS; frost_area++ ) {
//        cell[DRY][iveg][band].layer[i].ice[frost_area] /= (1.-(*lakefrac));
//        cell[WET][iveg][band].layer[i].ice[frost_area] /= (1.-(*lakefrac));
//      }
//#else	
//      cell[DRY][iveg][band].layer[i].ice /= (1.-(*lakefrac));
//      cell[WET][iveg][band].layer[i].ice /= (1.-(*lakefrac));
//#endif // SPATIAL_FROST
    }
  }  

  for(index=0; index<options.Nlayer; index++)
    moist[index] = cell[WET][iveg][band].layer[index].moist;
  ErrorFlag = distribute_node_moisture_properties(wland_energy[iveg][band].moist, wland_energy[iveg][band].ice,
                                      wland_energy[iveg][band].kappa_node, wland_energy[iveg][band].Cs_node,
                                      soil_con.dz_node, soil_con.Zsum_node, wland_energy[iveg][band].T,
                                      soil_con.max_moist_node,
#if QUICK_FS
                                      soil_con.ufwc_table_node,
#else
                                      soil_con.expt_node,
                                      soil_con.bubble_node,
#endif // QUICK_FS
#if EXCESS_ICE
                                      soil_con.porosity_node,
                                      soil_con.effective_porosity_node,
#endif // EXCESS_ICE
                                      moist, soil_con.depth,
                                      soil_con.soil_density,
                                      soil_con.bulk_density,
                                      soil_con.quartz, options.Nnode,
                                      options.Nlayer, soil_con.FS_ACTIVE);
  if (ErrorFlag == ERROR) return(ERROR);

  finstor = wland_snow[iveg][band].swq*(1.-(*lakefrac)) + lake_snow->swq*(*lakefrac)*(*fraci);
  finstor += veg_var[WET][iveg][band].Wdew*(1.-(*lakefrac));
  for(i=0; i<options.Nlayer; i++) {
#if EXCESS_ICE
    finstor += cell[WET][iveg][band].layer[i].moist*(1.-(*lakefrac)) + (*lakefrac)*soil_con.effective_porosity[i]*soil_con.depth[i]*1000.;
#else
    finstor += cell[WET][iveg][band].layer[i].moist*(1.-(*lakefrac)) + (*lakefrac)*soil_con.porosity[i]*soil_con.depth[i]*1000.;
#endif // EXCESS_ICE
  }
 
  return(0);

}
