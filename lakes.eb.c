#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vicNl.h>

#if LAKE_MODEL

void lakemain(atmos_data_struct  atmos, 
	      lake_var_struct   *lake, 
	      lake_con_struct    lake_con,
	      double             prec,
	      soil_con_struct    soil_con, 
	      int                dt, 
	      energy_bal_struct *lake_energy, 
	      snow_data_struct  *lake_snow, 
	      int                hour, 
	      int                rec, 
	      double             wind_h)
{
  /**********************************************************************
      lakes           Valentijn Pauwels             September 1998

  This subroutine sets up the boundary conditions for simulations with 
  the lake module and serves as the primary driver for the lake model.
  (Based on Hostetler and Bartlein, WRR, 26 (10), 2603-2612, 10/1990)
 
  Modifications:
  12/2000 Model was converted to C and incorporated into VIC.  LCB
  
  Parameters :

  soil_con.lat  	Latitude of the lake (degrees).
  soil_con.lng 	        Longitude of the lake (degrees).
  atmos.air_temp	Air temperature (C).
  atmos.vp	Air vapor pressure (Pa).
  atmos.vpd	Vapor pressure deficit (PA).
  atmos.pressure      	Air pressure (kPa).
  atmos.wind		Wind speed (m/s).
  atmos.longwave	Downwelling long wave radiation smixdepth(W/m2).
  atmos.shortwave       Incoming short wave radiation (W/m2).
  atmos.prec	        Precipitation (mm).
  atmos.density         Air density (kg/m^3).
  lake_energy.latent	Latent heat flux (W/m2).
  lake_energy.sensible	Sensible heat flux (W/m2).
  lake_con.eta_a	Decline of solar radiation input with depth (m-1).  
  lake_con.surface[numnod]	Area of the lake (m2).
  lake.temp[numnod]	Temperature of the lake water (K).
  lake.tempi	         Temperature of the lake ice (K).
  lake.hice           	Height of the lake ice (m).
  lake.fraci	         Fractional coverage of ice (-).
  lake_snow.sdepth	         Height of the snow on top of the lake (m).
  lake_snow.swe              SWE of snow on top of the lake (m).
  lake.numnod        	Number of nodes in the lake (-).
  global_param.dt		         Time step size (hrs).
**********************************************************************/

  /**********************************************************************
   * Solve the energy budget for the lake.
   **********************************************************************/

  solve_lake(prec, atmos.air_temp[hour], atmos.wind[hour],
	     atmos.vp[hour] / 1000., atmos.shortwave[hour], 
	     atmos.longwave[hour], atmos.vpd[hour] / 1000., 
	     atmos.pressure[hour] / 1000., atmos.density[hour], 
	     lake, lake_con, soil_con, dt, rec, lake_energy, lake_snow, 
	     wind_h);

  /**********************************************************************
   * Solve the water budget for the lake.
   **********************************************************************/

  water_balance(lake, lake_con, dt);

}  /*End of  lake main. */


void solve_lake(double             prec, 
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
		double             wind_h)
{

/**********************************************************************
* This subroutine solves the energy budget for open water bodies.
*
* Paremeter description :
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
**********************************************************************/

      double LWnetw,LWneti;
      double sw_water, sw_ice;
      double T[MAX_LAKE_NODES];   /* temp of the water column, open fraction. */
      double Tnew[MAX_LAKE_NODES];  
      double Ti[MAX_LAKE_NODES];   /* temp of the water column, ice fraction. */
      double water_density[MAX_LAKE_NODES], water_cp[MAX_LAKE_NODES];
      float albs, albi, albw;
      double Ts;
      double Tcutoff, Tcutk;  /* Lake freezing temperature (K). */
      double Qhw, Qhi;
      double Qgw, Qgi;
      double Qew, Qei;
      double eflux, eadd;
      double Tki, Tkw;     /* Surface temp. of ice and water in Kelvin. */
      double qbot, qw;
      int freezeflag;
      double mixmax;
      int mixdepth;
      double fracprv; /* Ice coverage at beginning of the timestep. */
      double fracice;
      double new_ice_fraction; /* Ice fraction formed by freezing. */
      int k, i;
      double tw1, tw2, r1, r2, rtu, rnet;
      double Ls, Le;
      double sumjoula, sumjoulb, sumjouli;
      double temp;
      double water_energy;
      double temphw, temphi;
      double SW_under_ice;
      double newsnow;
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
      double sw_out;

      /**********************************************************************
       * 1. Initialize and read in info from previous dt.
       **********************************************************************/
    
      lake_energy->refreeze_energy=0.0;
      lake_energy->deltaCC = 0.0;
      lake_energy->grnd_flux = 0.0;
      lake->snowmlt = 0.0;
      qbot=qw=0.0;
      fracprv=lake->fraci;
      new_ice_height = new_ice_fraction = 0.0;
      lake->evapw=0.0;
      energy_ice_formation = 0.0;
      lake_snow->vapor_flux=0.0;
      
      /* --------------------------------------------------------------------
       * Calculate the water freezing point.
       * -------------------------------------------------------------------- */
     
      rhoinit(&Tcutoff, pressure);          /* degrees C */

      /* --------------------------------------------------------------------
       * Initialize liquid water and ice temperature profiles to be used
       * later in the module.  */
   
      for ( k = 0; k < lake->activenod; k++ )
	{
	  T[k] = lake->temp[k];
	  Ti[k] = lake->temp[k];
	  water_density[k] = calc_density(T[k]);
	  water_cp[k] = specheat(T[k]);
	}

      energycalc( lake->temp, &sumjoulb, lake->activenod, lake_con.dz, 
		  lake->surface, water_cp, water_density);

      /**********************************************************************
       * 2. Calculate added precipitation and total snow height.
       **********************************************************************/

      /* --------------------------------------------------------------------
       * If the air temperature is colder than 2.2 C and there is ice on
       * the lake than the precipitation fallen is treated as snow and
       * added to the snow layer. - THIS IS DIFFERENT FROM REST OF VIC!!
       * ----------------------------------------------------------------- */

      if ( lake->fraci > 0.) {
	if ( tair < -1.5 ) {          /* Precip is falling as snow. */
	  snowfall= (prec);
	  rainfall = 0.0;
	  lake->volume = ( lake->volume + (prec/1000.)*(1-lake->fraci)
			   * lake->surface[0] );
	  newsnow = (prec/1000.);
	}
	else if(tair > 0.5){ /* Precip is falling as rain. */
	  newsnow=0.0;
	  snowfall = 0.0;
	  if(lake_snow->swq > 0.0) {    
	    rainfall = prec;
	    lake->volume += (prec/1000.)*(1-lake->fraci)*lake->surface[0];
	  }
	  else   { /* Assume rain falling on ice drains to lake. */
	    lake->volume += (prec/1000.)*lake->surface[0];
	    rainfall = 0.0; /* Because do not want it added to snow->surf_water */
	  }
	}
	else {
	  rainfall = prec*(tair - (-1.5))/(.5+1.5);
	  snowfall = prec - rainfall;
	  if(lake_snow->swq <=0.0) {
	    lake->volume += (rainfall/1000.)*(lake->fraci)*lake->surface[0];
	    rainfall = 0.0;
	  }
	  lake->volume += (prec/1000.)*(1-lake->fraci)*lake->surface[0];
	}
      }
      else { /* To preserve water balance, lake depth must be increased. */
	lake->volume += (prec/1000.)*lake->surface[0];
	rainfall = 0.0;
	snowfall = 0.0;
	newsnow  = 0.;
      }
    

      /**********************************************************************
       * 3. Calculate incoming solar radiation over water and ice.
       **********************************************************************/

      /* --------------------------------------------------------------------
       * Calculate the albedo of the lake depending on the solar zenith
       * angle for ice, snow and liquid water.  NOT TRUE, CONSTANT VALUES!
       * -------------------------------------------------------------------- */

      alblake(Tcutoff, tair, &albs, &albi, &albw, newsnow, 
	      lake_snow->coldcontent, dt, &lake_snow->last_snow, 
	      lake_snow->swq );
      lake_energy->AlbedoLake = (fracprv) * albi + (1. - fracprv) * albw;
      
      /* --------------------------------------------------------------------
       * Calculate the incoming solar radiaton for both the ice fraction
       * and the liquid water fraction of the lake.
       * -------------------------------------------------------------------- */
      lake_snow->albedo= albs;
      if (lake_snow->swq > SNOWCRIT *RHOSNOW/RHO_W)
	sw_ice = shortin * (1.-albs);
      else if (lake_snow->swq > 0. && lake_snow->swq <= SNOWCRIT*RHOSNOW/RHO_W) {
	albi=(albi+albs)/2.;
	sw_ice = shortin * (1.-albi);
      }
      else if (lake->hice > 0. && lake_snow->swq <= 0.) {
         sw_ice = shortin * (1.-albi);
      }
      else /* lake->hice <= 0 */
	sw_ice=0.0;

      sw_water = shortin * (1.-albw);

      /**********************************************************************
       * 4. Calculate initial energy balance over water.
       **********************************************************************/

      if(fracprv < .999999) {
	freezeflag=1;         /* Calculation for water, not ice. */
	windw = wind*log((2. + ZWATER)/ZWATER)/log(wind_h/ZWATER);  
	
	water_energy_balance( lake->activenod, lake->surface, &lake->evapw, 
			      dt, freezeflag, lake_con.dz, 
			      (double)soil_con.lat, Tcutoff, tair, windw, 
			      pressure, vp, air_density, longin, sw_water, 
			      sumjoulb, wind_h, &Qhw, &Qew, &LWnetw, T, 
			      water_density, &lake_energy->deltaH, 
			      &energy_ice_formation, fracprv, 
			      &new_ice_fraction, water_cp,  &new_ice_height);

	mixmax = 0.0;	 
	/* --------------------------------------------------------------------
	 * Do the convective mixing of the lake water.
	 * -------------------------------------------------------------------- */
	 
	mixdepth=0;        /* Set to zero for this time step. */
	tracer_mixer( T, &mixdepth, freezeflag, lake->surface, 
		      lake->activenod, lake_con.dz, water_cp );
    
	if ( mixdepth > mixmax ) mixmax = mixdepth;

      }           /* End of water fraction calculations. */
      else {
	// ice coveras 100% of lake, reset open water fluxes
	LWnetw = 0;
	Qew    = 0.;
	Qhw    = 0;
      }
 
 
      /**********************************************************************
       *  Calculate initial energy balance over ice.
       **********************************************************************/
 
      if ( fracprv > 0. ) {
	
      	freezeflag = 0;         /* Calculation for ice. */
      	Le = (677. - 0.07 * tair) * JOULESPCAL * GRAMSPKG; /* ice*/
      	windi = ( wind * log((2. + soil_con.snow_rough) / soil_con.snow_rough) 
		  / log(wind_h/soil_con.snow_rough) );
	if ( windi < 1.0 ) windi = 1.0;
	lake->aero_resist = (log((2. + soil_con.snow_rough) 
				 / soil_con.snow_rough) * 
			     log(wind_h/soil_con.snow_rough) / 
			     (von_K*von_K)) / windi;

	/* Calculate snow/ice temperature and change in ice thickness from 
	   surface melting. */
	ice_melt( wind_h+soil_con.snow_rough, lake->aero_resist, Le, 
		  lake_snow, lake, dt,  0.0, soil_con.snow_rough, 1.0, 
		  rainfall, snowfall,  windi, Tcutoff, tair, sw_ice, 
		  longin, air_density, pressure,  vpd,  vp, &lake->snowmlt, 
		  &lake_energy->advection, &lake_energy->deltaCC, 
		  &lake_energy->snow_flux,&Qei, &Qhi, &Qnet_ice, 
		  &lake_energy->refreeze_energy, &LWneti, fracprv);

	lake->tempi = lake_snow->surf_temp;  
	
	/**********************************************************************
	 *  Adjust temperatures of water column in ice fraction.
	 **********************************************************************/

	/* --------------------------------------------------------------------
	 * Calculate inputs to temp_area..
	 * -------------------------------------------------------------------- */

	water_under_ice( freezeflag, sw_ice, wind, Ti, water_density, 
			 (double)soil_con.lat, lake->activenod, lake_con.dz, 
			 Tcutoff, &qw, lake->surface, &temphi, water_cp, 
			 mixdepth, lake->hice, lake_snow->swq*RHO_W/RHOSNOW,
			(double)dt);     
	
	/**********************************************************************
	 *    Calculate change in ice thickness and fraction
	 *    within fraction that already has ice.
	 **********************************************************************/
	if(lake->hice > 0.0) {
	  lakeice(&lake->tempi, Tcutoff, sw_ice, &lake->hice, Ti[0], 
		  &lake->fraci, dt, lake_energy->snow_flux, qw, 
		  &energy_ice_melt_bot, lake_snow->swq*RHO_W/RHOSNOW, temphi);
	}
      }
      else {
	/* No Lake Ice */
	LWneti = 0;
	Qei    = 0.;
	Qhi    = 0;
      }
      
      /**********************************************************************
       * 11. Average ice and water columns.
       **********************************************************************/

      colavg( lake->temp, T, Ti, fracprv,lake->density, lake->activenod, 
	      lake_con.dz);

      /**********************************************************************
       * 12. Calculate the final water heat content and energy balance.
       **********************************************************************/

      lake_energy->AtmosLatent      = ( 1. - fracprv ) * Qew + fracprv * Qei; 
      lake_energy->refreeze_energy *= fracprv;
      lake_energy->refreeze_energy += energy_ice_formation*(1.-fracprv);
      lake_energy->advection       *= fracprv;
     
      lake_energy->AtmosSensible    = ( 1. - fracprv ) * Qhw + fracprv * Qhi;
      lake_energy->NetLongAtmos     = ( ( 1. - fracprv ) * LWnetw + fracprv 
					* LWneti );
      lake_energy->NetShortAtmos    = ( ( 1. - fracprv ) * sw_water 
					+ fracprv * sw_ice );
      lake_energy->snow_flux        = ( lake_energy->snow_flux * fracprv 
					- sw_ice * fracprv ); 
      lake_energy->deltaH          *= ( 1. - fracprv );
      lake_energy->deltaCC         *= ( - fracprv );
 
      lake_energy->error            = ( lake_energy->NetShortAtmos 
					+ lake_energy->NetLongAtmos 
					+ lake_energy->AtmosSensible 
					+ lake_energy->AtmosLatent 
					- lake_energy->deltaCC 
					+ lake_energy->snow_flux
					+ lake_energy->deltaH 
					+ lake_energy->refreeze_energy 
					+ lake_energy->advection );
      
      temphw=temphi=0.0;

      /**********************************************************************
       * 14. In the lake module latent and sensible heat fluxes are positive
       * downwards, to pass out of VIC the fluxes have to be
       * changed sign.  
       **********************************************************************/
       
      /* Adjust water and ice evaporation to be representative of 
	 entire lake. */
      lake->evapw           *= ( (1. - fracprv ) * dt * SECPHOUR ); // in mm
      lake_snow->vapor_flux *= fracprv; // in meters 
      lake->snowmlt         *= fracprv; // in mm 

      // update ice cover 
      if ( new_ice_fraction > 0.0 ) {
	 lake->hice = ( ( new_ice_fraction * new_ice_height 
			  + lake->fraci * lake->hice ) 
			/ ( lake->fraci + new_ice_fraction ) );
	 // remove ice from lake water content
/* 	 lake->volume -= ( new_ice_fraction * new_ice_height  */
/* 			   * lake->surface[0] * ice_density / RHO_W ); */
      }

      lake->fraci += new_ice_fraction;

      if(lake->fraci >= .999999) lake->fraci = 1.0;

      // adjust snowpack for changes in ice cover
      if(lake->fraci < fracprv) {
	lake->snowmlt += ( lake_snow->swq * 1000. 
			   * ( fracprv - lake->fraci ) * lake_con.Cl[0] );
      }
      else {
	if(lake->fraci > fracprv ) {
	  lake_snow->swq        *= fracprv / lake->fraci;
	  lake_snow->surf_water *= fracprv / lake->fraci;
	}
      }
      lake_snow->depth = lake_snow->swq * RHO_W / RHOSNOW;
      
}   /* End of lake function. */

void latsens (double Tsurf, double Tcutk, double hice, 
	       double tair, double wind, double pressure, double vp, double air_density, 
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
	      float  *albs, 
	      float  *albi, 
	      float  *albw, 
	      double  newsnow, 
	      double  qmelts, 
	      int     dt, 
	      int    *last_snow, 
	      double  swq)
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
 **********************************************************************/
  double albgl, albgs;


  if( (Tair-Tcutoff) > 0.0)
    {
      if( (Tair-Tcutoff) < 20.) {
	albgl=0.4 - 0.011*(Tair-Tcutoff);
	albgs=0.6 - 0.0245*(Tair-Tcutoff);
      }
      else {	
	albgl=0.4 - 0.011*20.;
	albgs = 0.6 - 0.0245*20.;
      }
    }
  else
    {
      albgl=0.4;
      albgs=0.6;
    }

  /**albi=0.5*albgs + 0.5*albgl; */
  *albi = 0.25;
  
 /* *albs=0.7; */
  if(newsnow > 0.0)
    *last_snow = 1;
  else if ( swq == 0 )
    (*last_snow) = 0;
  else
    (*last_snow)++;

  *albs = (float)snow_albedo(newsnow, swq, qmelts, dt, *last_snow);
  *albw = 0.05 / 0.15;

}


void colavg (double *finaltemp, double *T,double *Ti,float lakeprv,double *density, int numnod, double dz)
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
	   z=SURF;

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
	   double *de, double lat, int numnod, double dz)
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
      double zhalf[MAX_LAKE_NODES];  /* Originally set to 1000 - why? */
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
 * Calculate the distance between the first node and the surface.
 **********************************************************************/

      zhalf[0] = (SURF + dz) * 0.5;

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
      
	ks=6.6*pow(sin((double)lat*PI/180.),0.5)*pow(wind,-1.84);
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
	  z=SURF+((float)k)*dz;	 
	
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
	      double *fracadd, 
	      int     numnod, 
	      int     dt,  
	      double  dz, 
	      double *cp, 
	      double *surface, 
	      double *new_ice_height, 
	      double *water_density)
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
 **********************************************************************/
     double sum, extra;
     int j;
     double xfrac, di;

/**********************************************************************
 * Calculate the newly added ice to the ice layer.
 **********************************************************************/
     *qfusion = 0.0;
      sum=0.;

      for(j=0; j<numnod;j++) 
	{
	  if (T[j]  < Tcutoff) {

/* --------------------------------------------------------------------
 * Recalculate density and specific heat for frozen water.
 * -------------------------------------------------------------------- */

	    //  water_density = calc_density(T[j]);
           
/* --------------------------------------------------------------------
 * Calculate the ice growth (extra in J).
 * -------------------------------------------------------------------- */
	    
            extra = ( (Tcutoff-T[j]) * dz * (water_density[j]+1000.) * cp[j] 
		      * surface[j] * (1.0-fracprv) );

            if (j == 0) {
	      extra = ( (Tcutoff-T[j]) * SURF * (water_density[j]+1000.) 
			* cp[j] * surface[j] * (1.0-fracprv) );
	    }

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
     
      *qfusion=(sum/(dt*SECPHOUR*surface[0]*(1.0-fracprv))); /*W/m2*/

      di=sum/(Lf*RHOICE);    /* m^3 of ice formed */
      

/**********************************************************************
 * Calculate the added fractional coverage of ice, make sure the
 * total fractional coverage does not exceed 1.
 **********************************************************************/

      *new_ice_height = FRACMIN;
      *fracadd=di/(FRACMIN*surface[0]);

      /* IF entire lake is covered. */
      if ( *fracadd > (1.0 - fracprv)) {
	*new_ice_height = di/((1.-fracprv)*surface[0]);
         *fracadd=1.-fracprv;
      }
     
   }

void icerad (double sw,double hi,double hs,double *avgcond,
	     double *SWnet,double *SW_under_ice)
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
 * SWnet	Net short wave radiation at the top of the lake.
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

      a=(1.-exp(-lamssw*hs))/(CONDS*lamssw);
      b=exp(-lamssw*hs)*(1-exp(-lamisw*hi))/(CONDI*lamisw);
      c=(1.-exp(-lamslw*hs))/(CONDS*lamslw);
      d=exp(-lamslw*hs)*(1-exp(-lamilw*hi))/(CONDI*lamilw);

/* --------------------------------------------------------------------
 * Solar radiation at bottom of snow pack. RHS of Patterson and Hamblin, eq. 7
 * -------------------------------------------------------------------- */

      *SWnet=sw*a1*(a+b)+sw*a2*(c+d);

/* --------------------------------------------------------------------
 * Solar radiation at bottom of snow/ice layer.
 * Patterson and Hamblin eq. 8 (qf-qo)
 * -------------------------------------------------------------------- */

      *SW_under_ice=a1*sw*(1-exp(-(lamssw*hs+lamisw*hi)))
	+a2*sw*(1-exp(-(lamslw*hs+lamilw*hi)));
   }
  
void lakeice (double *tempi, double Tcutoff, double sw_ice, double *hice, 
	      double twater, double *fracice, int dt, double qf,
	      double qw, double *energy_ice_melt_bot, double sdepth,
	      double deltaH)
{
  /**********************************************************************
   * Calculate the growth and decrease in the lake ice cover.
   *
   * Parameters :
   *
   * longin       Downwelling long wave radiation (W/m2).
   * tempi	Temperature of the ice (K).
   * Qhi		Sensible heat flux (W/m2).
   * Qei		Latent heat flux (W/m2).
   * Tcutoff	Temperature at which water freezes (K).
   * sw_ice	Net short wave radiation over ice (W/m2).
   * hice		Ice height (m).
   * swe  	Snow water equivalent (m).
   * twater       Temperature of the top water layer (K).
   * *qbot	Incoming short wave radiation at bottom of the ice (W/m2).
   * *qw		Heat storage in the lake ice (J/m3).
   * evapi	Evaporation from the ice (m/s)..
   * *fracice	Revised fraction of the lake covered by ice due to growth/decrease of exisitng cover (-).
   * dt		Time step (s).
   * ds           Amount of snowmelt (m)
   * qnetice;	Heat flux absorbed into the ice (W/m2). 
   **********************************************************************/
  
  double evaps, evapi;     /* Evaporation from the snow, ice (m). */
  double condqw, evapl; 
  double H;                /* Net radiation for ice layer. */
  double q0;
  double avgcond;          /* Thermal conductivity of the ice and snow pack combined (W/mK). */
  double SWnet;            /* Net short wave radiation at the top of the lake. */
  double SW_under_ice;     /* Incoming short wave radiation at the bottom of the snow-ice layer. */
  double qmelts;           /* Energy available to melt ice. ? */
  double dibot;            /* change in ice surface at the bottom of the pack. */
  double delta_ice_frac;   /* Change in ice covered fraction. */
  double delta_ice_depth;
  double extraf;           /* Extra energy due to the melting of thin ice. */
  double xfrac;
  double tprev;
  double RefrozenWater;
 
  /**********************************************************************
   * Calculate the thermal conductivity of water.
   **********************************************************************/
  
  // condqw = RHO_W * CPW_ICE * SURF / QWTAU;
     
  /**********************************************************************
   * Calculate fluxes at the base of the ice.
   **********************************************************************/
  
  /* --------------------------------------------------------------------
   * Flux of heat in the ice at the ice/water interface (P & H, eq. 8)
   * -------------------------------------------------------------------- */

  //  qf = sw_ice + grnd_flux;
  
  /* --------------------------------------------------------------------
   * Flux of heat from the water to the ice - assumes no through flow. 
   * -------------------------------------------------------------------- */
 
  //qw=-condqw*(Tcutoff-twater)/(SURF/2.);      
  

  //  printf("qw = %f \n", qw);

  /* --------------------------------------------------------------------
   * Amount of heat used to melt the ice (positive means freezing).
   * -------------------------------------------------------------------- */
  
  *energy_ice_melt_bot=qf-qw;
  
  /* --------------------------------------------------------------------
   * Calculate the growth of the ice pack at the bottom (in meters).
   * -------------------------------------------------------------------- */

  dibot=(*energy_ice_melt_bot/(RHOICE*Lf))*dt*SECPHOUR;

  /**********************************************************************
   * Calculate the new height of the ice pack and the fractional ice
   * cover.
   **********************************************************************/

  if (*fracice >= 1.) {    
    /* --------------------------------------------------------------------
     * If the lake is fully covered with ice change only the thickness and
     * not the fractional cover.
 * -------------------------------------------------------------------- */

    *hice+=dibot;

    if (*hice < FRACMIN) {
      /* ....................................................................
       * If the ice height is lower than the minimum ice height increase the
       * height and decrease the fractional cover in order to conserve
       * numerical stability.  Ice covered fraction changes linearly up to FRACMIN.
       * ....................................................................*/

      delta_ice_frac=((FRACMIN - *hice)/FRACMIN)*(*fracice);
      *fracice=1.0-delta_ice_frac;

      if(*fracice > 0.0)
	*hice=FRACMIN;
      else
	*hice=0.0;
    }
  }
  else {

    /* --------------------------------------------------------------------
     * If the lake is not fully covered with ice change both the fractional
     * cover and the ice height.
     * -------------------------------------------------------------------- */

    delta_ice_frac=*fracice*(dibot)/FRACMIN;
    *fracice=*fracice+delta_ice_frac;
         
    if (*fracice > 1.0) {

      /* ....................................................................
       * If the fractioncal cover is larger than 100 % decrease the
       * fractional cover to 100 % and increase the ice depth.
       * ....................................................................*/

      delta_ice_depth=(*fracice - 1.0)*FRACMIN/1.0;
      *hice=*hice+delta_ice_depth;
      *fracice=1.0;
    }

    if ( (*fracice < 0.01) && (delta_ice_frac < 0.) ) {
      
      /* ....................................................................
       * If the ice thickness is lower than the minimal ice thickness
       * and melting has occured set the ice and snow layer thicknesses to zero.
       * ....................................................................*/
      *energy_ice_melt_bot += -1*(*fracice*FRACMIN/1.0)*RHOICE*Lf*(1./(dt*SECPHOUR));
      *fracice=0.0;
      *hice=0.0;
    }
  }
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
	       int dt, double *surface, int numnod, double dz, double *temph, 
	       double *cp, double *energy_out_bottom)
    {
/********************************************************************** 				       
 * Calculate the water temperature for different levels in the lake.
 *
 * Parameters :
 *
 * sw_visible   Shortwave rad in visible band entering top of water column
 * sw_nir       Shortwave rad in near infrared band entering top of water column
 * surface_force The remaining rerms i nthe top layer energy balance 
 * T		Lake water temperature at different levels (K).
 
 * water_density		Water density at different levels (kg/m3).
 * de		Diffusivity of water (or ice) (m2/d).
 * dt		Time step size (s).
 * surface	Area of the lake at different levels (m2).
 * numnod	Number of nodes in the lake (-).
 * dz        Thickness of the lake layers. 
 **********************************************************************/

      double z[MAX_LAKE_NODES], zhalf[MAX_LAKE_NODES];
      double a[MAX_LAKE_NODES], b[MAX_LAKE_NODES], c[MAX_LAKE_NODES];
      double d[MAX_LAKE_NODES];
      double told[MAX_LAKE_NODES];
     
      double dist12;
      int k;
      double surface_1, surface_2, T1;
      double cnextra;
      double swtop; /* The solar radiation at the top of the water column. */
      double top, bot; /* The depth of the top and the bottom of the current water layer. */
      double water_density_new;
      double jouleold, joulenew;
      double energyinput;
      double energymixed;
      double term1, term2;

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
	   z[k] = SURF;
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

/* --------------------------------------------------------------------
 * First calculate d for the surface layer of the lake.
 * -------------------------------------------------------------------- */
      
	
      surface_1 = surface[0];
      surface_2 =(surface[0] + surface[1]) / 2.;

       T1 = (sw_visible*(1*surface_1-surface_2*exp(-lamwsw*SURF)) + 
      	    sw_nir*(1*surface_1-surface_2*exp(-lamwlw*SURF)))/surface[0] 
      	+ (surface_force);          /* W/m2 */

       energyinput +=T1;
       cnextra = 0.5*(surface_2/surface[0])*(de[0]/zhalf[0])*((T[1]-T[0])/z[0]);
      
       energymixed += cnextra;

      d[0]= T[0]+(T1*dt*SECPHOUR)/((1.e3+water_density[0])*cp[0]*z[0])+cnextra*dt*SECPHOUR;
    
/* --------------------------------------------------------------------
 * Calculate d for the remainder of the column.
 * --------------------------------------------------------------------*/

/* ....................................................................
 * All nodes but the deepest node.
 * ....................................................................*/

      for(k=1; k<numnod-1; k++) {
	top = (SURF+(k-1)*dz);
        bot = (SURF+(k)*dz);

         surface_1 =( surface[k-1] + surface[k])/ 2.; 
         surface_2 =( surface[k]  + surface[k+1]) / 2.;


	  T1 = (sw_visible*(surface_1*exp(-lamwsw*top)-surface_2*exp(-lamwsw*bot)) + 
	       sw_nir*(surface_1*exp(-lamwlw*top)-surface_2*exp(-lamwlw*bot)))/surface[k]; 

	  energyinput +=T1;
	  term1 = 0.5 *(1./surface[k])*((de[k]/zhalf[k])*((T[k+1]-T[k])/z[k]))*surface_2;
	  term2 = 0.5 *(-1./surface[k])*((de[k-1]/zhalf[k-1])*((T[k]-T[k-1])/z[k]))*surface_1;
	 
	  cnextra = term1 + term2;
	
	 energymixed += (term1+term2);
         d[k]= T[k]+(T1*dt*SECPHOUR)/((1.e3+water_density[k])*cp[k]*z[k])+cnextra*dt*SECPHOUR;
      }
      
/* ....................................................................
 * Calculation for the deepest node.
 * ....................................................................*/

      surface_1 =( surface[k-1] + surface[k]) / 2.;

      k=numnod-1;
      top = (SURF+(k-1)*dz);
      bot = (SURF+(k)*dz);

       T1 = (sw_visible*(exp(-lamwsw*top)) + 
      	    sw_nir*(exp(-lamwlw*top)))*surface_1/surface[k]; 
      energyinput +=T1;
      cnextra = 0.5 * (-1.*surface_1/surface[k])*((de[k-1]/zhalf[k-1])*((T[k]-T[k-1])/z[k]));
      *energy_out_bottom = (sw_visible*(exp(-lamwsw*bot)) + sw_nir*(exp(-lamwlw*bot)));

      energymixed += cnextra;
    
      d[k] = T[k]+(T1*dt*SECPHOUR)/((1.e3+water_density[k])*cp[k]*z[k])+cnextra*dt*SECPHOUR;


/**********************************************************************
 * Calculate arrays for tridiagonal matrix.
 **********************************************************************/

/* --------------------------------------------------------------------
 * Top node of the column.
 * --------------------------------------------------------------------*/

     
      surface_2 =( surface[0]  + surface[1]) / 2.;

      b[0] = -0.5 * ( de[0] / zhalf[0] ) *
	(  dt*SECPHOUR / z[0] ) * surface_2/surface[0];
      a[0] = 1. - b[0];

/* --------------------------------------------------------------------
 * Second to second last node of the column.
 * --------------------------------------------------------------------*/

      for(k=1;k<numnod-1;k++) {
         surface_1 =( surface[k-1] + surface[k]) / 2.;
         surface_2 =( surface[k]  + surface[k+1]) / 2.;

         b[k] = -0.5 * ( de[k] / zhalf[k] ) *
	   (  dt*SECPHOUR / z[k] )*surface_2/surface[k];
         c[k] = -0.5 * ( de[k-1] / zhalf[k-1] ) *
	   (  dt*SECPHOUR / z[k] )*surface_1/surface[k];
         a[k] = 1. - b[k] - c[k];

	   }
/* --------------------------------------------------------------------
 * Deepest node of the column.
 * --------------------------------------------------------------------*/

      surface_1 =( surface[numnod-1] + surface[numnod-2]) / 2.;
      c[numnod-1] = -0.5 * ( de[numnod-1] / zhalf[numnod-1] ) *
	(  dt*SECPHOUR / z[numnod-1] ) * surface_1/surface[numnod-1];
      a[numnod-1] = 1. - c[numnod-1];

/**********************************************************************
 * Solve the tridiagonal matrix.
 **********************************************************************/
     
      tridia(numnod,c,a,b,d,Tnew);

/**********************************************************************
 * Adjust energy fluxes for change in density -> should be fixed by
 * moving to lagrangian scheme
 **********************************************************************/
    
      energycalc(Tnew, &joulenew, numnod,dz, surface, cp, water_density);
     
      *temph=0.0;
      //for(k=0; k<numnod;k++) {        
      //	water_density_new = calc_density(Tnew[k]);
      //*temph +=z[k]*(water_density[k]-water_density_new)*Tnew[k]*cp[k]*surface[k];
      //}

      *temph = joulenew;

    }

void tracer_mixer (double *T, int *mixdepth,
		   int freezeflag,double *surface,int numnod, double dz, double *cp)
  {
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

    int k,j,m;             /* Counter variables. */
    int mixprev;
    double avet, avev;
    double vol; /* Volume of the surface layer (formerly vol_tr). */
    double heatcon; /*( Heat content of the surface layer (formerly vol). */ 
    double Tav, densnew;
    double rho_max;
    double water_density[MAX_LAKE_NODES];

    for(k=0; k<numnod; k++)
      water_density[k]=calc_density(T[k]);

/**********************************************************************
 * Initialize the top depth of local instability.
 **********************************************************************/

      mixprev = 0;
      
      for(k=0; k<numnod-1; k++) {
/**********************************************************************
 * Check for instability at each slice in water column.
 **********************************************************************/

	avet=0.0;
	avev=0.0;
	
	if (water_density[k] > water_density[k+1]) {

/* --------------------------------------------------------------------
 * If there is instability apply the mixing scheme.
 * --------------------------------------------------------------------*/

	  if (mixprev == 0 && (k+1) > *mixdepth){

/* ....................................................................
 * Correct the top depth of local instability.
 * ....................................................................*/

	    *mixdepth=k+1;
	  }

/*----------------------------------------------------------------------
 * Mix from mixprev to k+1  
 *----------------------------------------------------------------------*/
	   
	  for(m=mixprev; m<=k+1;m++) {

/* --------------------------------------------------------------------
 * Apply the mixing scheme from the previous depth to the instability
 * up to this node.
 * --------------------------------------------------------------------*/

	    if (m == 0) {

/* ....................................................................
 * Calculate the heat content and volume of the surface layer.
 * ....................................................................*/
	      
	      heatcon = SURF*(1.e3+water_density[m])*cp[m]*surface[m];
	      vol = SURF * surface[m];
	    }

	    else
	      {
/* ....................................................................
 * Calculate the heat content and volume of all layers but the surface
 * layer.
 * ....................................................................*/

		heatcon = dz*(1.e3+water_density[m])*cp[m]*surface[m];
		vol = dz * surface[m];
	      }

/* ....................................................................
 * Calculate the volumetric weighed average lake temperature.
 * ....................................................................*/

	    avet=avet+T[m]*heatcon;
	    avev=avev+heatcon;

	  }

	  Tav=avet/avev;

/* --------------------------------------------------------------------
 * Calculate the density of the surface layer.
 * --------------------------------------------------------------------*/

	  densnew = calc_density(Tav);

/* --------------------------------------------------------------------
 * Calculate the maximum density up to the local level.
 * --------------------------------------------------------------------*/

	  rho_max = 0.0;
	  
	  for(j=0;j<mixprev; j++) {
	    
	    if ( (1000.+water_density[j]) > rho_max ) {
	      
	      rho_max=1000.+water_density[j];
	    }
	  }

/* --------------------------------------------------------------------
 * Adjust temperatures and density in the mixed part of column.
 * --------------------------------------------------------------------*/

	  for(j = mixprev;j <= k+1; j++)
	    {
	      T[j]=Tav;
	      water_density[j]=densnew;
	    }

/* --------------------------------------------------------------------
 * Check to make sure that the mixing has not generated new instabilities
 * above the previous depth to the instabilities.
 * --------------------------------------------------------------------*/

	  if (rho_max > (1000.+densnew)) {

/* ....................................................................
 * If there are still instabilities iterate again.
 * ....................................................................*/

	    mixprev = 0;
	   k=-1;
	  }
	}

	else
	  {
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

      for(k=0; k<numnod; k++) {
         water_density[k] = calc_density(T[k]);
      }
  }      

 void tridia (int ne, double *a, double *b, double *c, double *y, double *x)
   {
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

    
 void energycalc (double *finaltemp, double *sumjoule, int numnod, double dz, double *surface, double *cp, double *density)
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
	   energy=finaltemp[k]*SURF*surface[k]*(1000.+density[k])*cp[k];
	 else
	   energy=finaltemp[k]*dz*surface[k]*(1000.+density[k])*cp[k];

	 *sumjoule+=energy;
       }
   }

/**********************************************************************
 * This routine calculates the water balance of the lake
 **********************************************************************/

void water_balance (lake_var_struct *lake, lake_con_struct lake_con, int dt)
{
  int isave_n;
  int i_surf;
  double d_area, d_level, d_volume;
  double runoff_volume;
  double tempvolume,surfacearea,ldepth;
  double remain;
  double i;
  int k;

  /**********************************************************************
   * 1. convert runoff input/output volume to rate
   **********************************************************************/

  isave_n = lake->activenod;   /* save initial depth for later */

  i_surf = (int)floor((lake_con.maxdepth - lake->ldepth)/lake_con.dz); /* basin index of water surface */  

  d_area = lake_con.basin[i_surf] - lake_con.basin[i_surf+1];  /* area slope. */

  lake->sarea = ( ( lake->ldepth - (lake->activenod -.5) * lake_con.dz ) 
		  * (d_area / lake_con.dz) + lake_con.basin[i_surf] ); 

  /* convert from mm over grid cell to m3. */
  runoff_volume = ( ( lake->runoff_in * lake_con.rpercent / 1000. ) 
		    * lake_con.cell_area );  
      
/**********************************************************************
 * 2. calculate change in lake level
 *     grid cell runoff (m3/TS for non-lake area)
 *     snow meltwater (mm)
 *     open water evaporation (mm)
 *     precip (added in solvelake)  
 **********************************************************************/

  lake->volume += ( runoff_volume + ((lake->snowmlt-lake->evapw) 
				     * lake_con.Cl[0] 
				     * lake_con.cell_area / 1000. ) );
  tempvolume = lake->volume;
  ldepth = 0.0;
  for(k=lake_con.numnod -1; k>= 0; k--)
    {
      if(tempvolume >= lake_con.basin[k]*lake_con.dz)
	{
	  tempvolume -= lake_con.basin[k]*lake_con.dz;
	  ldepth += lake_con.dz;
	  surfacearea = lake_con.basin[k];
	}
      else if(tempvolume > 0.0) {
	ldepth += tempvolume/lake_con.basin[k];
	tempvolume=0.0;
	surfacearea=lake_con.basin[k];
      }
    }

/**********************************************************************
 * 3. calculate runoff from lake.  Runoff out is in meters over lake surface.
 * Therefore, maxrate is currently in meters.
 **********************************************************************/

  if(ldepth < lake_con.mindepth)
    lake->runoff_out = 0.0;
  else if(ldepth > lake_con.maxdepth) {
    if((ldepth - lake_con.maxdepth) > lake_con.maxrate)
      lake->runoff_out = (ldepth - lake_con.maxdepth);
    else
      lake->runoff_out = lake_con.maxrate;
  }
  else
    lake->runoff_out = ( lake_con.maxrate * ( ldepth - lake_con.mindepth ) 
			 / ( lake_con.maxdepth - lake_con.mindepth ) );

  // extract runoff volume from the lake
  lake->volume -= ( lake->runoff_out * surfacearea );

  // adjust runoff to lake fractional coverage and mm from m
  lake->runoff_out *= lake_con.Cl[0] * 1000.;

  /* Add in runoff that did not pass through lake. */
  lake->runoff_out += (1. - lake_con.rpercent) * lake->runoff_in;
  
  tempvolume = lake->volume;
  lake->ldepth = 0.0;
  for ( k = lake_con.numnod - 1; k >= 0; k-- ) {
    if(tempvolume >= lake_con.basin[k]*lake_con.dz) {
      tempvolume -= lake_con.basin[k]*lake_con.dz;
      lake->ldepth += lake_con.dz;
      lake->sarea=lake_con.basin[k];
    }
    else if(tempvolume > 0.0) {
      lake->ldepth += tempvolume/lake_con.basin[k];
      tempvolume=0.0;
      lake->sarea=lake_con.basin[k];
    }
  }
  
  if ( lake->ldepth < 0.0 ) lake->ldepth = 0.0;
  
  if (lake->ldepth > (lake_con.maxdepth + 0.001)) {
    printf("volume = %f, depth = %f,  surfacea = %f\n",lake->volume, lake->ldepth, lake->surface[0]);
    printf("runoff_in = %f, runoff_out = %f\n",lake->runoff_in, lake->runoff_out);
    nrerror("Lake is overflowing the basin.");
  }

/**********************************************************************
 *  4. Adjust the activenodes and lake area array. 
 **********************************************************************/

  remain = modf(lake->ldepth/lake_con.dz, &i);   /* remainder over whole dzs */ 
  lake->activenod=i;
  /* Temporary fix. */
  lake->activenod = lake_con.numnod;

  for ( k = 0; k < lake->activenod ; k++ ) {
    lake->surface[k] = lake_con.basin[k + lake_con.numnod - lake->activenod];
  }

  if ( lake->activenod > isave_n ) {
    for ( k = isave_n; k < lake->activenod; k++ )
      lake->temp[k] = lake->temp[isave_n];
  }
  lake->baseflow_out = lake->baseflow_in;
}

#endif // LAKE_MODEL
