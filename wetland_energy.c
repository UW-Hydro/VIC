#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vicNl.h>

static char vcid[] = "$Id$";

#if LAKE_MODEL

void wetland_energy(int                  rec,
		    atmos_data_struct   *atmos,
		    dist_prcp_struct    *prcp,
		    dmy_struct          *dmy,
		    global_param_struct *gp,
		    soil_con_struct     *soil_con,
		    int                  iveg,
		    int                  band,
		    double               lake_frac,
		    lake_con_struct      lake_con)

/**********************************************************************
	wetland_energy	Laura Bowling		May 12, 2002


**********************************************************************/
{
  extern veg_lib_struct *veg_lib;
  extern option_struct   options;
#if LINK_DEBUG
  extern debug_struct    debug;
#endif
  char                   overstory;
  char                   SOLVE_SURF_ENERGY;
  int                    i, j, k;
  int                    Ndist;
  int                    dist;
  int                    Nveg;
  int                    veg_class;
  int                    Nbands;
  int                    hour;
  double                 rad_atten;
  double                 LAI;
  double                 out_prec[2*MAX_BANDS];
  double                 tmp_surf_temp;
  double                 last_T1;
  double                 out_short=0;
  double                 inshort;
  double                 inlong;
  double                 dp;
  double                 ice0;
  double                 moist;
  double                 surf_atten;
  double                 Tsurf;
  double                 Tgrnd;
  double                 Tend_surf;
  double                 Tend_grnd;
  double                 rainfall[2]; 
  double                 wind_h;
  double                 height;
  double                 displacement[3];
  double                 roughness[3];
  double                 ref_height[3];
  double                 Cv;
  double                 Le;
  double                 Evap;
  double                 Melt[2*MAX_BANDS];
  double                 bare_albedo;
  double                 snow_inflow[MAX_BANDS];
  double                 step_rad;
  double                 step_net_short;
  double                 tmp_aero_resist;
  double                 tmp_throughfall[2][MAX_BANDS];
  double                 tmp_wind[3];
  double                 tmp_melt[MAX_BANDS*2];
  double                 tmp_vapor_flux[MAX_BANDS];
  double                 tmp_canopy_vapor_flux[MAX_BANDS];
  double                 tmp_canopyevap[2][MAX_BANDS];
  double                 tmp_snow_energy;
  double                 tmp_Wdew[2];
  double                 tmp_mu;
  double                 tmp_layerevap[2][MAX_BANDS][MAX_LAYERS];
  double                 tmp_Tmin;
  float root[MAX_LAYERS];
  double                 gauge_correction[2];
  veg_var_struct         tmp_veg_var[2];
  cell_data_struct    ***cell;
  veg_var_struct      ***veg_var;
  energy_bal_struct    **energy;
  energy_bal_struct     *ptr_energy;
  snow_data_struct     **snow;
  snow_data_struct      *tmp_snow;
  veg_var_struct        *tmp_veg[2];
  veg_var_struct        *wet_veg_var;
  veg_var_struct        *dry_veg_var;
  veg_var_struct         empty_veg_var;

  /* set local pointers */
  cell    = prcp->cell;
  energy  = prcp->energy;
  snow    = prcp->snow;
  veg_var = prcp->veg_var;

  /* set variables for distributed precipitation */
  if(options.DIST_PRCP) 
    Ndist = 2;
  else 
    Ndist = 1;

  /* No snow bands in lake model. */
  Nbands = 1;

  /* Set number of vegetation types */
  Nveg      =  iveg-1;

  /** Set Damping Depth **/
  dp        = soil_con->dp;

  /* Compute gauge undercatch correction factors 
     - this assumes that the gauge is free of vegetation effects, so gauge
     correction is constant for the entire grid cell */
  if( options.CORRPREC && atmos->prec[NR] > 0 ) 
    correct_precip(gauge_correction, atmos->wind[NR], gp->wind_h, 
		   soil_con->rough, soil_con->snow_rough);
  else {
    gauge_correction[0] = 1;
    gauge_correction[1] = 1;
  }

  /** Solve Wetland only if Lake Coverage Less than 100% **/
  if (lake_frac  < 0.999) {

    /* fraction of lake area. */
    Cv = 1. - lake_frac;

    /**************************************************
        Initialize Model Parameters
    **************************************************/
      
    /* Initialize energy balance variables */
    energy[iveg][band].shortwave = 0;
    energy[iveg][band].longwave  = 0.;

    /* Initialize snow variables */
    snow[iveg][band].vapor_flux        = 0.;
    snow[iveg][band].canopy_vapor_flux = 0.;
    snow_inflow[band]                  = 0.;
    Melt[band*2]                       = 0.;
  
    /* Initialize precipitation storage */
    for ( j = 0; j < 2; j++ ) out_prec[j] = 0;

    wind_h = gp->wind_h;

    /** Compute Surface Attenuation due to Vegetation Coverage **/
    rad_atten = 0.4;
    LAI = 3.;
    surf_atten = exp(-rad_atten * LAI);

        
    /* Initialize soil thermal properties for the top two layers */
    if(options.FULL_ENERGY || options.FROZEN_SOIL) {
      prepare_full_energy(iveg, Nveg, options.Nnode, prcp, 
			  soil_con, &moist, &ice0);
    }

    /** Compute Bare Soil (free of snow) Albedo **/
    bare_albedo = 0.2;       
    overstory = FALSE;
  
    /*************************************
	Compute the aerodynamic resistance 
    *************************************/

    displacement[0] = .218;
    roughness[0]    = .04;

    /* Initialize wind speeds */
    tmp_wind[0] = atmos->wind[NR];
    tmp_wind[1] = -999.;
    tmp_wind[2] = -999.;
 
    /* Estimate vegetation height */
    height = calc_veg_height(displacement[0]);

    /* Estimate reference height */
    if(displacement[0] < wind_h) 
      ref_height[0] = wind_h;
    else 
      ref_height[0] = displacement[0] + wind_h + roughness[0];
    
    /* Compute aerodynamic resistance over various surface types */
    CalcAerodynamic(overstory, height, 0.2, 
		    soil_con->snow_rough, soil_con->rough, 
		    .4, gp->wind_h, cell[WET][iveg][0].aero_resist, tmp_wind, 
		    displacement, ref_height, roughness, 
		    Nveg, iveg);

    /******************************
        Solve ground surface fluxes 
    ******************************/
    
    wet_veg_var = &(veg_var[WET][iveg][band]);
    dry_veg_var = &(veg_var[DRY][iveg][band]);

    for(i=0; i<MAX_LAYERS; i++)
      root[i] = 0.0;
    root[0] = 1.;

    veg_class=0;

    surface_fluxes(overstory, bare_albedo, height, ice0, moist, 
		   prcp->mu[iveg], surf_atten, &(Melt[band*2]), &Le, 
		   cell[WET][iveg][0].aero_resist, 
		   &(cell[DRY][iveg][band].baseflow), 
		   &(cell[WET][iveg][band].baseflow), displacement, 
		   gauge_correction, &(cell[DRY][iveg][band].inflow), 
		   &(cell[WET][iveg][band].inflow), &out_prec[band*2], 
		   ref_height, roughness, 
		   &(cell[DRY][iveg][band].runoff), 
		   &(cell[WET][iveg][band].runoff), &snow_inflow[band], 
		   tmp_wind, root, Nbands, Ndist, 
		   options.Nlayer, Nveg, band, dp, iveg, rec, veg_class, 
		   atmos, dmy, &(energy[iveg][band]), gp, 
		   cell[DRY][iveg][band].layer, 
		   cell[WET][iveg][band].layer, &(snow[iveg][band]), 
		   soil_con, dry_veg_var, wet_veg_var, .8,
		   .0001, 500.);

    atmos->out_prec += out_prec[band*2] * Cv * lake_con.Cl[0];
  } /** end current vegetation type **/
  
}

#endif // LAKE_MODEL

