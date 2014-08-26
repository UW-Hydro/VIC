#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vicNl.h>

static char vcid[] = "$Id$";

double func_surf_energy_bal(double Ts, va_list ap)
/**********************************************************************
	func_surf_energy_bal	Keith Cherkauer		January 3, 1996

  This subroutine computes the surface energy balance for bare soil
  and vegetation uncovered by snow.  It computes outgoing longwave,
  sensible heat flux, ground heat flux, and storage of heat in the thin
  upper layer, based on the given surface temperature.

  The Energy Balance Equation used comes from Xu Liang's Paper 
  "Insights of the Ground Heat Flux in Land Surface Parameterization
  Schemes."

  Modifications:
  04-14-98 modified to compute evapotranspiration within this routine
           in the hopes of reducing the number of iteration 
	  needed to find a solution surface temperature.       KAC
  07-13-98 modified to include elevation bands for vegetation 
           and snow                                             KAC
  01-20-00 modified to work with the updated radiation estimation
           routines, as well as the simplified frozen soil moisture
           storage                                              KAC
  01-08-01 sensible heat is now set to 0 W/m^2 when the ground is
           fully covered by snow.                                KAC
  04-12-01 fixed error where sensible heat flux for partial bare
           ground was not multiplied by the snow cover fraction. KAC
  04-29-02 moved calculation of sensible heat so that it is computed
           even in water balance mode.  This assures that it is set
           to 0 in water balance mode and does not yield the 
           cumulative sum of sensible heat from the snowpack.    KAC
  11-18-02 modified to compute the effects of blowing snow on the
           surface energy balance.                               LCB
  10-May-04 Added check that both FS_ACTIVE and FROZEN_SOIL are true
	    before computing *fusion.  This is just a safety measure;
	    ice and ice0 should both be 0 if FS_ACTIVE is FALSE.TJB
  16-Jul-04 Renamed VaporMassFlux, BlowingMassFlux, and SurfaceMassFlux
	    to vapor_flux, blowing_flux, and surface_flux, respectively,
	    to denote fact that their units are m/timestep rather than
	    kg/m2s.  Created new variables VaporMassFlux, BlowingMassFlux,
	    and SurfaceMassFlux with units of kg/m2s.  The addresses of
	    the *MassFlux variables are passed to latent_heat_from_snow()
	    where values for the variables are computed.  After these
	    values are computed, vapor_flux, blowing_flux and surface_flux
	    are derived from them by unit conversion.  vapor_flux,
	    blowing_flux, and surface_flux are the variables that are
	    passed in/out of this function.			TJB
  16-Jul-04 Changed the type of the last few variables (lag_one, Nveg,
	    etc) in the va_list to be double.  For some reason, passing
	    them as float or int caused them to become garbage.  This may
	    have to do with the fact that they followed variables of type
	    (double *) in va_list, which may have caused memory alignment
	    problems.						TJB
  05-Aug-04 Removed iveg, LastSnow, dt, SnowDepth, lag_one, sigma_slope,
	    fetch, and Nveg from the this function's argument list,
	    since these variables were only used in the call to
	    latent_heat_from_snow() which no longer needs them.	TJB
  28-Sep-04 Added Ra_used to store the aerodynamic resistance used in
	    flux calculations.					TJB
  2007-Apr-11 Modified to handle grid cell errors by returning to the
              main subroutine, rather than ending the simulation.	GCT
  2007-Apr-24 Removed (1.-snow_coverage) from three equations where it did not
              belong: for calculating LongBareOut in the second two cases and for
              calculating NetBareRad in the third case.			JCA
  2007-Apr-24 Features included for IMPLICIT frozen soils option.	JCA
              (including passing in nrec, nrecs, and iveg)
              (including passing in bulk_density, soil_density, and quartz)
              (including counting cases when IMPLICIT fails and involes EXPLICIT)
  2007-Apr-24 Features included for EXP_TRANS frozen soils option.	JCA
  2007-Apr-24 Passing in Zsum_node.					JCA
  2007-Aug-07 Moved Implicit error counting above call for
              solve_T_profile.						JCA
  2008-Aug-08 Added option for EXCESS_ICE.				JCA
              including: passing entire soil_con structure to
              calc_surf_energy_bal
  2007-Aug-24 Modified to use arno_evap rather than canopy_evap if LAI
              is 0, e.g. winter cropland.				KAC via TJB
  2008-Mar-01 Fixed typo in declaration of ufwc_table_layer.		TJB
  2009-Feb-09 Modified to remove dz_node.				KAC via TJB
  2009-May-20 Corrected the deltaH and fusion terms to account for
	      surf_atten, as in Liang et al, 1999.
	      Added options.GRND_FLUX_TYPE to allow backwards-compatibility
	      with versions 4.0.6 and 4.1.0.				TJB
  2009-Jun-19 Added T flag to indicate whether TFALLBACK occurred.	TJB
  2009-Jul-26 Removed the special logic for the water balance mode, in
	      which net longwave is stored in the "longwave" variable.	TJB
  2009-Sep-19 Added T fbcount to count TFALLBACK occurrences.		TJB
  2009-Nov-15 Changed definitions of D1 and D2 to work for arbitrary
	      node spacing.						TJB
  2010-Apr-24 Replaced ra_under with Ra_used[0].			TJB
  2010-Apr-28 Removed net_short, displacement, roughness, and ref_height
	      from arg list of arno_evap() as they are no longer used.	TJB
  2011-May-31 Removed options.GRND_FLUX.  Soil temperatures and ground
	      flux are now always computed.				TJB
  2011-Jun-03 Added options.ORGANIC_FRACT.  Soil properties now take
	      organic fraction into account.				TJB
  2011-Aug-09 Now method used for estimating soil temperatures depends only
	      on QUICK_FLUX setting.					TJB
  2012-Jan-16 Removed LINK_DEBUG code					BN
  2012-Jan-28 Changed ground flux etc calculations for case of exponential
	      node distribution to be over the same control volume as for
	      the linear node distribution and quick flux cases.	TJB
  2012-Jan-28 Removed AR_COMBO and GF_FULL.				TJB
  2013-Jul-25 Added photosynthesis terms.				TJB
  2013-Dec-26 Removed EXCESS_ICE option.				TJB
  2013-Dec-27 Removed QUICK_FS option.					TJB
  2014-Mar-28 Removed DIST_PRCP option.					TJB
  2014-Apr-25 Added non-climatological veg params.			TJB
  2014-Apr-25 Added partial veg cover fraction, bare soil evap between
	      the plants, and re-scaling of LAI & plant fluxes from
	      global to local and back.					TJB
**********************************************************************/
{
  extern option_struct options;
  extern veg_lib_struct *veg_lib;

  /* define routine input variables */

  /* general model terms */
  int i;
  int rec;
  int nrecs;
  int month;
  int VEG;
  int veg_class;
  int iveg;
  int Error;
  
  //error counting variables for IMPLICIT option
  static int error_cnt0, error_cnt1;  

  double delta_t;

  /* soil layer terms */
  double Cs1;
  double Cs2;
  double D1;
  double D2;
  double T1_old;
  double T2;
  double Ts_old;
  double *Told_node;
  double b_infilt;
  double bubble;
  double dp;
  double expt;
  double ice0;
  double kappa1;
  double kappa2;
  double max_infil;
  double max_moist;
  double moist;

  double *Wmax;
  double *Wcr;
  double *Wpwp;
  double *depth;
  double *resid_moist;
  double *bulk_dens_min;
  double *soil_dens_min;
  double *quartz;
  double *bulk_density;
  double *soil_density;
  double *organic;

  float *root;
  double *CanopLayerBnd;

  /* meteorological forcing terms */
  int UnderStory;
  int overstory;

  double NetShortBare;  // net SW that reaches bare ground
  double NetShortGrnd;  // net SW that penetrates snowpack
  double NetShortSnow;  // net SW that reaches snow surface
  double Tair;          // temperature of canopy air or atmosphere
  double atmos_density;
  double atmos_pressure;
  double elevation;
  double emissivity;
  double LongBareIn; // incoming LW to snow-free surface
  double LongSnowIn; // incoming LW to snow surface - if INCLUDE_SNOW
  double surf_atten;
  double vp;
  double vpd;
  double shortwave;
  double Catm;
  double *dryFrac;

  double *Wdew;
  double *displacement;
  double *ra;
  double *Ra_used;
  double rainfall;
  double *ref_height;
  double *roughness;
  double *wind;
 
  /* latent heat terms */
  double  Le;

  /* snowpack terms */
  double Advection;
  double OldTSurf;
  double TPack;
  double Tsnow_surf;
  double kappa_snow; // snow conductance / depth
  double melt_energy; // energy consumed in reducing the snowpack coverage
  double snow_coverage; // snowpack coverage fraction
  double snow_density;
/*   double snow_depth; */
  double snow_swq;
  double snow_water;

  double *deltaCC;
  double *refreeze_energy;
  double *vapor_flux;
  double *blowing_flux;
  double *surface_flux;

  /* soil node terms */
  int     Nnodes;

  double *Cs_node;
  double *T_node;
  double *Tnew_node;
  char   *Tnew_fbflag;
  int    *Tnew_fbcount;
  double *alpha;
  double *beta;
  double *bubble_node;
  double *Zsum_node;
  double *expt_node;
  double *gamma;
  double *ice_node;
  double *kappa_node;
  double *max_moist_node;
  double *moist_node;

  /* spatial frost terms */
  double *frost_fract;

  /* model structures */
  soil_con_struct *soil_con;
  layer_data_struct *layer;
  veg_var_struct *veg_var;

  /* control flags */
  int INCLUDE_SNOW;
  int FS_ACTIVE;
  int NOFLUX;
  int EXP_TRANS;
  int SNOWING;

  int *FIRST_SOLN;

  /* returned energy balance terms */
  double *NetLongBare; // net LW from snow-free ground
  double *NetLongSnow; // net longwave from snow surface - if INCLUDE_SNOW
  double *T1;
  double *deltaH;
  double *fusion;
  double *grnd_flux;
  double *latent_heat;
  double *latent_heat_sub;
  double *sensible_heat;
  double *snow_flux;
  double *store_error;

  /* Define internal routine variables */
  double Evap;		/** Total evap in m/s **/
  double LongBareOut; // outgoing LW from snow-free ground
  double NetBareRad;
  double TMean;
  double Tmp;
  double error;
  double ice;
/*   double             kappa_snow; */
  double out_long;
  double temp_latent_heat;
  double temp_latent_heat_sub;
  double VaporMassFlux;
  double BlowingMassFlux;
  double SurfaceMassFlux;
  double T11;
  double T1_minus;
  double T1_plus;
  double D1_minus;
  double D1_plus;
  double *transp;
  double Ra_bare[3];
  double tmp_wind[3];
  double tmp_height;
  double tmp_displacement[3];
  double tmp_roughness[3];
  double tmp_ref_height[3];

  /************************************
    Read variables from variable list 
  ************************************/

  /* general model terms */
  rec                     = (int) va_arg(ap, int);
  nrecs                    = (int) va_arg(ap, int);
  month                   = (int) va_arg(ap, int);
  VEG                     = (int) va_arg(ap, int);
  veg_class               = (int) va_arg(ap, int);
  iveg                    = (int) va_arg(ap, int);
  delta_t                 = (double) va_arg(ap, double);

  /* soil layer terms */
  Cs1                     = (double) va_arg(ap, double);
  Cs2                     = (double) va_arg(ap, double);
  D1                      = (double) va_arg(ap, double);
  D2                      = (double) va_arg(ap, double);
  T1_old                  = (double) va_arg(ap, double);
  T2                      = (double) va_arg(ap, double);
  Ts_old                  = (double) va_arg(ap, double);
  Told_node               = (double *) va_arg(ap, double *);
  bubble                  = (double) va_arg(ap, double);
  dp                      = (double) va_arg(ap, double);
  expt                    = (double) va_arg(ap, double);
  ice0                    = (double) va_arg(ap, double);
  kappa1                  = (double) va_arg(ap, double);
  kappa2                  = (double) va_arg(ap, double);
  max_moist               = (double) va_arg(ap, double);
  moist                   = (double) va_arg(ap, double);

  root                    = (float  *) va_arg(ap, float  *);
  CanopLayerBnd           = (double *) va_arg(ap, double *);

  /* meteorological forcing terms */
  UnderStory              = (int) va_arg(ap, int);
  overstory               = (int) va_arg(ap, int);

  NetShortBare            = (double) va_arg(ap, double);
  NetShortGrnd            = (double) va_arg(ap, double);
  NetShortSnow            = (double) va_arg(ap, double);
  Tair                    = (double) va_arg(ap, double);
  atmos_density           = (double) va_arg(ap, double);
  atmos_pressure          = (double) va_arg(ap, double);
  emissivity              = (double) va_arg(ap, double);
  LongBareIn              = (double) va_arg(ap, double);
  LongSnowIn              = (double) va_arg(ap, double);
  surf_atten              = (double) va_arg(ap, double);
  vp                      = (double) va_arg(ap, double);
  vpd                     = (double) va_arg(ap, double);
  shortwave               = (double) va_arg(ap, double);
  Catm                    = (double) va_arg(ap, double);
  dryFrac                 = (double *) va_arg(ap, double *);

  Wdew                    = (double *) va_arg(ap, double *);
  displacement            = (double *) va_arg(ap, double *);
  ra                      = (double *) va_arg(ap, double *);
  Ra_used                 = (double *) va_arg(ap, double *);
  rainfall                = (double) va_arg(ap, double);
  ref_height              = (double *) va_arg(ap, double *);
  roughness               = (double *) va_arg(ap, double *);
  wind                    = (double *) va_arg(ap, double *);

  /* latent heat terms */
  Le                      = (double) va_arg(ap, double);

  /* snowpack terms */
  Advection               = (double) va_arg(ap, double);
  OldTSurf                = (double) va_arg(ap, double);
  TPack                   = (double) va_arg(ap, double);
  Tsnow_surf              = (double) va_arg(ap, double);
  kappa_snow              = (double) va_arg(ap, double);
  melt_energy             = (double) va_arg(ap, double);
  snow_coverage           = (double) va_arg(ap, double);
  snow_density            = (double) va_arg(ap, double);
  snow_swq                = (double) va_arg(ap, double);
  snow_water              = (double) va_arg(ap, double);
    
  deltaCC                 = (double *) va_arg(ap, double *);
  refreeze_energy         = (double *) va_arg(ap, double *);
  vapor_flux              = (double *) va_arg(ap, double *);
  blowing_flux            = (double *) va_arg(ap, double *);
  surface_flux            = (double *) va_arg(ap, double *);

  /* soil node terms */
  Nnodes                  = (int) va_arg(ap, int);

  Cs_node                 = (double *) va_arg(ap, double *);
  T_node                  = (double *) va_arg(ap, double *);
  Tnew_node               = (double *) va_arg(ap, double *);
  Tnew_fbflag             = (char *) va_arg(ap, char *);
  Tnew_fbcount            = (int *) va_arg(ap, int *);
  alpha                   = (double *) va_arg(ap, double *);
  beta                    = (double *) va_arg(ap, double *);
  bubble_node             = (double *) va_arg(ap, double *);
  Zsum_node               = (double *) va_arg(ap, double *);
  expt_node               = (double *) va_arg(ap, double *);
  gamma                   = (double *) va_arg(ap, double *);
  ice_node                = (double *) va_arg(ap, double *);
  kappa_node              = (double *) va_arg(ap, double *);
  max_moist_node          = (double *) va_arg(ap, double *);
  moist_node              = (double *) va_arg(ap, double *);

  /* model structures */
  soil_con                = (soil_con_struct *) va_arg(ap, soil_con_struct *);
  layer               = (layer_data_struct *) va_arg(ap, layer_data_struct *);
  veg_var             = (veg_var_struct *) va_arg(ap, veg_var_struct *);

  /* control flags */
  INCLUDE_SNOW            = (int) va_arg(ap, int);
  NOFLUX                  = (int) va_arg(ap, int);
  EXP_TRANS               = (int) va_arg(ap, int);
  SNOWING                 = (int) va_arg(ap, int);

  FIRST_SOLN              = (int *) va_arg(ap, int *);

  /* returned energy balance terms */
  NetLongBare             = (double *) va_arg(ap, double *);
  NetLongSnow             = (double *) va_arg(ap, double *);
  T1                      = (double *) va_arg(ap, double *);
  deltaH                  = (double *) va_arg(ap, double *);
  fusion                  = (double *) va_arg(ap, double *);
  grnd_flux               = (double *) va_arg(ap, double *);
  latent_heat             = (double *) va_arg(ap, double *);
  latent_heat_sub         = (double *) va_arg(ap, double *);
  sensible_heat           = (double *) va_arg(ap, double *);
  snow_flux               = (double *) va_arg(ap, double *);
  store_error             = (double *) va_arg(ap, double *);

  /* take additional variables from soil_con structure */
  b_infilt = soil_con->b_infilt;
  max_infil = soil_con->max_infil;
  Wmax = soil_con->max_moist;
  Wcr = soil_con->Wcr;
  Wpwp = soil_con->Wpwp;
  depth = soil_con->depth;
  resid_moist = soil_con->resid_moist;
  elevation = (double)soil_con->elevation;
  frost_fract = soil_con->frost_fract;
  FS_ACTIVE = soil_con->FS_ACTIVE;
  /* more soil layer terms for IMPLICIT option*/
  bulk_dens_min = soil_con->bulk_dens_min;
  soil_dens_min = soil_con->soil_dens_min;
  quartz = soil_con->quartz;
  bulk_density = soil_con->bulk_density;
  soil_density = soil_con->soil_density;
  organic = soil_con->organic;


  /***************
    MAIN ROUTINE
  ***************/

  Error = 0;
  if(rec==0){
    error_cnt0=0;
    error_cnt1=0;
  }

  TMean = Ts;
  Tmp = TMean + KELVIN;

  transp = (double *) calloc(options.Nlayer,sizeof(double));
  for (i=0; i<options.Nlayer; i++) {
    transp[i] = 0;
  }

  /**********************************************
    Compute Surface Temperature at Half Time Step
  **********************************************/

  if ( snow_coverage > 0 && !INCLUDE_SNOW ) {

    /****************************************
      Compute energy flux through snow pack
    ****************************************/

    *snow_flux = ( kappa_snow * (Tsnow_surf - TMean) );

  }
  else if ( INCLUDE_SNOW ) {
    *snow_flux = 0;
    Tsnow_surf = TMean;
  }
  else *snow_flux = 0;

  /***************************************************************
    Estimate soil temperatures for ground heat flux calculations
  ***************************************************************/

  if ( options.QUICK_FLUX ) {
    /**************************************************************
      Use Liang et al. 1999 Equations to Calculate Ground Heat Flux
      NOTE: T2 is not the temperature of layer 2, nor of node 2, nor at depth dp;
      T2 is the constant temperature at depths much greater than dp to which
      the soil temperature profile is asymptotic.
    **************************************************************/
    *T1 = estimate_T1(TMean, T1_old, T2, D1, D2, kappa1, kappa2, Cs1, Cs2, dp, delta_t);
    
    /*****************************************************
      Compute the Ground Heat Flux from the Top Soil Layer
    *****************************************************/
    if (options.GRND_FLUX_TYPE == GF_406) {
      *grnd_flux = (snow_coverage + (1. - snow_coverage) * surf_atten) * (kappa1 / D1 * ((*T1) - TMean));
    }
    else {
      *grnd_flux = (snow_coverage + (1. - snow_coverage) * surf_atten) * (kappa1 / D1 * ((*T1) - TMean)
				       + (kappa2 / D2 * ( 1. - exp( -D1 / dp )) * (T2 - (*T1)))) / 2.;
    }
 
  }
  else {
    /*************************************************************
      Use Finite Difference Method to Solve Ground Heat
      Flux at Soil Thermal Nodes (Cherkauer and Lettenmaier, 1999)
    *************************************************************/
    T_node[0] = TMean;
      
    /* IMPLICIT Solution */
    if(options.IMPLICIT) {
      Error = solve_T_profile_implicit(Tnew_node, T_node, Tnew_fbflag, Tnew_fbcount, Zsum_node, kappa_node, Cs_node, 
				       moist_node, delta_t, max_moist_node, bubble_node, expt_node, 
				       ice_node, alpha, beta, gamma, dp, Nnodes, 
				       FIRST_SOLN, FS_ACTIVE, NOFLUX, EXP_TRANS, veg_class,
				       bulk_dens_min, soil_dens_min, quartz, bulk_density, soil_density, organic, depth);
      
      /* print out error information for IMPLICIT solution */
      if(Error==0)
        error_cnt0++;
      else
        error_cnt1++;
      if(FIRST_SOLN[1]){
        FIRST_SOLN[1] = FALSE;
#if VERBOSE
        if ( iveg == 0 && rec == nrecs - 1) 
          fprintf(stderr,"The implicit scheme failed %d instances (%.1f%c of attempts).\n",error_cnt1,100.0*(float)error_cnt1/((float)error_cnt0+(float)error_cnt1),'%');
#endif
      }
    }

    /* EXPLICIT Solution, or if IMPLICIT Solution Failed */
    if(!options.IMPLICIT || Error == 1) {
      if(options.IMPLICIT)
        FIRST_SOLN[0] = TRUE;
      Error = solve_T_profile(Tnew_node, T_node, Tnew_fbflag, Tnew_fbcount, Zsum_node, kappa_node, Cs_node, 
			      moist_node, delta_t, max_moist_node, bubble_node, 
			      expt_node, ice_node, alpha, beta, gamma, dp, depth, 
			      Nnodes, FIRST_SOLN, FS_ACTIVE, NOFLUX, EXP_TRANS, veg_class);
    }
      
    if ( (int)Error == ERROR ) {
      fprintf(stderr, "ERROR: func_surf_energy_bal calling solve_T_profile\n");
      return( ERROR ); 
    }
    /* Compute temperatures for calculations of ground heat flux, delta_H, and fusion */
    if (!options.EXP_TRANS) {
      *T1 = Tnew_node[1];
      T11 = Tnew_node[2];
    }
    else {
      i=0;
      while (soil_con->Zsum_node[i] < D1) { i++; }
      D1_minus = soil_con->Zsum_node[i-1];
      D1_plus = soil_con->Zsum_node[i];
      T1_minus = Tnew_node[i-1];
      T1_plus = Tnew_node[i];
      *T1 = T1_minus + (T1_plus - T1_minus)*(D1-D1_minus)/(D1_plus-D1_minus);
    }
      

    /*****************************************************
      Compute the Ground Heat Flux between layers 0 and 1
    *****************************************************/
    if (!options.EXP_TRANS) {
      if (options.GRND_FLUX_TYPE == GF_406) {
        *grnd_flux = (snow_coverage + (1. - snow_coverage) * surf_atten) * (kappa1 / D1 * ((*T1) - TMean));
      }
      else {
        *grnd_flux = (snow_coverage + (1. - snow_coverage) * surf_atten) * (kappa1 / D1 * ((*T1) - TMean) 
		       + (kappa2 / D2 * (T11 - (*T1)))) / 2.;
      }
    }
    else {
        *grnd_flux = (snow_coverage + (1. - snow_coverage) * surf_atten) * (kappa1 / (D1-D1_minus) * ((*T1) - T1_minus) 
		       + (kappa2 / (D1_plus-D1) * (T1_plus - (*T1)))) / 2.;
    }

  }

  /******************************************************
    Compute the change in heat storage in the region between layers 0 and 1
  ******************************************************/
  if (!options.EXP_TRANS) {
    *deltaH = (Cs1 * ((Ts_old + T1_old) - (TMean + *T1)) * D1 / delta_t / 2.);
  }
  else {
    *deltaH = 0;
    i=0;
    while (soil_con->Zsum_node[i+1] < D1) {
      *deltaH += (Cs1 * ((Told_node[i] + Told_node[i+1]) - (Tnew_node[i] + Tnew_node[i+1])) * (soil_con->Zsum_node[i+1]-soil_con->Zsum_node[i]) / delta_t / 2.);
      i++;
    }
    *deltaH += (Cs1 * ((Told_node[i] + T1_old) - (Tnew_node[i] + *T1)) * (D1-soil_con->Zsum_node[i]) / delta_t / 2.);
  }

  /******************************************************
    Compute the change in heat due to solid-liquid phase changes in the region between layers 0 and 1
  ******************************************************/
  if (FS_ACTIVE && options.FROZEN_SOIL) {

    if (!options.EXP_TRANS) {
      if((TMean+ *T1)/2.<0.) {
        ice = moist - maximum_unfrozen_water((TMean+ *T1)/2.,
					     max_moist,bubble,expt);
        if(ice<0.) ice=0.;
      }
      else ice=0.;
      *fusion = (-ice_density * Lf * (ice0 - ice) * D1 / delta_t);
    }
    else {
      *fusion = 0;
      i=0;
      while (soil_con->Zsum_node[i+1] < D1) {
        if((Told_node[i]+Told_node[i+1])/2.<0.) {
          ice0 = moist - maximum_unfrozen_water((Told_node[i]+Told_node[i+1])/2.,
					       max_moist,bubble,expt);
          if(ice0<0.) ice0=0.;
        }
        else ice0=0.;
        if((Tnew_node[i]+Tnew_node[i+1])/2.<0.) {
          ice = moist - maximum_unfrozen_water((Tnew_node[i]+Tnew_node[i+1])/2.,
					       max_moist,bubble,expt);
          if(ice<0.) ice=0.;
        }
        else ice=0.;
        *fusion += (-ice_density * Lf * (ice0 - ice) * (soil_con->Zsum_node[i+1]-soil_con->Zsum_node[i]) / delta_t);
        i++;
      }
      if((Told_node[i]+T1_old)/2.<0.) {
        ice0 = moist - maximum_unfrozen_water((Told_node[i]+T1_old)/2.,
					     max_moist,bubble,expt);
        if(ice0<0.) ice0=0.;
      }
      else ice0=0.;
      if((Tnew_node[i]+ *T1)/2.<0.) {
        ice = moist - maximum_unfrozen_water((Tnew_node[i]+ *T1)/2.,
					     max_moist,bubble,expt);
        if(ice<0.) ice=0.;
      }
      else ice=0.;
      *fusion += (-ice_density * Lf * (ice0 - ice) * (D1-soil_con->Zsum_node[i]) / delta_t);
    }

  }
    
  /* if thin snowpack, compute the change in energy stored in the pack */
  if ( INCLUDE_SNOW ) {
    if ( TMean > 0 )
      *deltaCC = CH_ICE * (snow_swq - snow_water) * (0 - OldTSurf) / delta_t;
    else
      *deltaCC = CH_ICE * (snow_swq - snow_water) * (TMean - OldTSurf) / delta_t;
    *refreeze_energy = (snow_water * Lf * snow_density) / delta_t;
    *deltaCC *= snow_coverage; // adjust for snow cover fraction
    *refreeze_energy *= snow_coverage; // adjust for snow cover fraction
  }
    
  /** Compute net surface radiation of snow-free area for evaporation estimates **/
  LongBareOut = STEFAN_B * Tmp * Tmp * Tmp * Tmp;
  if ( INCLUDE_SNOW ) { // compute net LW at snow surface
    (*NetLongSnow) = (LongSnowIn - snow_coverage * LongBareOut);
  }
  (*NetLongBare)   = (LongBareIn - (1. - snow_coverage) * LongBareOut); // net LW snow-free area
  NetBareRad = (NetShortBare + (*NetLongBare) + *grnd_flux + *deltaH + *fusion);
    
  /** Compute atmospheric stability correction **/
  /** CHECK THAT THIS WORKS FOR ALL SUMMER SITUATIONS **/
  if ( wind[UnderStory] > 0.0 && overstory && SNOWING )
    Ra_used[0] = ra[UnderStory] 
      / StabilityCorrection(ref_height[UnderStory], 0.f, TMean, Tair, 
			    wind[UnderStory], roughness[UnderStory]);
  else if ( wind[UnderStory] > 0.0 )
    Ra_used[0] = ra[UnderStory] 
      / StabilityCorrection(ref_height[UnderStory], displacement[UnderStory], 
			    TMean, Tair, wind[UnderStory], roughness[UnderStory]);
  else
    Ra_used[0] = HUGE_RESIST;
  
  /*************************************************
    Compute Evapotranspiration if not snow covered

    Should evapotranspiration be active when the 
    ground is only partially covered with snow????

    Use Arno Evap if LAI is set to zero (e.g. no
    winter crop planted).
  *************************************************/
  if ( VEG && !SNOWING && veg_var->vegcover > 0 ) {
    Evap = canopy_evap(layer, veg_var, TRUE, 
		       veg_class, month, Wdew, delta_t, NetBareRad, vpd, 
		       NetShortBare, Tair, Ra_used[1], 
		       displacement[1], roughness[1], ref_height[1], 
		       elevation, rainfall, depth, Wmax, Wcr, Wpwp, frost_fract,
		       root, dryFrac, shortwave, Catm, CanopLayerBnd);
    if (veg_var->vegcover < 1) {
      for (i=0; i<options.Nlayer; i++) {
        transp[i] = layer[i].evap;
        layer[i].evap = 0;
      }
      tmp_wind[0] = wind[0];
      tmp_wind[1] = -999.;
      tmp_wind[2] = -999.;
      tmp_height = soil_con->rough/0.123;
      tmp_displacement[0] = calc_veg_displacement(tmp_height);
      tmp_roughness[0] = soil_con->rough;
      tmp_ref_height[0] = 10;
      Error = CalcAerodynamic(0,0,0,soil_con->snow_rough,soil_con->rough,0,Ra_bare,tmp_wind,tmp_displacement,tmp_ref_height,tmp_roughness);
      Ra_bare[0] /= StabilityCorrection(tmp_ref_height[0], tmp_displacement[0], TMean, Tair, tmp_wind[0], tmp_roughness[0]);
      Evap *= veg_var->vegcover;
      Evap += (1-veg_var->vegcover)
	       * arno_evap(layer, surf_atten*NetBareRad, Tair, vpd, 
		       depth[0], max_moist * depth[0] * 1000., 
		       elevation, b_infilt, Ra_bare[0], delta_t, 
		       resid_moist[0], frost_fract);
      for (i=0; i<options.Nlayer; i++) {
        layer[i].evap = veg_var->vegcover*transp[i] + (1-veg_var->vegcover)*layer[i].evap;
        if (layer[i].evap > 0)
          layer[i].bare_evap_frac = 1 - (veg_var->vegcover*transp[i])/layer[i].evap;
        else
          layer[i].bare_evap_frac = 0;
      }
      veg_var->throughfall = (1-veg_var->vegcover)*rainfall + veg_var->vegcover*veg_var->throughfall;
      veg_var->canopyevap *= veg_var->vegcover;
      veg_var->Wdew *= veg_var->vegcover;
    }
    else {
      for (i=0; i<options.Nlayer; i++) {
        layer[i].bare_evap_frac = 0;
      }
    }
  }
  else if(!SNOWING) {
    Evap = arno_evap(layer, NetBareRad, Tair, vpd, 
		     depth[0], max_moist * depth[0] * 1000., 
		     elevation, b_infilt, Ra_used[0], delta_t, 
		     resid_moist[0], frost_fract);
    for (i=0; i<options.Nlayer; i++) {
      layer[i].bare_evap_frac = 1;
    }
  }
  else Evap = 0.;

  free(transp);
  
  /**********************************************************************
    Compute the Latent Heat Flux from the Surface and Covering Vegetation
  **********************************************************************/
  *latent_heat  = -RHO_W * Le * Evap;
  *latent_heat_sub = 0.;

  /** Compute the latent heat flux from a thin snowpack if present **/
  if (INCLUDE_SNOW) {

    /* Convert sublimation terms from m/timestep to kg/m2s */
    VaporMassFlux = *vapor_flux * ice_density / delta_t;
    BlowingMassFlux = *blowing_flux * ice_density / delta_t;
    SurfaceMassFlux = *surface_flux * ice_density / delta_t;

    latent_heat_from_snow(atmos_density, vp, Le, atmos_pressure, 
			  Ra_used[0], TMean, vpd, &temp_latent_heat, 
			  &temp_latent_heat_sub, &VaporMassFlux,
                          &BlowingMassFlux, &SurfaceMassFlux);
    *latent_heat += temp_latent_heat * snow_coverage;
    *latent_heat_sub = temp_latent_heat_sub * snow_coverage;

    /* Convert sublimation terms from kg/m2s to m/timestep */
    *vapor_flux = VaporMassFlux * delta_t / ice_density;
    *blowing_flux = BlowingMassFlux * delta_t / ice_density;
    *surface_flux = SurfaceMassFlux * delta_t / ice_density;

  }
  else *latent_heat *= (1. - snow_coverage);

  /************************************************
    Compute the Sensible Heat Flux from the Surface
  ************************************************/
  if ( snow_coverage < 1 || INCLUDE_SNOW ) {
    *sensible_heat = atmos_density * Cp * (Tair - (TMean)) / Ra_used[0];
    if ( !INCLUDE_SNOW ) (*sensible_heat) *= (1. - snow_coverage);
  }
  else *sensible_heat = 0.;

  /*************************************
    Compute Surface Energy Balance Error
  *************************************/
  error = (NetBareRad // net radiation on snow-free area
	   + NetShortGrnd + NetShortSnow // net surface SW 
	   + emissivity * (*NetLongSnow)) // net surface LW 
	   + *sensible_heat // surface sensible heat
	   + (*latent_heat + *latent_heat_sub) // surface latent heats
	   /* heat flux through snowpack - for snow covered fraction */
	   + *snow_flux * snow_coverage
	   /* energy used in reducing snow coverage area */
	   + melt_energy 
	   /* snow energy terms - values are 0 unless INCLUDE_SNOW */
	   + Advection - *deltaCC;
    
  if ( INCLUDE_SNOW ) {
    if (Tsnow_surf == 0.0 && error > -(*refreeze_energy)) {
      *refreeze_energy = -error;
      error = 0.0;
    }
    else {
      error += *refreeze_energy;
    }
  }

  *store_error = error;

  return error;

}

