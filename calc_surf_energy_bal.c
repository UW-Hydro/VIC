#include <stdio.h>
#include <stdlib.h>
#include "vicNl.h"

double calc_surf_energy_bal(int rec,
			    int iveg,
			    int nlayer,
			    int Nveg,
			    int snow,
			    double ice0,
			    double moist,
			    double dp,
			    double surf_atten,
			    double T0,
			    double *T1,
			    double *aero_resist,
			    double Le,
			    double Ls,
			    double mu,
			    double *x,
			    double displacement,
			    double roughness,
			    double ref_height,
			    double snow_energy,
			    double *wind,
			    double vapor_flux,
			    double rainfall,
			    atmos_data_struct *atmos,
			    veg_var_struct *veg_var,
			    energy_bal_struct *energy,
			    global_param_struct gp,
			    layer_data_struct *layer,
			    soil_con_struct soil_con,
			    int veg_class,
			    dmy_struct dmy)

/**************************************************************
  calc_surf_energy_bal.c  Greg O'Donnell and Keith Cherkauer  Sept 9 1997
  
  This function calculates the surface temperature, in the
  case of no snow cover.  Evaporation is computed using the
  previous ground heat flux, and then used to comput latent heat
  in the energy balance routine.  Surface temperature is found
  using the Root Brent method (Numerical Recipies).
  
  Units:
	layer[].evap 		[mm/time step]
	veg_var.canopyevap 	[mm/time step]
	surf_temp		[C]
	atmos->rad		[W/m^2]
	energy->grnd_flux	[W/m^2]
	energy->sensible	[W/m^2]
	energy->latent		[W/m^2]
  
	rec		N/A	Current Record Number
	iveg		N/A	Current Vegetation Index
	nlayer		N/A	Number of soil layers
	Nveg		N/A	Number of vegetation types in current cell
	snow		N/A	TRUE = Snow was present this time step
	ice0		mm/mm	Ice content of top layer
	moist		mm/mm	Water content of top soil layer
	dp		m	Soil Thermal Damping Depth 
	surf_atten	fract	Attenuation of radiation due to vegetation
	T0		C	Surface Temperature
	*T1		C	Temperature between top two layers
	*aero_resist	s/m	Aerodynamic resistance
	Le		J/kg	Latent heat of evaporation
	Ls		J/kg	Latent heat of sublimation
	mu		fract	Fraction of grid cell receiving precipitation
	*x
	displacement	m	Displacement height of vegetation type
	roughness	m	Roughness length of surface or vegetation
	ref_height	m	Reference height at which wind measured
	snow_energy	W/m^2	Internal energy used to melt snow pack
	*wind		m/s	Wind speed
	vapor_flux	mm/TS	Vapor flux from snow pack
	rainfall	mm/TS	Rain
	*atmos		N/A
	*veg_var	N/A
	*energy		N/A
	gp		N/A
	*layer		N/A
	soil_con	N/A
	veg_class	N/A
	dmy		N/A
  Outstanding:
  emissivity of veg
***************************************************************/
{
  extern veg_lib_struct *veg_lib;
  extern option_struct options;

  double T2;
  double Ts_old;
  double T1_old;
  double Tair;
  double ra;
  double atmos_density;
  double shortwave;
  double longwave;
  double albedo;
  double emissivity;
  double kappa1;
  double kappa2;
  double Cs1;
  double Cs2;
  double D1;
  double D2;
  double delta_t;
  double Vapor;
  double max_moist;
  double bubble;
  double expt;
  double dH_height;

  char   CALC_EVAP;
  char   VEG;
  int    i;
  double surf_temp;
  double U;
  double error;
  double C1, C2, C3;
  double ice;
  double Wdew;

  /**************************************************
    Correct Aerodynamic Resistance for Stability
  **************************************************/
  ra = aero_resist[0];
  U = wind[0];
  if (U > 0.0)
    ra /= StabilityCorrection(ref_height, displacement, T0,
          atmos->air_temp, U, roughness);
  else
    ra = HUGE_RESIST;


  /**************************************************
    Compute Net Surface Radiation
  **************************************************/
  if(iveg!=Nveg) albedo = veg_lib[veg_class].albedo[dmy.month-1];
  else albedo = atmos->albedo;
  energy->albedo = albedo;

  /**************************************************
    Compute Evaporation and Transpiration 
  **************************************************/
  if(iveg!=Nveg) {
    if(veg_lib[veg_class].LAI[dmy.month-1] > 0.0) VEG = TRUE;
    else VEG = FALSE;
  }
  else VEG = FALSE;
  if(!snow) CALC_EVAP=TRUE;
  else CALC_EVAP=FALSE;
  Vapor = vapor_flux / (double)gp.dt / 3600.;

  /**************************************************
    Set All Variables For Use
  **************************************************/
  T2 = energy->T[gp.Ulayer+gp.Llayer+1];
  Ts_old = energy->T[0];
  T1_old = energy->T[1];
  Tair = atmos->air_temp;
  atmos_density = atmos->density;
  shortwave = atmos->shortwave;
  longwave = atmos->longwave;
  emissivity = 1.;
  kappa1 = energy->kappa[0];
  kappa2 = energy->kappa[1];
  Cs1 = energy->Cs[0];
  Cs2 = energy->Cs[1];
  D1 = soil_con.depth[0];
  D2 = soil_con.depth[1];
  delta_t = (double)gp.dt*3600.;
  max_moist = soil_con.max_moist[0]/(soil_con.depth[0]*1000.);
  bubble = soil_con.bubble;
  expt = soil_con.expt[0];
  Wdew = veg_var[0].Wdew;

  /**************************************************
    Find Surface Temperature Using Root Brent Method
  **************************************************/
  surf_temp = root_brent(0.5*(energy->T[0]+atmos->air_temp)+25.,
			 0.5*(energy->T[0]+atmos->air_temp)-25.,
			 func_surf_energy_bal,T2,Ts_old,T1_old,Tair,ra,
			 atmos_density,
			 shortwave,longwave,albedo,emissivity,
			 kappa1,kappa2,Cs1,Cs2,D1,D2,dp,
			 delta_t,Le,Ls,Vapor,moist,ice0,
			 max_moist,bubble,expt,surf_atten,U,displacement,
			 roughness,
			 ref_height,soil_con.elevation,
			 soil_con.b_infilt,soil_con.max_infil,
			 (double)gp.dt,atmos->vpd,rainfall,Wdew,
			 snow_energy,
			 &energy->grnd_flux,T1,
			 &energy->latent,
			 &energy->sensible,&energy->deltaH,
			 &energy->error,&energy->Trad[0],&atmos->rad,
			 soil_con.depth,
			 soil_con.Wcr,soil_con.Wpwp,layer,veg_var,
			 VEG,CALC_EVAP,
			 veg_class,dmy.month);

  if(surf_temp<=-9998) 
    error_calc_surf_energy_bal(surf_temp,T2,Ts_old,T1_old,Tair,ra,
			       atmos_density,shortwave,longwave,albedo,
			       emissivity,kappa1,kappa2,Cs1,Cs2,D1,D2,dp,
			       delta_t,Le,Ls,Vapor,moist,ice0,
			       max_moist,bubble,expt,surf_atten,U,
			       displacement,roughness,
			       ref_height,soil_con.elevation,
			       soil_con.b_infilt,soil_con.max_infil,
			       (double)gp.dt,atmos->vpd,rainfall,Wdew,
			       snow_energy,
			       &energy->grnd_flux,T1,&energy->latent,
			       &energy->sensible,&energy->deltaH,
			       &energy->error,&energy->Trad[0],&atmos->rad,
			       soil_con.depth,
			       soil_con.Wcr,soil_con.Wpwp,layer,veg_var,
			       VEG,CALC_EVAP,
			       veg_class,dmy.month,iveg);

  /**************************************************
    Recalculate Energy Balance Terms for Final Surface Temperature
  **************************************************/
  error = solve_surf_energy_bal(surf_temp,T2,Ts_old,T1_old,Tair,ra,
				atmos_density,shortwave,longwave,albedo,
				emissivity,kappa1,kappa2,Cs1,Cs2,D1,D2,dp,
				delta_t,Le,Ls,Vapor,moist,ice0,
				max_moist,bubble,expt,surf_atten,U,
				displacement,roughness,
				ref_height,soil_con.elevation,
				soil_con.b_infilt,soil_con.max_infil,
				(double)gp.dt,atmos->vpd,rainfall,Wdew,
				snow_energy,
				&energy->grnd_flux,T1,&energy->latent,
				&energy->sensible,&energy->deltaH,
				&energy->error,&energy->Trad[0],&atmos->rad,
				soil_con.depth,
				soil_con.Wcr,soil_con.Wpwp,layer,veg_var,
				VEG,CALC_EVAP,
				veg_class,dmy.month);

  energy->error = error;

  return (surf_temp);
    
}

double solve_surf_energy_bal(double Tsurf, ...) {

  va_list ap;

  double error;

  va_start(ap,Tsurf);
  error = func_surf_energy_bal(Tsurf, ap);
  va_end(ap);

  return error;

}

double error_calc_surf_energy_bal(double Tsurf, ...) {

  va_list ap;

  double error;

  va_start(ap,Tsurf);
  error = error_print_surf_energy_bal(Tsurf, ap);
  va_end(ap);

  return error;

}

double error_print_surf_energy_bal(double Ts, va_list ap) {

  extern option_struct options;

  /** Thermal Properties **/
  double T2;		/** average soil temperature (C) **/
  double Ts_old;	/** last temperature (C) **/
  double T1_old;	/** last layer 1 soil temperature (C) **/
  double Tair;		/** Air Temperature **/
  double ra;		/** aerodynamic reisistance (s/m) **/
  double atmos_density;	/** atmospheric density (kg/m^3) **/
  double shortwave;
  double longwave;
  double albedo;
  double emissivity;
  double kappa1;	/** thermal conductivity of 1st layer */
  double kappa2;	/** thermal conductivity of 2nd layer */
  double Cs1; 		/** volumetric heat capacity of 1st layer **/
  double Cs2; 		/** volumetric heat capacity of 2nd layer **/
  double D1;		/** thickness of 1st layer (m) **/
  double D2;		/** thickness of 2nd layer (m) **/
  double dp;		/** depth to constant temperature (m) */
  double delta_t;	/** Time Step in Seconds **/
  double Le;		/** Latent heat of vapoization **/
  double Ls;		/** Latent heat of sublimation **/
  double Vapor;		/** Total vapor flux from snow in m/s **/
  double moist;		/** layer moisture content in m/m **/
  double ice0;		/** layer ice content in m/m **/
  double max_moist;	/** layer maximum moisture content in m/m **/
  double bubble;	/** bubbling pressure in cm **/
  double             expt;
  double             surf_atten;
  double             wind;
  double             displacement;
  double             roughness;
  double             ref_height;
  double             elevation;
  double             b_infilt;
  double             max_infil;
  double             dt;
  double             vpd;
  double             rainfall;
  double             Wdew;
  double             snow_energy;
  double            *grnd_flux;
  double            *T1;
  double            *latent_heat;
  double            *sensible_heat;
  double            *deltaH;
  double            *store_error;
  double            *TMean;
  double            *rad;
  double            *depth;
  double            *Wcr;
  double            *Wpwp;
  layer_data_struct *layer;
  veg_var_struct    *veg_var;
  char               VEG;
  char               CALC_EVAP;
  int                veg_class;
  int                month;
  int                iveg;

  /* Initialize Variables */
  T2            = (double) va_arg(ap, double);
  Ts_old        = (double) va_arg(ap, double);
  T1_old        = (double) va_arg(ap, double);
  Tair          = (double) va_arg(ap, double);
  ra            = (double) va_arg(ap, double);
  atmos_density = (double) va_arg(ap, double);
  shortwave     = (double) va_arg(ap, double);
  longwave      = (double) va_arg(ap, double);
  albedo        = (double) va_arg(ap, double);
  emissivity    = (double) va_arg(ap, double);
  kappa1        = (double) va_arg(ap, double);
  kappa2        = (double) va_arg(ap, double);
  Cs1           = (double) va_arg(ap, double);
  Cs2           = (double) va_arg(ap, double);
  D1            = (double) va_arg(ap, double);
  D2            = (double) va_arg(ap, double);
  dp            = (double) va_arg(ap, double);
  delta_t       = (double) va_arg(ap, double);
  Le            = (double) va_arg(ap, double);
  Ls            = (double) va_arg(ap, double);
  Vapor         = (double) va_arg(ap, double);
  moist         = (double) va_arg(ap, double);
  ice0          = (double) va_arg(ap, double);
  max_moist     = (double) va_arg(ap, double);
  bubble        = (double) va_arg(ap, double);
  expt          = (double) va_arg(ap, double);
  surf_atten    = (double) va_arg(ap, double);
  wind          = (double) va_arg(ap, double);
  displacement  = (double) va_arg(ap, double);
  roughness     = (double) va_arg(ap, double);
  ref_height    = (double) va_arg(ap, double);
  elevation     = (double) va_arg(ap, double);
  b_infilt      = (double) va_arg(ap, double);
  max_infil     = (double) va_arg(ap, double);
  dt            = (double) va_arg(ap, double);
  vpd           = (double) va_arg(ap, double);
  rainfall      = (double) va_arg(ap, double);
  Wdew          = (double) va_arg(ap, double);
  snow_energy   = (double) va_arg(ap, double);
  grnd_flux     = (double *) va_arg(ap, double *);
  T1            = (double *) va_arg(ap, double *);
  latent_heat   = (double *) va_arg(ap, double *);
  sensible_heat = (double *) va_arg(ap, double *);
  deltaH        = (double *) va_arg(ap, double *);
  store_error   = (double *) va_arg(ap, double *);
  TMean         = (double *) va_arg(ap, double *);
  rad           = (double *) va_arg(ap, double *);
  depth         = (double *) va_arg(ap, double *);
  Wcr           = (double *) va_arg(ap, double *);
  Wpwp          = (double *) va_arg(ap, double *);
  layer         = (layer_data_struct *) va_arg(ap, layer_data_struct *);
  veg_var       = (veg_var_struct *) va_arg(ap, veg_var_struct *);
  VEG           = (char) va_arg(ap, char);
  CALC_EVAP     = (char) va_arg(ap, char);
  veg_class     = (int)  va_arg(ap, int);
  month         = (int)  va_arg(ap, int);
  iveg          = (int)  va_arg(ap, int);

  /* Print Variables */
  fprintf(stderr,"T2 = %lf\n",T2);
  fprintf(stderr,"Ts_old = %lf\n",Ts_old);
  fprintf(stderr,"T1_old = %lf\n",T1_old);
  fprintf(stderr,"Tair = %lf\n",Tair);
  fprintf(stderr,"ra = %lf\n",ra);
  fprintf(stderr,"atmos_density = %lf\n",atmos_density);
  fprintf(stderr,"shortwave = %lf\n",shortwave);
  fprintf(stderr,"longwave = %lf\n",longwave);
  fprintf(stderr,"albedo = %lf\n",albedo);
  fprintf(stderr,"emissivity = %lf\n",emissivity);
  fprintf(stderr,"kappa1 = %lf\n",kappa1);
  fprintf(stderr,"kappa2 = %lf\n",kappa2);
  fprintf(stderr,"Cs1 = %lf\n",Cs1);
  fprintf(stderr,"Cs2 = %lf\n",Cs2);
  fprintf(stderr,"D1 = %lf\n",D1);
  fprintf(stderr,"D2 = %lf\n",D2);
  fprintf(stderr,"dp = %lf\n",dp);
  fprintf(stderr,"delta_t = %lf\n",delta_t);
  fprintf(stderr,"Le = %lf\n",Le);
  fprintf(stderr,"Ls = %lf\n",Ls);
  fprintf(stderr,"Vapor = %lf\n",Vapor);
  fprintf(stderr,"moist = %lf\n",moist);
  fprintf(stderr,"ice0 = %lf\n",ice0);
  fprintf(stderr,"max_moist = %lf\n",max_moist);
  fprintf(stderr,"bubble = %lf\n",bubble);
  fprintf(stderr,"expt = %lf\n",expt);
  fprintf(stderr,"surf_atten = %lf\n",surf_atten);
  fprintf(stderr,"wind = %lf\n",wind);
  fprintf(stderr,"displacement = %lf\n",displacement);
  fprintf(stderr,"roughness = %lf\n",roughness);
  fprintf(stderr,"ref_height = %lf\n",ref_height);
  fprintf(stderr,"elevation = %lf\n",elevation);
  fprintf(stderr,"b_infilt = %lf\n",b_infilt);
  fprintf(stderr,"max_infil = %lf\n",max_infil);
  fprintf(stderr,"dt = %lf\n",dt);
  fprintf(stderr,"vpd = %lf\n",vpd);
  fprintf(stderr,"rainfall = %lf\n",rainfall);
  fprintf(stderr,"Wdew = %lf\n",Wdew);
  fprintf(stderr,"snow_energy = %lf\n",snow_energy);
  fprintf(stderr,"grnd_flux = %lf\n",grnd_flux[0]);
  fprintf(stderr,"T1 = %lf\n",T1[0]);
  fprintf(stderr,"latent_heat = %lf\n",latent_heat[0]);
  fprintf(stderr,"sensible_heat = %lf\n",sensible_heat[0]);
  fprintf(stderr,"deltaH = %lf\n",deltaH[0]);
  fprintf(stderr,"store_error = %lf\n",store_error[0]);
  fprintf(stderr,"TMean = %lf\n",TMean[0]);
  fprintf(stderr,"rad = %lf\n",rad[0]);
  fprintf(stderr,"depth = %lf\n",depth[0]);
  fprintf(stderr,"Wcr = %lf\n",Wcr[0]);
  fprintf(stderr,"Wpwp = %lf\n",Wpwp[0]);
  fprintf(stderr,"VEG = %c\n",VEG);
  fprintf(stderr,"CALC_EVAP = %c\n",CALC_EVAP);
  fprintf(stderr,"veg_class = %i\n",veg_class);
  fprintf(stderr,"veg_class = %i\n",month);
  write_layer(layer,iveg,options.Nlayer,depth);
  write_vegvar(veg_var[0],iveg);

  vicerror("Finished writing calc_surf_energy_bal variables");

}
