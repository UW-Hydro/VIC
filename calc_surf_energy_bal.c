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
