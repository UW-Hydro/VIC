/******************************************************************************
// $Id$
  Modifications:
  2006-Nov-07 Removed LAKE_MODEL option.				TJB
  2007-Aug-16 Made return value of initialize_prcp an int.		JCA
  2007-Aug-21 Added features for EXCESS_ICE option.			JCA
  2007-Oct-24 Changed get_sarea, get_volume, and get_depth to return exit
	      status so that errors can be trapped and communicated up the
	      chain of function calls.					KAC via TJB
  2007-Oct-24 Changed lakeice() to return exit status.			KAC via TJB
  2007-Nov-06 New lake physics parameters.  Modified argument lists for
	      various functions.  Moved get_dist() to vicNl.h.		LCB via TJB
  2008-Apr-21 Added argument to alblake.				LCB via TJB
  2008-Sep-10 Updated values of CONDS and lamwlw to match Laura
	      Bowling's lake work.					LCB via TJB
  2009-Jul-31 Removed lakemain() and wetland_energy(); initialize_lake
	      no longer takes a snow structure as input.		TJB
  2009-Sep-28 Removed initialize_prcp and update_prcp.  Modified
	      argument list of initialize_lake.				TJB
  2009-Sep-30 Miscellaneous fixes for lake model.			TJB
  2009-Oct-05 Added functions for updating/rescaling lake and wetland
	      fluxes and storages when lake area changes.		TJB
  2010-Nov-11 Added skip_hydro flag to initialize_lake() arg list.
	      Removed rescale_lake_fluxes().				TJB
  2010-Nov-26 Changed the argument list of water_balance().		TJB
  2010-Dec-28 Added latitude to alblake() arglist.			TJB
  2011-Mar-01 Added rescale_snow_storage().  Added terms to argument
	      list of initialize_lake().				TJB
******************************************************************************/

//#ifndef LAKE_SET

#define LAKE_SET
#define TMELT 0.0
#define EMICE 0.97      /* Ice emissivity */
#define EMH2O .98
#define RHOSNOW   250.  /* densities of water and snow */
#define RHOICE   917.   /* ice density*/
#define rhosurf 1.275   /* surface air density */
#define MAX_SURFACE_LAKE   .6  /* max. surface layer thickness for E-B (m) */
#define BETA 0.001       /* Curve shape parameter for lake profile. */
#define FRACMIN  0.10   /* min ice thickness in meters */
#define FRACLIM   0.02  /* lower limit on fractional ice cover */
#define CPW_ICE   4200. /* specific heat of ice */
#define DM 1.38889E-07  /* molecular diffusivity of water */
#define SNOWCRIT   0.05  /* for albedo, in m */
//#define G 9.80616
#define ZWATER 0.0045    // 0.004 - original value
#define ZSNOW 0.005
#define CONDI 2.3        /* thermal conductivity of ice */
#define CONDS 0.7       /* thermal conductivity of snow */ 

// attenuation of short and longwave radiation through ice (1/m)
#define lamisw 1.5 // 1.5 in Patterson & Hamblin
#define lamilw 20  // 20.0 in Patterson & Hamblin
// attenuation of short and longwave radiation through snow (1/m)
#define lamssw 6.0 // 6.0 in Patterson & Hamblin
#define lamslw 20  // 20.0 in Patterson & Hamblin
// attenuation of short and longwave radiation through water (1/m)
#define lamwsw .3  // San Fran Bay data: 0.31 - 29.9 1/m (visible)
#define lamwlw 1.4 // Hostetler and Bartlein assume 0.85 1/m (total)
#define  a1 0.7        /* Percent of radiation in visible band. */
#define  a2 0.3        /* Percent of radiation in infrared band. */
#define QWTAU 86400./2.   /* D. Pollard sub-ice time constant. */
#define RADIUS 6371.228 /* Earth radius in km. */

//#endif // LAKE_SET

/*** Subroutine prototypes ***/

void advect_soil_veg_storage(double, double, double, double *, soil_con_struct *, veg_con_struct *, cell_data_struct *, veg_var_struct *, lake_con_struct);
void advect_snow_storage(double, double, double, snow_data_struct *);
void alblake(double, double, double *, double *, float *, float *, double, double, 
	     int, int *, double, double, char *, int, double);
double calc_density(double);
double CalcIcePackEnergyBalance(double Tsurf, ...);
void colavg (double *, double *, double *, float, double *, int, double, double);
float dragcoeff(float, double, double);
void eddy (int, double, double * , double *, double *, double, int, double, double);
void energycalc(double *, double *, int, double, double,double *, double *, double *);
double ErrorIcePackEnergyBalance(double Tsurf, ...);
double ErrorPrintIcePackEnergyBalance(double, va_list);
int get_depth(lake_con_struct, double, double *);
int get_sarea(lake_con_struct, double, double *);
int get_volume(lake_con_struct, double, double *);
void iceform (double *,double *,double ,double,double *,int, int, double, double, double *, double *, double *, double *, double *, double);
void icerad(double,double ,double,double *, double *,double *);
int ice_melt(double, double, double *, double, snow_data_struct *, lake_var_struct *, int, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double *, double *, double *, double *, double *, double *, double *, double *, double *, double);
double IceEnergyBalance(double, va_list);
int initialize_lake(lake_var_struct *, lake_con_struct, soil_con_struct *, cell_data_struct *, double, int);
int lakeice(double *, double, double, double, double, int, 
	    double, double, double *, double, double, int, dmy_struct, double *, double *, double, double);
void latsens(double,double, double, double, double, double, double, double,
	     double *, double *, double);
float lkdrag(float, double, double, double, double);
lake_con_struct read_lakeparam(FILE *, soil_con_struct, veg_con_struct *);
void rescale_soil_veg_fluxes(double, double, cell_data_struct *, veg_var_struct *);
void rescale_snow_energy_fluxes(double, double, snow_data_struct *, energy_bal_struct *);
void rescale_snow_storage(double, double, snow_data_struct *);
void rhoinit(double *, double);
int solve_lake(double, double, double, double, double, double, double, double, 
		double, double, lake_var_struct *, lake_con_struct, 
		soil_con_struct, int, int, double, dmy_struct, double);
double specheat (double);
void temp_area(double, double, double, double *, double *, double *, double *, int, double *, int, double, double, double*, double *, double *);
void tracer_mixer(double *, int *, int, double*, int, double, double, double *);
void tridia(int, double *, double *, double *, double *, double *);
int water_balance (lake_var_struct *, lake_con_struct, int, dist_prcp_struct *, int, int, int, double, soil_con_struct,
#if EXCESS_ICE
		    veg_con_struct, int, double);
#else
		    veg_con_struct);
#endif		    
int  water_energy_balance(int, double*, double*, int, int, double, double, double, double, double, double, double, double, double, double, double, double, double, double *, double *, double *, double*, double *, double *, double *, double, double *, double *, double *, double *, double *, double);
int water_under_ice(int, double,  double, double *, double *, double, int, double, double, double, double *, double *, double *, double *, int, double, double, double, double *);
