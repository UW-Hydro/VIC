// $Id$
#if LAKE_MODEL

//#ifndef LAKE_SET

#define LAKE_SET
#define TMELT 0.0
#define EMICE 0.97      /* Ice emissivity */
#define EMH2O 0.98
#define RHOSNOW   250.  /* densities of water and snow */
#define RHOICE   917.   /* ice density*/
#define rhosurf 1.275   /* surface air density */
#define JOULESPCAL 4.1868
#define GRAMSPKG 1000.
#define SURF   0.09     /* surface thickness for E-B computation (m) */
#define FRACMIN  0.10   /* min ice thickness in meters */
#define FRACLIM   0.02  /* lower limit on fractional ice cover */
#define CPW_ICE   4200. /* specific heat of ice */
#define DM 1.38889E-07  /* molecular diffusivity of water */
#define SNOWCRIT   0.05  /* for albedo, in m */
//#define G 9.80616
#define ZWATER 0.0045    // 0.004 - original value
#define ZSNOW 0.005
#define CONDI 2.3        /* thermal conductivity of ice */
#define CONDS 0.31       /* thermal conductivity of snow */ 

// attenuation of short and longwave radiation through ice (W/m/K)
#define lamisw 4.5 // 3.5  // should be 1.5
#define lamilw 20   // should be 20.0
// attenuation of short and longwave radiation through snow (W/m/K)
#define lamssw 8 // 12.0 // 1.6  // should be 6.0 
#define lamslw 20 // 1.5  // should be 20.0
// attenuation of short and longwave radiation through water (W/m/K)
#define lamwsw 1.1  // 2.3 1.2
#define lamwlw .15 // .005  .05
#define  a1 0.7        /* Percent of radiation in visible band. */
#define  a2 0.3        /* Percent of radiation in infrared band. */
#define QWTAU 86400./2.   /* D. Pollard sub-ice time constant. */
#define RADIUS 6371.228 /* Earth radius in km. */

//#endif // LAKE_SET

/*** Subroutine prototypes ***/

double adjflux(double, double, double ,double, double, double, double,
	       double, double, double, double *, double *);
void alblake(double, double, float *, float *, float *, double, double, 
	     int, int *, double, char *);
void alloc_atmos(int, atmos_data_struct **);
double calc_density(double);
void colavg (double *, double *, double *, float, double *, int, double);
float dragcoeff(float, double, double);
void eddy (int, double, double * , double *, double *, double, int, double);
void energycalc(double *, double *, int, double, double *, double *, double *);
double get_dist(double, double, double, double);
void iceform (double *,double *,double ,double,double *,int, int, double, double *, double *, double *, double *);
void icerad(double,double ,double,double *, double *,double *);
void initialize_lake(lake_var_struct *, lake_con_struct, snow_data_struct *, 
		     double);
void lakeice(double *, double, double, double *, double, double *,int, 
	     double, double, double *, double, double, int, dmy_struct, double *);
void lakemain(atmos_data_struct *, lake_con_struct, double, double,
	      soil_con_struct *, int, dist_prcp_struct *, 
	      int, int, double, global_param_struct *, dmy_struct *, int, int);
void latsens(double,double, double, double, double, double, double, double,
	     double *, double *, double);
float lkdrag(float, double, double, double, double);
lake_con_struct read_lakeparam(FILE *, soil_con_struct, float, double *);
void rhoinit(double *, double);
void solve_lake(double, double, double, double, double, double, double, double, 
		double, double, double, lake_var_struct *, lake_con_struct, 
		soil_con_struct, int, int, energy_bal_struct *, 
		snow_data_struct *, double, dmy_struct);
double specheat (double);
void temp_area(double, double, double, double *, double *, double *, double *, int, double *, int, double, double*, double *, double *);
void tracer_mixer(double *, int *, int, double*, int, double, double *);
void tridia(int, double *, double *, double *, double *, double *);
void water_balance (lake_var_struct *, lake_con_struct, int, dist_prcp_struct *, int, int, double, soil_con_struct, double, double, double, double);
double func_lake_energy_balance(double, va_list);
double solve_surf_energy_bal(double Tsurf, ...);
void ice_melt(double, double, double *, double, snow_data_struct *, lake_var_struct *, int, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double *, double *, double *, double *, double *, double *, double *, double *, double *, double);
double CalcIcePackEnergyBalance(double Tsurf, ...);
double func_lake_energy_bal(double , va_list);
double IceEnergyBalance(double, va_list);
double ErrorPrintIcePackEnergyBalance(double, va_list);
double ErrorIcePackEnergyBalance(double Tsurf, ...);
void water_energy_balance(int, double*, double*, int,int,double, double ,double ,double ,double ,double ,double ,double ,double,double,double ,double ,double *,double *,double *,double*,double *, double *,double *,double ,double *,double *, double *, double *);
void water_under_ice(int, double,  double, double *, double *, double, int, double, double, double *, double *, double *, double *, int, double, double, double, double, double *);
void write_lake_var(lake_var_struct, int);
double get_sarea(lake_con_struct, double);
void wetland_energy(int, atmos_data_struct *, dist_prcp_struct *, dmy_struct *,
		    global_param_struct *, soil_con_struct  *, int,
		    int, double, lake_con_struct);
void update_prcp(dist_prcp_struct *, energy_bal_struct *,  snow_data_struct *, double, int, int, double, soil_con_struct); 
void initialize_prcp(dist_prcp_struct *, energy_bal_struct *,  snow_data_struct *, double, int, int, double, soil_con_struct, lake_var_struct *, int, lake_con_struct); 

#endif // LAKE_MODEL
