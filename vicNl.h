#include <math.h>
#include <vicNl_def.h>

// $Id$

#if LAKE_MODEL
#include <LAKE.h>
#endif // LAKE_MODEL

/*** SubRoutine Prototypes ***/

double advected_sensible_heat(double, double, double, double, double);
void alloc_atmos(int, atmos_data_struct **);
double arno_evap(layer_data_struct *, layer_data_struct *, double, double, 
		 double, double, double, double, double, double, double, 
		 double, double, double, double, double, 
#if SPATIAL_FROST
		 double, double *);
#else
		 double);
#endif // SPATIAL_FROST
unsigned char average_moisture_for_storm(double *, double *, double, double);

void   CalcAerodynamic(char, double, double, double, double, double, double, 
	  	       double *, double *, double *, double *, double *, 
                       int, int);
void   calc_cloud_cover_fraction(atmos_data_struct *, dmy_struct *, int,
				 int, int, double *);
void   calc_energy_balance_error(int, double, double, double, double, double);
#if OUTPUT_FORCE_STATS
void   calc_forcing_stats(int, atmos_data_struct *);
#endif // OUTPUT_FORCE_STATS
void   calc_longwave(double *, double, double, double);
void   calc_netlongwave(double *, double, double, double);
double calc_netshort(double, int, double, double *);
double calc_rainonly(double,double,double,double,double);
void   calc_root_fractions(veg_con_struct *, soil_con_struct *);
double calc_snow_coverage(int *, double, double, double, double, double, 
                          double, double, double *, double *, double *, 
                          double *, double *);
double calc_snow_ground_flux(int, int, int, int, double, double, double, 
			     double, double, double *, double *, double *, 
			     double *, energy_bal_struct *, 
			     snow_data_struct *, layer_data_struct *,
                             layer_data_struct *, soil_con_struct *, char *);
#if QUICK_FS
int    calc_soil_thermal_fluxes(int, double *, double *, double *, double *, 
				double *, double *, double *, double *, 
				double *, double *, double *, double *, 
				double *, double *, double ***, int, int);
#else
int    calc_soil_thermal_fluxes(int, double *, double *, double *, double *, 
				double *, double *, double *, double *, 
				double *, double *, double *, double *, 
				double *, double *, int, int);
#endif // QUICK_FS
double CalcSnowPackEnergyBalance(double Tsurf, ...);
double CalcBlowingSnow(double, double, int, double, double, double, double, 
                       double, double, double *, double, double, float, 
                       float, double, int, int, float, double, double, double *); 
double calc_atmos_energy_bal(double, double, double, double, double, double, 
                             double, double, double, double, double, double, 
                             double, double, double, double, 
                             double *, double *, double *, double *, 
                             double *, double *, double *, double *);
double calc_surf_energy_bal(double, double, double, double, double, double, 
                            double, double, double, double, double, double, 
                            double, double, double, double, double, double, 
                            double, double, double, double, double, double, 
                            double, double, double, 
                            double *, double *, double *, double *, double *, 
                            double *, double *, double *, double *, double *, 
                            float *, int, int, 
                            int, int, int, int, int, int, int, int, int, int, 
                            atmos_data_struct *, dmy_struct *, 
                            energy_bal_struct *, layer_data_struct *, 
                            layer_data_struct *, 
                            snow_data_struct *, soil_con_struct *, 
                            veg_var_struct *, veg_var_struct *);
double calc_trans(double, double);
double calc_veg_displacement(double);
double calc_veg_height(double);
double calc_veg_roughness(double);
void   calc_water_balance_error(int, double, double, double);
double canopy_evap(layer_data_struct *, layer_data_struct *,
		   veg_var_struct *, veg_var_struct *, char, int, int, 
		   double, double *, double, double, double, double, 
		   double, double, double, double, double, double, 
		   double *, double *, double *, double *, 
#if SPATIAL_FROST
                   double *, float *);
#else
                   float *);
#endif
void   check_files(infiles_struct *, filenames_struct *);
FILE  *check_state_file(char *, dmy_struct *, global_param_struct *, int, int, 
                        int *);
void   close_files(infiles_struct *, outfiles_struct *, filenames_struct *);
filenames_struct cmd_proc(int argc, char *argv[]);
void   compress_files(char string[]);
void   compute_dz(double *, double *, int, double);
void   compute_penman_constants(double, double, double, double, double, 
				double, double, float, float, double *, 
				double *, double *, double *, double *);
void   correct_precip(double *, double, double, double, double);
void   compute_soil_layer_thermal_properties(layer_data_struct *, double *,
					     double *, double *, double *, 
#if SPATIAL_FROST
                                             double *,
#endif
					     int);
void   compute_treeline(atmos_data_struct *, dmy_struct *, double *, char *);
void display_current_settings(int, filenames_struct *, global_param_struct *);
void   dist_prec(atmos_data_struct *,dist_prcp_struct *,soil_con_struct *,
		 veg_con_struct *,dmy_struct *,global_param_struct *,
#if LAKE_MODEL
		 lake_con_struct *,
#endif /* LAKE_MODEL */
		 outfiles_struct *, int, int, char, char, char *, int *);
#if QUICK_FS
void distribute_node_moisture_properties(double *, double *, double *,
					 double *, double *, double *,
					 double *, double ***, 
					 double *, double *, double *,
					 double *, double *, int, int, char);
#else
void distribute_node_moisture_properties(double *, double *, double *,
					 double *, double *, double *,
					 double *, double *, double *,
					 double *, double *, double *,
					 double *, double *, int, int, char);
#endif
void   distribute_soil_property(double *,double,double,
				double **l_param,
				int, int, double *, double *);

double error_calc_atmos_energy_bal(double Tcanopy, ...);
double error_calc_atmos_moist_bal(double , ...);
double error_calc_canopy_energy_bal(double Tsurf, ...);
double error_calc_snow_ground_flux(double Tsurf, ...);
double error_calc_surf_energy_bal(double Tsurf, ...);
double ErrorSnowPackEnergyBalance(double Tsurf, ...);
double error_print_atmos_energy_bal(double, va_list);
double error_print_atmos_moist_bal(double, va_list);
double error_print_canopy_energy_bal(double, va_list);
double error_print_snow_ground_flux(double, va_list);
double ErrorPrintSnowPackEnergyBalance(double, va_list);
double error_print_solve_T_profile(double, va_list);
double error_print_surf_energy_bal(double, va_list);
double error_solve_T_profile(double Tsurf, ...);
double estimate_dew_point(double, double, double, double, double);
#if QUICK_FS
void estimate_layer_ice_content(layer_data_struct *, double *, double *,
				double *, double ***, double *,
				double *, double ***, 
#if SPATIAL_FROST
                                double *, double,
#endif // SPATIAL_FROST
                                double *, double *, double *, float **, 
                                int, int, char);
#else
void estimate_layer_ice_content(layer_data_struct *, double *, double *,
				double *, double *, double *, double *,
				double *, double *, double *, 
#if SPATIAL_FROST
                                double *, double, 
#endif // SPATIAL_FROST
                                double *, double *, double *, float **, 
                                int, int, char);
#endif
double estimate_T1(double, double, double, double, double, double, double, 
		   double, double, double, double);
double exp_interp(double,double,double,double,double);

double f(double, double, double, double, double, double, double, double,
         double, double, int, double *, double, double, double, double *,
         double *, double *, double *, double *, double *);
void   find_0_degree_fronts(energy_bal_struct *,double *, double *, int);
layer_data_struct find_average_layer(layer_data_struct *, layer_data_struct *,
				     double, double);
void   find_sublayer_temperatures(layer_data_struct *, double *, double *,
				  double *, double, double, int, int);
void   finish_frozen_soil_calcs(energy_bal_struct *, layer_data_struct *,
				layer_data_struct *, layer_data_struct *,
				soil_con_struct *, int, int, double, 
				double *, double *, double *, double *);
void   free_atmos(int nrecs, atmos_data_struct **atmos);
void   free_dist_prcp(dist_prcp_struct *, int);
void   free_vegcon(veg_con_struct **);
void   full_energy(char, int, int, atmos_data_struct *, dist_prcp_struct *,
		   dmy_struct *, global_param_struct *, 
#if LAKE_MODEL
		   lake_con_struct *,
#endif // LAKE_MODEL 
                   soil_con_struct *, veg_con_struct *);
double func_aero_resist(double,double,double,double,double);
double func_atmos_energy_bal(double, va_list);
double func_atmos_moist_bal(double, va_list);
double func_canopy_energy_bal(double, va_list);
double func_snow_ground_flux(double, va_list);
double func_surf_energy_bal(double, va_list);

double get_avg_temp(double, double, double *, double *, int);
void   get_force_type(char *, int, int *);
global_param_struct get_global_param(filenames_struct *, FILE *);
void   get_next_time_step(int *, int *, int *, int *, int *, int);

double hermint(double, int, double *, double *, double *, double *, double *);
void   hermite(int, double *, double *, double *, double *, double *);
void   HourlyT(int, int, int *, double *, int *, double *, double *);

void   initialize_atmos(atmos_data_struct *, dmy_struct *, FILE **, double, 
			double, double, double, double, double, double, 
                        double *, 
#if OUTPUT_FORCE
			char *, outfiles_struct *);
#else
			char *);
#endif
void   initialize_global();
void   initialize_model_state(double, int, int, int, int, dist_prcp_struct *, 
                              dmy_struct, global_param_struct *, 
                              infiles_struct, 
#if LAKE_MODEL
                              lake_con_struct,
#endif // LAKE_MODEL
                              soil_con_struct *, veg_con_struct *, char **, 
                              int **);
void   initialize_new_storm(cell_data_struct ***, veg_var_struct ***,
			    int, int, int, double, double);
void   initialize_snow(snow_data_struct **,int,FILE *,int);
void   initialize_soil(cell_data_struct **, soil_con_struct *, int);
void   initialize_veg( veg_var_struct **, veg_con_struct *,
		       global_param_struct *, int);
void   latent_heat_from_snow(double, double, double, double, double, 
                             double, double, double *, double *, 
                             double *, double *, double *);
double linear_interp(double,double,double,double,double);

cell_data_struct **make_cell_data(int, int);
dist_prcp_struct make_dist_prcp(int, int *);
dmy_struct *make_dmy(global_param_struct *);
energy_bal_struct **make_energy_bal(int, int *);
filenames_struct make_in_and_outfiles(infiles_struct *, filenames_struct *, 
				      soil_con_struct *, outfiles_struct *);
out_data_struct *make_out_data(int);
snow_data_struct **make_snow_data(int);
veg_var_struct **make_veg_var(int);
void   MassRelease(double *,double *,double *,double *);
double maximum_unfrozen_water(double, double, double, double);
#if QUICK_FS
double maximum_unfrozen_water_quick(double, double, double **);
#endif
double modify_Ksat(double);
void mtclim42_wrapper(int, int, double, double, double, double, 
		      global_param_struct *, dmy_struct *, double *, 
		      double *, double *, double *, double *, double *);

double new_snow_density(double);
void   nrerror(char *);

void   open_debug();
FILE  *open_file(char string[], char type[]);
#if SAVE_STATE
FILE  *open_state_file(global_param_struct *, int, int);
#endif

double penman(double, double, double, double, double, double, double, 
	      double, double, float, float);
void   prepare_full_energy(int, int, int, dist_prcp_struct *, 
			   soil_con_struct *, double *, double *);
double priestley(double, double);
void   put_data(char *, double *, double *, double *, 
#if SPATIAL_FROST
                double *, double,
#endif // SPATIAL_FROST
                double, int, int, int, int, atmos_data_struct *, 
                dist_prcp_struct *, dmy_struct *,
#if LAKE_MODEL
                lake_con_struct *, 
#endif // LAKE_MODEL
 	        outfiles_struct *, veg_con_struct *); 

double quick_penman(double, double, double, double, double, double, 
		    double, double);
double read_arcinfo_value(char *, double, double);
int    read_arcinfo_info(char *, double **, double **, int **);
void   read_atmos_data(FILE *, global_param_struct, int, int, double **);
double **read_forcing_data(FILE **, global_param_struct);
void   read_initial_model_state(FILE *, dist_prcp_struct *, 
				global_param_struct *, int, int, int, 
				soil_con_struct *, int, char *, int *);
void   read_PILPS2c(atmos_data_struct *, FILE *, int *, int, int, int);
void   read_rosemount(atmos_data_struct *, FILE *, int *, int, int, int, int);
void   read_sawd(atmos_data_struct *, FILE *, int *, int, int, int);
void   read_sawd_binary(atmos_data_struct *, FILE *, int *, int, int, int);
void   read_snowband(FILE *, int, double, double **, double **, double **, 
                     char **);
void   read_snowmodel(atmos_data_struct *, FILE *, int, int, int, int);
soil_con_struct read_soilparam(FILE *, int);
soil_con_struct read_soilparam_arc(FILE *, char *, int *, int *, int);
veg_lib_struct *read_veglib(FILE *, int *);
veg_con_struct *read_vegparam(FILE *, int, int);
void   redistribute_during_storm(cell_data_struct ***, veg_var_struct ***,
				 int, int, int, double, double, double, 
				 double *);
void   redistribute_moisture(layer_data_struct *, double *, double *,
			     double *, double *, double *, int);
unsigned char redistribute_moisture_for_storm(double *, double *, double, 
					      double, double);
double root_brent(double, double, char *, double (*Function)(double, va_list), ...);
void   runoff(layer_data_struct *, layer_data_struct *, energy_bal_struct *, 
	      soil_con_struct *, double *, double *, double *, double *, 
	      double *, 
#if SPATIAL_FROST
              double *, 
#endif
              double, int, int, int, int, int);

void set_max_min_hour(double *, int, int *, int *);
void set_node_parameters(double *, double *, double *, double *, double *,
			 double *, double *, double *, double *, double *,
			 double *, double *, float **,
#if QUICK_FS
			 double ***,
#endif
			 int, int, char);
void   setup_frozen_soil(soil_con_struct *, layer_data_struct  *,
			 layer_data_struct *, layer_data_struct *,
			 energy_bal_struct, int, int, int, double,
			 double *, double *, double *);
double shrad(double,double,double,double,double,int,double);
double snow_albedo(double, double, double, double, int, char);
double snow_density(int, double, double, double, double, double, double, 
		    double);
void   snow_intercept(double, double, double, double, double, double, double,
                      double, double, double, double, double, double, double, 
                      double, double, 
                      double *, double *, double *, double *, double *, 
                      double *, double *, double *, double *, double *, 
                      double *, double *, double *, double *, double *, 
                      double *, double *, double *, double *, double *, 
                      double *, double *, double *, float *, int, int, int, 
                      int, int, int, int, layer_data_struct *, 
                      layer_data_struct *, soil_con_struct *, 
                      veg_var_struct *, veg_var_struct *);
void   snow_melt(double, double, double, double, double *, double, double *, double, 
		 double, double, double, double, double, double, double, 
                 double, double, double, double, double, double, 
                 double *, double *, double *, double *, double *, double *, 
                 double *, double *, double *, double *, double *, double *, 
                 int, int, int, snow_data_struct *, soil_con_struct *);
double SnowPackEnergyBalance(double, va_list);
double soil_conductivity(double, double, double, double, double);
void   soil_thermal_calc(soil_con_struct *, layer_data_struct *,
			 energy_bal_struct, double *, double *, double *,
			 int, int);
double soil_thermal_eqn(double, va_list);
double solve_snow(char, double, double, double, double, double, double, 
                  double, double, double, double, double, double, double, 
                  double, double, double, double, double, double, double, 
                  double *, double *,
                  double *, double *, double *, double *, double *, double *, 
                  double *, double *, double *, double *, double *, double *, 
                  double *, double *, double *, double *, double *, double *, 
                  double *, double *, double *, double *, double *, double *,
                  float *, 
                  int, int, int, int, int, int, int, int, int, 
                  int, int, int, int, int *, energy_bal_struct *, 
                  layer_data_struct *, layer_data_struct *, 
                  snow_data_struct *, soil_con_struct *, 
                  veg_var_struct *, veg_var_struct *);
double solve_atmos_energy_bal(double Tcanopy, ...);
double solve_atmos_moist_bal(double , ...);
double solve_canopy_energy_bal(double Tfoliage, ...);
double solve_snow_ground_flux(double Tsurf, ...);
double solve_surf_energy_bal(double Tsurf, ...);
#if QUICK_FS
void   solve_T_profile(double *, double *, double *, double *,double *,
		       double *, double, double *, double *, double *,
		       double *, double *, double *, double *, double ***, 
		       int, int *, int, int, int);
#else
void   solve_T_profile(double *, double *, double *, double *,double *,
		       double *, double, double *, double *, double *,
		       double *, double *, double *, double *, int, int *, 
		       int, int, int);
#endif
double StabilityCorrection(double, double, double, double, double, double);
void   store_moisture_for_debug(int,int,double *,cell_data_struct ***,
				veg_var_struct ***,snow_data_struct **,
				soil_con_struct *);
void   surface_fluxes(char, double, double, double, double, double, double, 
                      double *, double *, double *, double *, double *, 
                      double *, double *, double *, double *, double *, 
                      double *, double *, double *, double *, double *, 
                      double *, double *, float *, int, int, int, int, int, 
                      int, int, int, int, atmos_data_struct *, dmy_struct *, 
                      energy_bal_struct *, global_param_struct *, 
                      layer_data_struct *, layer_data_struct *, 
                      snow_data_struct *, soil_con_struct *, 
                      veg_var_struct *, veg_var_struct *, float, float, float);
double svp(double);
double svp_slope(double);

void transpiration(layer_data_struct *, int, int, double, double, double, 
		   double, double, double, double, double, double, double, 
		   double *, double *, double *, double *, double *, double *,
#if SPATIAL_FROST
                   double *,
#endif
                   float *);
void tridag(double *,double *,double *,double *,double *,int);

void usage(char *);

void   vicerror(char *);
double volumetric_heat_capacity(double,double,double);

void write_atmosdata(atmos_data_struct *, int);
void write_data(out_data_struct *, outfiles_struct *, dmy_struct *, int);
void write_debug(atmos_data_struct *, soil_con_struct *, cell_data_struct *,
                 energy_bal_struct *, snow_data_struct *, veg_var_struct *,
                 dmy_struct *, global_param_struct *,
                 double, double, int, int, int, int, int, char);
void write_dist_prcp(dist_prcp_struct *);
#if OUTPUT_FORCE
void write_forcing_file(atmos_data_struct *, int, outfiles_struct *);
#endif
void write_layer(layer_data_struct *, int, int, 
#if SPATIAL_FROST
                 double *,
#endif
                 double *);
#if SAVE_STATE
void write_model_state(dist_prcp_struct *, global_param_struct *, int, 
		       int, outfiles_struct *, soil_con_struct *, char *, int *);
#endif
void write_snow_data(snow_data_struct, int, int);
void write_soilparam(soil_con_struct *);
void write_vegparam(veg_con_struct *);
void write_vegvar(veg_var_struct *, int);

