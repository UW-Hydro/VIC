/* header file for vic_driver_classic routines */
#include <math.h>
#include <vicNl_def.h>

/* ordered alphabetically (case-insensitive) by function name */
// lake: advect_carbon_storage
// lake: advect_snow_storage
// lake: advect_soil_veg_storage
// vic_run: double advected_sensible_heat(double, double, double, double, double);
void alloc_atmos(int, atmos_data_struct **);
// vic_run: double arno_evap(layer_data_struct *, layer_data_struct *, double, double, 
//		 double, double, double, double, double, double, double, double, 
//		 double, double *);
// mtclim_constants_vic.h: atm_pres
unsigned char average_moisture_for_storm(double *, double *, double, double);
// vic_run: double calc_atmos_energy_bal(double, double, double, double, double, double, 
//                             double, double, double, double, double, double, 
//                             double, double, double, double, 
//                             double *, double *, double *, double *, 
//                             double *, double *, double *, double *, char *, int *);
// lake: calc_density
// vic_run: double calc_energy_balance_error(int, double, double, double, double, double);
// vic_run: int    calc_layer_average_thermal_props(energy_bal_struct *, layer_data_struct *,
//          layer_data_struct *, layer_data_struct *,
//          soil_con_struct *, int, int, double *);
void   calc_longwave(double *, double, double, double);
void   calc_netlongwave(double *, double, double, double);
double calc_netshort(double, int, double, double *);
// vic_run: void calc_Nscale_factors(char, double *, double, double, double, double,
//                         dmy_struct, double *);
// mtclim_constants_vic.h: calc_pet
// mtclim_constants_vic.h: calc_prcp
double calc_rainonly(double,double,double,double,double);
// vic_run: double calc_rc(double,double,float,double,double,double,double,char);
// vic_run: void calc_rc_ps(char, double, double, double, double *, double,
//                double, double *, double, double, double *,
//                double, double, double, double *, double *);
void   calc_root_fractions(veg_con_struct *, soil_con_struct *);
// vic_run: double calc_snow_coverage(int *, double, double, double, double, double, 
//                          double, double, double *, double *, double *, 
//                          double *, double *);
// remove: double calc_snow_ground_flux(int, int, int, int, double, double, double, 
//           double, double, double *, double *, double *, 
//           double *, energy_bal_struct *, 
//           snow_data_struct *, layer_data_struct *,
//                             layer_data_struct *, soil_con_struct *, char *);
// vic_run: int    calc_soil_thermal_fluxes(int, double *, double *, char *, int *, double *, double *, 
//        double *, double *, double *,double *, 
//        double *, double *, double *, 
//        double *, double *, double *, 
//        int, int, int, int);
// mtclim_constants_vic.h: calc_srad_humidity_iterative
// vic_run: double calc_surf_energy_bal(double, double, double, double, double, double,
//                            double, double, double, double, double, double,
//                            double, double, double, double, double, double,
//                            double, double, double, double, double, double,
//                            double, double, double,
//                            double *, double *, double *, double *, double *,
//                            double *, double *, double *, double *, double *,
//                            float *, int, int,
//                            int, int, int, int, int, int, int, int, int, int,
//                            double *, double *,
//                            atmos_data_struct *, dmy_struct *,
//                            energy_bal_struct *, layer_data_struct *,
//                            layer_data_struct *,
//                            snow_data_struct *, soil_con_struct *,
//                            veg_var_struct *, veg_var_struct *, int);
// mtclim_constants_vic.h: calc_tair
// remove: double calc_trans(double, double);
double calc_veg_height(double);
// vic_run: double calc_water_balance_error(int, double, double, double);
// vic_run: int    CalcAerodynamic(char, double, double, double, double, double,
//	  	       double *, double *, double *, double *, double *);
// vic_run: double CalcBlowingSnow(double, double, int, double, double, double, double, 
//                       double, double, double, double, double, float, 
//                       float, double, int, int, float, double, double, double *); 
// lake: CalcIcePackEnergyBalance
// vic_run: double CalcSnowPackEnergyBalance(double Tsurf, ...);
// CalcBlowingSnow: CalcSubFlux
// vic_run: void canopy_assimilation(char, double, double, double, double *, double,
//                         double, double *, double, double, double *,
//                         double, char *, double *, double *,
//                         double *, double *, double *, double *,
//                         double *, double *, double *, double *);
// vic_run: double canopy_evap(layer_data_struct *, layer_data_struct *,
//		   veg_var_struct *, veg_var_struct *, char, int, int, 
//		   double, double *, double, double, double, double, 
//		   double, double, double, double, double, double, 
//		   double *, double *, double *, double *, double *, 
//                   double *, float *, double *, double, double, double *);
void   check_files(filep_struct *, filenames_struct *);
FILE  *check_state_file(char *, dmy_struct *, global_param_struct *, int, int, 
                        int *);
void   close_files(filep_struct *, out_data_file_struct *, filenames_struct *);
filenames_struct cmd_proc(int argc, char *argv[]);
// vic_run: void   collect_eb_terms(energy_bal_struct, snow_data_struct, cell_data_struct,
//                        int *, int *, int *, int *, int *, double, double, double,
//                        int, int, double, int, int, double *, double *,
//                        double *, double, out_data_struct *);
// vic_run: void   collect_wb_terms(cell_data_struct, veg_var_struct, snow_data_struct, lake_var_struct,
//                        double, double, double, double, int, int, double, int, double *,
//                        double *, out_data_struct *);
void   compress_files(char string[]);
// vic_run: double compute_coszen(double, double, double, dmy_struct);
// vic_run: void   compute_pot_evap(int, dmy_struct *, int, int, double, double , double, double, double, double **, double *);
// vic_run: void   compute_runoff_and_asat(soil_con_struct *, double *, double, double *, double *);
// vic_run: void   compute_soil_resp(int, double *, double, double, double *, double *,
//                         double, double, double, double *, double *, double *);
// vic_run: void   compute_soil_layer_thermal_properties(layer_data_struct *, double *,
//                                             double *, double *, double *, 
//                                             double *, double *, double *, 
//                                             double *, int);
// mtclim_constants_vic.h: compute_srad_humidity_onetime
void   compute_treeline(atmos_data_struct *, dmy_struct *, double, double *, char *);
// vic_run: double compute_zwt(soil_con_struct *, int, double);
void   correct_precip(double *, double, double, double, double);
out_data_struct *create_output_list();
// vic_run: double darkinhib(double);
// mtclim_constants_vic.h: data_alloc
// mtclim_constants_vic.h: data_free
void   display_current_settings(int, filenames_struct *, global_param_struct *);
// vic_run: int    dist_prec(atmos_data_struct *,dist_prcp_struct *,soil_con_struct *,
//		 veg_con_struct *, lake_con_struct *,
//		 dmy_struct *,global_param_struct *,
//		 filep_struct *, out_data_file_struct *,
//		 out_data_struct *, save_data_struct *,
//		 int, int, char, char, char *, int *);
// vic_run: int  distribute_node_moisture_properties(double *, double *, double *, 
//					 double *, double *, double *,
//					 double *, double *, double *,
//					 double *, double *, double *, double *, double *,
//					 double *, double *, double *, int, int, char);
// remove: void   distribute_soil_property(double *,double,double,
//				double **l_param,
//				int, int, double *, double *);
// vic_run: double error_calc_atmos_energy_bal(double Tcanopy, ...);
// vic_run: double error_calc_atmos_moist_bal(double , ...);
// vic_run: double error_calc_canopy_energy_bal(double Tsurf, ...);
// remove: double error_calc_snow_ground_flux(double Tsurf, ...);
// vic_run: double error_calc_surf_energy_bal(double Tsurf, ...);
// vic_run: double error_print_atmos_energy_bal(double, va_list);
// vic_run: double error_print_atmos_moist_bal(double, va_list);
// vic_run: double error_print_canopy_energy_bal(double, va_list);
// remove: double error_print_snow_ground_flux(double, va_list);
// vic_run: double error_print_solve_T_profile(double, va_list);
// vic_run: double error_print_surf_energy_bal(double, va_list);
// vic_run: double error_solve_T_profile(double Tsurf, ...);
// lake: ErrorIcePackEnergyBalance
// lake: ErrorPrintIcePackEnergyBalance
// vic_run: double ErrorPrintSnowPackEnergyBalance(double, va_list);
// vic_run: double ErrorSnowPackEnergyBalance(double Tsurf, ...);
// remove: double estimate_dew_point(double, double, double, double, double);
// vic_run: int estimate_layer_ice_content(layer_data_struct *, double *, double *,
//			       double *, double *, double *, double *,
//			       double *, double *, double *, 
//			       double *, double, int, int, char);
// vic_run: int estimate_layer_ice_content_quick_flux(layer_data_struct *, double *,
//					  double, double, double, double,
//					  double *, double *, double *,
//					  double *, double, char);
// vic_run: double estimate_T1(double, double, double, double, double, double, double, 
//		   double, double, double, double);
// vic_run: double exp_interp(double,double,double,double,double);
// remove: double f(double, double, double, double, double, double, double, double,
//         double, double, int, double *, double, double, double, double *,
//         double *, double *, double *, double *, double *);
// vic_run: void   faparl(double *, double, double, double, double, double *, double *);
// vic_run: void   fda_heat_eqn(double *, double *, int, int, ...);
// vic_run: void   fdjac3(double *, double *, double *, double *, double *,
//            void (*vecfunc)(double *, double *, int, int, ...), 
//            int);
// vic_run: void   find_0_degree_fronts(energy_bal_struct *, double *, double *, int);
// vic_run: layer_data_struct find_average_layer(layer_data_struct *, layer_data_struct *,
//				     double, double);
// remove: void   find_sublayer_temperatures(layer_data_struct *, double *, double *,
//				  double *, double, double, int, int);
void   free_atmos(int nrecs, atmos_data_struct **atmos);
void   free_dist_prcp(dist_prcp_struct *, int);
void   free_dmy(dmy_struct **dmy);
void   free_out_data_files(out_data_file_struct **);
void   free_out_data(out_data_struct **);
void   free_vegcon(veg_con_struct **);
void   free_veglib(veg_lib_struct **);
// vic_run: int    full_energy(char, int, int, atmos_data_struct *, dist_prcp_struct *,
//		   dmy_struct *, global_param_struct *, lake_con_struct *,
//                   soil_con_struct *, veg_con_struct *);
// remove: double func_aero_resist(double,double,double,double,double);
// vic_run: double func_atmos_energy_bal(double, va_list);
// vic_run: double func_atmos_moist_bal(double, va_list);
// vic_run: double func_canopy_energy_bal(double, va_list);
// remove: double func_snow_ground_flux(double, va_list);
// vic_run: double func_surf_energy_bal(double, va_list);
// remove: double get_avg_temp(double, double, double *, double *, int);
// lake: get_depth
// lake: get_depth_from_sarea
double get_dist(double, double, double, double);
void   get_force_type(char *, int, int *);
global_param_struct get_global_param(filenames_struct *, FILE *);
void   get_next_time_step(int *, int *, int *, int *, int *, int);
// CalcBlowingSnow: get_prob
// lake: get_sarea
// CalcBlowingSnow: get_shear
// CalcBlowingSnow: get_thresh
// lake: get_volume
double hermint(double, int, double *, double *, double *, double *, double *);
void   hermite(int, double *, double *, double *, double *, double *);
// vic_run: double hiTinhib(double);
void   HourlyT(int, int, int *, double *, int *, double *, double *);
// lake: ice_depth
// lake: ice_melt
// lake: IceEnergyBalance
void   init_output_list(out_data_struct *, int, char *, int, float);
void   initialize_atmos(atmos_data_struct *, dmy_struct *, FILE **,
                        soil_con_struct *, out_data_file_struct *,
                        out_data_struct *);
void   initialize_global();
int   initialize_model_state(dist_prcp_struct *, dmy_struct,
                             global_param_struct *, filep_struct, 
                             int, int, int, int, 
                             double, soil_con_struct *,
                             veg_con_struct *, lake_con_struct,
                             char **, int **);
// vic_run: int    initialize_new_storm(cell_data_struct ***, veg_var_struct ***,
//			    int, int, int, double, double);
void   initialize_snow(snow_data_struct **, int, int);
void   initialize_soil(cell_data_struct **, soil_con_struct *, veg_con_struct *, int);
void   initialize_veg(veg_var_struct **, veg_con_struct *,
                      global_param_struct *, int);
// vic_run: void   latent_heat_from_snow(double, double, double, double, double, 
//                             double, double, double *, double *, 
//                             double *, double *, double *);
// vic_run: double linear_interp(double,double,double,double,double);
cell_data_struct **make_cell_data(int, int);
dist_prcp_struct make_dist_prcp(int);
dmy_struct *make_dmy(global_param_struct *);
energy_bal_struct **make_energy_bal(int);
void make_in_and_outfiles(filep_struct *, filenames_struct *, 
                          soil_con_struct *, out_data_file_struct *);
// remove: out_data_struct *make_out_data(int);
snow_data_struct **make_snow_data(int);
veg_var_struct **make_veg_var(int);
// vic_run: void   MassRelease(double *,double *,double *,double *);
// vic_run: double maximum_unfrozen_water(double, double, double, double);
// vic_run: double modify_Ksat(double);
// mtclim_wrapper: void mtclim_init
// mtclim_wrapper: void mtclim_to_vic
void mtclim_wrapper(int, int, double, double, double, double,
                    double, double, double, double,
                    int, dmy_struct *, double *,
                    double *, double *, double *, double *, double *, double *);
// vic_run: double new_snow_density(double);
// vic_run: int    newt_raph(void (*vecfunc)(double *, double *, int, int, ...), 
//               double *, int);
void   nrerror(char *);
FILE  *open_file(char string[], char type[]);
FILE  *open_state_file(global_param_struct *, filenames_struct, int, int);
void parse_output_info(filenames_struct *, FILE *, out_data_file_struct **, out_data_struct *);
// vic_run: double penman(double, double, double, double, double, double, double);
// vic_run: void photosynth(char, double, double, double, double, double, double,
//                double, double, double, char *, double *, double *,
//                double *, double *, double *);
// CalcBlowingSnow: polint
// vic_run: void   prepare_full_energy(int, int, int, dist_prcp_struct *, 
//			   soil_con_struct *, double *, double *); 
// remove: double priestley(double, double);
// mtclim_constants_vic: int pulled_boxcar
// vic_run: int    put_data(dist_prcp_struct *, atmos_data_struct *,
//		soil_con_struct *, veg_con_struct *,
//                lake_con_struct *, out_data_file_struct *,
//		out_data_struct *, save_data_struct *,
// 	        dmy_struct *, int); 
// CalcBlowingSnow: qromb
void   read_atmos_data(FILE *, global_param_struct, int, int, double **);
double **read_forcing_data(FILE **, global_param_struct);
void   read_initial_model_state(FILE *, dist_prcp_struct *, 
                                global_param_struct *, int, int, int, 
                                soil_con_struct *, int, char *,
                                int *, lake_con_struct);
// lake: lake_con_struct read_lakeparam
void   read_snowband(FILE *, soil_con_struct *);
// remove: void   read_snowmodel(atmos_data_struct *, FILE *, int, int, int, int);
soil_con_struct read_soilparam(FILE *, char *, char *);
veg_lib_struct *read_veglib(FILE *, int *);
veg_con_struct *read_vegparam(FILE *, int, int);
// vic_run: int    redistribute_during_storm(cell_data_struct ***, veg_var_struct ***,
//				 int, int, int, double, double, double, 
//				 double *);
// remove: void   redistribute_moisture(layer_data_struct *, double *, double *,
//			     double *, double *, double *, int);
// vic_run: unsigned char redistribute_moisture_for_storm(double *, double *, double, 
//					      double, double);
// lake: rescale_snow_energy_fluxes
// lake: rescale_snow_storage
// lake: rescale_soil_veg_fluxes
// lake: rhoinit
// vic_run: double root_brent(double, double, char *, double (*Function)(double, va_list), ...);
// CalcBlowingSnow: rtnewt
// vic_run: int    runoff(cell_data_struct *, cell_data_struct *,
//              energy_bal_struct *, soil_con_struct *, double *,
//              double *, double, int, int, int, int, int);
void set_max_min_hour(double *, int, int *, int *);
// vic_run: void set_node_parameters(double *, double *, double *, double *, double *, double *,
//			 double *, double *, double *, double *, double *,
//			 double *, double *, int, int, char);
out_data_file_struct *set_output_defaults(out_data_struct *);
int set_output_var(out_data_file_struct *, int, int, out_data_struct *, char *, int, char *, int, float);
// CalcBlowingSnow: shear_stress
// vic_run: double snow_albedo(double, double, double, double, double, double, int, char);
// vic_run: double snow_density(snow_data_struct *, double, double, double, double, double);
// vic_run: int    snow_intercept(double, double, double, double, double, double,
//                      double, double, double, double, double, double, 
//                      double *, double *, double *, double *, double *, 
//                      double *, double *, double *, double *, double *, 
//                      double *, double *, double *, double *, double *, 
//                      double *, char *, int *, double *, double *, double *, 
//                      double *, double *, double *, float *,
//                      int, int, int, int, int, int, int, int,
//                      double *, double *,
//                      atmos_data_struct *, layer_data_struct *, 
//                      layer_data_struct *, soil_con_struct *, 
//                      veg_var_struct *, veg_var_struct *);
// vic_run: int    snow_melt(double, double, double, double, double *, double, double *, double, 
//		 double, double, double, double, double, double, double, 
//                 double, double, double, double, double, double, 
//                 double *, double *, double *, double *, double *, double *, 
//                 double *, double *, double *, double *, double *, double *, 
//                 int, int, int, int, snow_data_struct *, soil_con_struct *);
// mtclim_constants_vic.h: snowpack
// vic_run: double SnowPackEnergyBalance(double, va_list);
// vic_run: void   soil_carbon_balance(soil_con_struct *, energy_bal_struct *,
//                           cell_data_struct *, veg_var_struct *);
// vic_run: double soil_conductivity(double, double, double, double, double, double, double, double);
// remove: void   soil_thermal_calc(soil_con_struct *, layer_data_struct *,
//			 energy_bal_struct, double *, double *, double *,
//			 int, int);
// vic_run: double soil_thermal_eqn(double, va_list);
// vic_run: double solve_atmos_energy_bal(double Tcanopy, ...);
// vic_run: double solve_atmos_moist_bal(double , ...);
// vic_run: double solve_canopy_energy_bal(double Tfoliage, ...);
// lake: solve_lake
// vic_run: double solve_snow(char, double, double, double, double, double, double,
//                  double, double, double, double, double, double,
//                  double *, double *, double *, double *, double *,
//                  double *, double *, double *, double *, double *,
//                  double *, double *, double *, double *, double *, double *,
//                  double *, double *, double *, double *, double *, double *,
//                  double *, double *, double *, double *, double *, double *,
//                  float *, int, int, int, int, int, int, int, int, int, int *,
//                  double *, double *,
//                  dmy_struct *, atmos_data_struct *, energy_bal_struct *,
//                  layer_data_struct *, layer_data_struct *,
//                  snow_data_struct *, soil_con_struct *,
//                  veg_var_struct *, veg_var_struct *);
// remove: double solve_snow_ground_flux(double Tsurf, ...);
// vic_run: double solve_surf_energy_bal(double Tsurf, ...);
// vic_run: int    solve_T_profile(double *, double *, char *, int *, double *, double *,double *, 
//		       double *, double, double *, double *, double *,
//		       double *, double *, double *, double *, double, double *,
//		       int, int *, int, int, int, int);
// vic_run: int   solve_T_profile_implicit(double *, double *, char *, int *, double *, double *, double *,
//			       double *, double, double *, double *, double *,
//			       double *, double *, double *, double *, double, int, int *,
//			       int, int, int, int, 
//			       double *, double *, double *, double *, double *, double *, double *);
// vic_run: double StabilityCorrection(double, double, double, double, double, double);
// CalcBlowingSnow: sub_with_height
// vic_run: int    surface_fluxes(char, double, double, double, double, 
//		      double, double, double *, double *, double **,
//                      double *, double *, double *, double *, 
//                      double *, double *, double *, double *, double *,
//		      float *, int, int, int, int, int, 
//                      int, int, int, int, atmos_data_struct *, dmy_struct *, 
//                      energy_bal_struct *, global_param_struct *, 
//                      cell_data_struct *, cell_data_struct *, 
//                      snow_data_struct *, soil_con_struct *, 
//                      veg_var_struct *, veg_var_struct *, float, float, float, double *);
// vic_run: double svp(double);
// vic_run: double svp_slope(double);
// lake: temp_area
// vic_run: void transpiration(layer_data_struct *, int, int, double, double, double, 
//		   double, double, double, double, double, double, 
//		   double *, double *, double *, double *, double *,
//                 double *, float *, double *, double, double *,
//                   double, double *, double *, double *, double *);
// CalcBlowingSnow: transport_with_height
// CalcBlowingSnow" trapzd
// remove: void tridag(double *,double *,double *,double *,double *,int);
// vic_run: void tridiag(double *, double *, double *, double *, unsigned);
int update_thermal_nodes(dist_prcp_struct *, 
                         int, int, int, soil_con_struct *, veg_con_struct *);
// read_soilparam: ttrim
// read_vegparam: ttrim
void usage(char *);
void vicerror(char *);
// vic_run: double volumetric_heat_capacity(double,double,double,double);
// lake: water_energy_balance
// lake: water_under_ice
// vic_run: void wrap_compute_zwt(soil_con_struct *, cell_data_struct *);
void write_data(out_data_file_struct *, out_data_struct *, dmy_struct *, int);
// remove: void write_dist_prcp(dist_prcp_struct *);
void write_forcing_file(atmos_data_struct *, int, out_data_file_struct *, out_data_struct *);
void write_header(out_data_file_struct *, out_data_struct *, dmy_struct *, global_param_struct);
void write_layer(layer_data_struct *, int, int, 
                 double *, double *);
void write_model_state(dist_prcp_struct *, global_param_struct *, int, 
                       int, filep_struct *, soil_con_struct *, char *,
                       int *, lake_con_struct);
void write_vegvar(veg_var_struct *, int);
void zero_output_list(out_data_struct *);

