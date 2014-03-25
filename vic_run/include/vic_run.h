/* header file for vic_run routines */
#include <math.h>
#include <vicNl_def.h>
#include <LAKE.h>

/* ordered alphabetically (case-insensitive) by function name */
// lake: advect_carbon_storage
// lake: advect_snow_storage
// lake: advect_soil_veg_storage

double advected_sensible_heat(double, double, double, double, double);
// lake: alblake
// driver: void alloc_atmos(int, atmos_data_struct **);
double arno_evap(layer_data_struct *, layer_data_struct *, double, double,
		 double, double, double, double, double, double, double, double,
		 double, double *);
// driver: unsigned char average_moisture_for_storm(double *, double *, double, double);
double calc_atmos_energy_bal(double, double, double, double, double, double,
                             double, double, double, double, double, double,
                             double, double, double, double,
                             double *, double *, double *, double *,
                             double *, double *, double *, double *, char *, int *);
// lake: calc_density
double calc_energy_balance_error(int, double, double, double, double, double);
int    calc_layer_average_thermal_props(energy_bal_struct *, layer_data_struct *,
          layer_data_struct *, layer_data_struct *,
          soil_con_struct *, int, int, double *);
// driver: void   calc_longwave(double *, double, double, double);
// driver: void   calc_netlongwave(double *, double, double, double);
// driver: double calc_netshort(double, int, double, double *);
void calc_Nscale_factors(char, double *, double, double, double, double,
                         dmy_struct, double *);
// driver: double calc_rainonly(double,double,double,double,double);
double calc_rc(double,double,float,double,double,double,double,char);
void calc_rc_ps(char, double, double, double, double *, double,
                double, double *, double, double, double *,
                double, double, double, double *, double *);
// driver: void   calc_root_fractions(veg_con_struct *, soil_con_struct *);
double calc_snow_coverage(int *, double, double, double, double, double,
                          double, double, double *, double *, double *,
                          double *, double *);
// remove: double calc_snow_ground_flux(int, int, int, int, double, double, double,
//           double, double, double *, double *, double *,
//           double *, energy_bal_struct *,
//           snow_data_struct *, layer_data_struct *,
//                             layer_data_struct *, soil_con_struct *, char *);
int    calc_soil_thermal_fluxes(int, double *, double *, char *, int *, double *, double *,
        double *, double *, double *,double *,
        double *, double *, double *,
        double *, double *, double *,
        int, int, int, int);
double calc_surf_energy_bal(double, double, double, double, double, double,
                            double, double, double, double, double, double,
                            double, double, double, double, double, double,
                            double, double, double, double, double, double,
                            double, double, double,
                            double *, double *, double *, double *, double *,
                            double *, double *, double *, double *, double *,
                            float *, int, int,
                            int, int, int, int, int, int, int, int, int, int,
                            double *, double *,
                            atmos_data_struct *, dmy_struct *,
                            energy_bal_struct *, layer_data_struct *,
                            layer_data_struct *,
                            snow_data_struct *, soil_con_struct *,
                            veg_var_struct *, veg_var_struct *, int);
// remove: double calc_trans(double, double);
// driver: double calc_veg_height(double);
double calc_water_balance_error(int, double, double, double);
int    CalcAerodynamic(char, double, double, double, double, double,
	  	       double *, double *, double *, double *, double *);
double CalcBlowingSnow(double, double, int, double, double, double, double,
                       double, double, double, double, double, float,
                       float, double, int, int, float, double, double, double *);
// lake: CalcIcePackEnergyBalance
double CalcSnowPackEnergyBalance(double Tsurf, ...);
// CalcBlowingSnow: CalcSubFlux
void canopy_assimilation(char, double, double, double, double *, double,
                         double, double *, double, double, double *,
                         double, char *, double *, double *,
                         double *, double *, double *, double *,
                         double *, double *, double *, double *);
double canopy_evap(layer_data_struct *, layer_data_struct *,
		   veg_var_struct *, veg_var_struct *, char, int, int,
		   double, double *, double, double, double, double,
		   double, double, double, double, double, double,
		   double *, double *, double *, double *, double *,
                   double *, float *, double *, double, double, double *);
// driver: void   check_files(filep_struct *, filenames_struct *);
// driver: FILE  *check_state_file(char *, dmy_struct *, global_param_struct *, int, int,
//                        int *);
// driver: void   close_files(filep_struct *, out_data_file_struct *, filenames_struct *);
// driver: filenames_struct cmd_proc(int argc, char *argv[]);
// lake: colavg
void   collect_eb_terms(energy_bal_struct, snow_data_struct, cell_data_struct,
                        int *, int *, int *, int *, int *, double, double, double,
                        int, int, double, int, int, double *, double *,
                        double *, double, out_data_struct *);
void   collect_wb_terms(cell_data_struct, veg_var_struct, snow_data_struct, lake_var_struct,
                        double, double, double, double, int, int, double, int, double *,
                        double *, out_data_struct *);
// driver: void   compress_files(char string[]);
double compute_coszen(double, double, double, dmy_struct);
// driver: void   correct_precip(double *, double, double, double, double);
void   compute_pot_evap(int, dmy_struct *, int, int, double, double , double, double, double, double **, double *);
void   compute_runoff_and_asat(soil_con_struct *, double *, double, double *, double *);
void   compute_soil_resp(int, double *, double, double, double *, double *,
                         double, double, double, double *, double *, double *);
void   compute_soil_layer_thermal_properties(layer_data_struct *, double *,
                                             double *, double *, double *,
                                             double *, double *, double *,
                                             double *, int);
// driver: void   compute_treeline(atmos_data_struct *, dmy_struct *, double, double *, char *);
double compute_zwt(soil_con_struct *, int, double);
// driver: out_data_struct *create_output_list();
double darkinhib(double);
// driver: void   display_current_settings(int, filenames_struct *, global_param_struct *);
int    dist_prec(atmos_data_struct *,dist_prcp_struct *,soil_con_struct *,
		 veg_con_struct *, lake_con_struct *,
		 dmy_struct *,global_param_struct *,
		 filep_struct *, out_data_file_struct *,
		 out_data_struct *, save_data_struct *,
		 int, int, char, char, char *, int *);
int  distribute_node_moisture_properties(double *, double *, double *,
					 double *, double *, double *,
					 double *, double *, double *,
					 double *, double *, double *, double *, double *,
					 double *, double *, double *, int, int, char);
// remove: void   distribute_soil_property(double *,double,double,
//				double **l_param,
//				int, int, double *, double *);
// lake: dragcoeff
// lake: eddy
// lake: energycalc
double error_calc_atmos_energy_bal(double Tcanopy, ...);
double error_calc_atmos_moist_bal(double , ...);
double error_calc_canopy_energy_bal(double Tsurf, ...);
// remove: double error_calc_snow_ground_flux(double Tsurf, ...);
double error_calc_surf_energy_bal(double Tsurf, ...);
double error_print_atmos_energy_bal(double, va_list);
double error_print_atmos_moist_bal(double, va_list);
double error_print_canopy_energy_bal(double, va_list);
// remove: double error_print_snow_ground_flux(double, va_list);
double error_print_solve_T_profile(double, va_list);
double error_print_surf_energy_bal(double, va_list);
double error_solve_T_profile(double Tsurf, ...);
// lake: ErrorIcePackEnergyBalance
// lake: ErrorPrintIcePackEnergyBalance
double ErrorPrintSnowPackEnergyBalance(double, va_list);
double ErrorSnowPackEnergyBalance(double Tsurf, ...);
// remove: double estimate_dew_point(double, double, double, double, double);
int estimate_layer_ice_content(layer_data_struct *, double *, double *,
			       double *, double *, double *, double *,
			       double *, double *, double *,
			       double *, double, int, int, char);
int estimate_layer_ice_content_quick_flux(layer_data_struct *, double *,
					  double, double, double, double,
					  double *, double *, double *,
					  double *, double, char);
double estimate_T1(double, double, double, double, double, double, double,
		   double, double, double, double);
double exp_interp(double,double,double,double,double);
// remove: double f(double, double, double, double, double, double, double, double,
//         double, double, int, double *, double, double, double, double *,
//         double *, double *, double *, double *, double *);
void   faparl(double *, double, double, double, double, double *, double *);
void   fda_heat_eqn(double *, double *, int, int, ...);
void   fdjac3(double *, double *, double *, double *, double *,
            void (*vecfunc)(double *, double *, int, int, ...),
            int);
void   find_0_degree_fronts(energy_bal_struct *, double *, double *, int);
layer_data_struct find_average_layer(layer_data_struct *, layer_data_struct *,
				     double, double);
// remove: void   find_sublayer_temperatures(layer_data_struct *, double *, double *,
//				  double *, double, double, int, int);
// driver: void   free_atmos(int nrecs, atmos_data_struct **atmos);
// driver: void   free_dist_prcp(dist_prcp_struct *, int);
// driver: void   free_dmy(dmy_struct **dmy);
// driver: void   free_vegcon(veg_con_struct **);
// driver: void   free_veglib(veg_lib_struct **);
// driver: void   free_out_data_files(out_data_file_struct **);
// driver: void   free_out_data(out_data_struct **);
int    full_energy(char, int, int, atmos_data_struct *, dist_prcp_struct *,
		   dmy_struct *, global_param_struct *, lake_con_struct *,
                   soil_con_struct *, veg_con_struct *);
// remove: double func_aero_resist(double,double,double,double,double);
double func_atmos_energy_bal(double, va_list);
double func_atmos_moist_bal(double, va_list);
double func_canopy_energy_bal(double, va_list);
// remove: double func_snow_ground_flux(double, va_list);
double func_surf_energy_bal(double, va_list);
// remove: double get_avg_temp(double, double, double *, double *, int);
// lake: get_depth
// driver: double get_dist(double, double, double, double);
// driver: void   get_force_type(char *, int, int *);
// driver: global_param_struct get_global_param(filenames_struct *, FILE *);
// driver: void   get_next_time_step(int *, int *, int *, int *, int *, int);
// CalcBlowingSnow: get_prob
// lake: get_sarea
// CalcBlowingSnow: get_shear
// CalcBlowingSnow: get_thresh
// lake: get_volume -- should be moved to driver
// driver: double hermint(double, int, double *, double *, double *, double *, double *);
// driver: void   hermite(int, double *, double *, double *, double *, double *);
double hiTinhib(double);
// driver: void   HourlyT(int, int, int *, double *, int *, double *, double *);
// lake: ice_melt
// lake: IceEnergyBalance
// lake: iceform
// lake: icerad
// driver: void   init_output_list(out_data_struct *, int, char *, int, float);
// driver: void   initialize_atmos(atmos_data_struct *, dmy_struct *, FILE **,
//			soil_con_struct *, out_data_file_struct *, out_data_struct *);
// driver: void   initialize_global();
// lake: initialize_lake
// driver int   initialize_model_state(dist_prcp_struct *, dmy_struct,
//			      global_param_struct *, filep_struct,
//			      int, int, int, int,
//			      double, soil_con_struct *,
//                             veg_con_struct *, lake_con_struct,
//			      char **, int **);
int    initialize_new_storm(cell_data_struct ***, veg_var_struct ***,
			    int, int, int, double, double);
// driver: void   initialize_snow(snow_data_struct **, int, int);
// driver: void   initialize_soil(cell_data_struct **, soil_con_struct *, veg_con_struct *, int);
// driver: void   initialize_veg( veg_var_struct **, veg_con_struct *,
//		       global_param_struct *, int);
// lake: lakeice
void   latent_heat_from_snow(double, double, double, double, double,
                             double, double, double *, double *,
                             double *, double *, double *);
// lake: latsens
double linear_interp(double,double,double,double,double);
// lake: lkdrag
// driver: cell_data_struct **make_cell_data(int, int);
// driver: dist_prcp_struct make_dist_prcp(int);
// driver: dmy_struct *make_dmy(global_param_struct *);
// driver: energy_bal_struct **make_energy_bal(int);
// driver: void make_in_and_outfiles(filep_struct *, filenames_struct *,
//			  soil_con_struct *, out_data_file_struct *);
// driver: out_data_struct *make_out_data(int);
// driver: snow_data_struct **make_snow_data(int);
// driver: veg_var_struct **make_veg_var(int);
void   MassRelease(double *,double *,double *,double *);
double maximum_unfrozen_water(double, double, double, double);
double modify_Ksat(double);
// driver: void mtclim_wrapper(int, int, double, double, double, double,
//                      double, double, double, double,
//                      int, dmy_struct *, double *,
//                      double *, double *, double *, double *, double *, double *);
double new_snow_density(double);
int    newt_raph(void (*vecfunc)(double *, double *, int, int, ...),
               double *, int);
// driver: void   nrerror(char *);
// driver: FILE  *open_file(char string[], char type[]);
// driver: FILE  *open_state_file(global_param_struct *, filenames_struct, int, int);
// driver: void parse_output_info(filenames_struct *, FILE *, out_data_file_struct **, out_data_struct *);
double penman(double, double, double, double, double, double, double);
void photosynth(char, double, double, double, double, double, double,
                double, double, double, char *, double *, double *,
                double *, double *, double *);
// CalcBlowingSnow: polint
void   prepare_full_energy(int, int, int, dist_prcp_struct *,
			   soil_con_struct *, double *, double *);
// remove: double priestley(double, double);
int    put_data(dist_prcp_struct *, atmos_data_struct *,
		soil_con_struct *, veg_con_struct *,
                lake_con_struct *, out_data_file_struct *,
		out_data_struct *, save_data_struct *,
 	        dmy_struct *, int);
// CalcBlowingSnow: qromb
// driver: void   read_atmos_data(FILE *, global_param_struct, int, int, double **);
// driver: double **read_forcing_data(FILE **, global_param_struct);
// driver: void   read_initial_model_state(FILE *, dist_prcp_struct *,
//				global_param_struct *, int, int, int,
//				soil_con_struct *, int, char *,
//				int *, lake_con_struct);
// lake: read_lakeparam
// driver: void   read_snowband(FILE *, soil_con_struct *);
// driver: void   read_snowmodel(atmos_data_struct *, FILE *, int, int, int, int);
// driver: soil_con_struct read_soilparam(FILE *, char *, char *);
// driver: veg_lib_struct *read_veglib(FILE *, int *);
// driver: veg_con_struct *read_vegparam(FILE *, int, int);
int    redistribute_during_storm(cell_data_struct ***, veg_var_struct ***,
				 int, int, int, double, double, double,
				 double *);
// remove: void   redistribute_moisture(layer_data_struct *, double *, double *,
//			     double *, double *, double *, int);
unsigned char redistribute_moisture_for_storm(double *, double *, double,
					      double, double);
// lake: rescale_snow_energy_fluxes
// lake: rescale_snow_storage
// lake: rescale_soil_veg_fluxes
// lake: rhoinit
double root_brent(double, double, char *, double (*Function)(double, va_list), ...);
// CalcBlowingSnow: rtnewt
int    runoff(cell_data_struct *, cell_data_struct *,
              energy_bal_struct *, soil_con_struct *, double *,
              double *, double, int, int, int, int, int);
// driver: void set_max_min_hour(double *, int, int *, int *);
void set_node_parameters(double *, double *, double *, double *, double *, double *,
			 double *, double *, double *, double *, double *,
			 double *, double *, int, int, char);
// driver: out_data_file_struct *set_output_defaults(out_data_struct *);
// driver: int set_output_var(out_data_file_struct *, int, int, out_data_struct *, char *, int, char *, int, float);
// CalcBlowingSnow: shear_stress
double snow_albedo(double, double, double, double, double, double, int, char);
double snow_density(snow_data_struct *, double, double, double, double, double);
int    snow_intercept(double, double, double, double, double, double,
                      double, double, double, double, double, double,
                      double *, double *, double *, double *, double *,
                      double *, double *, double *, double *, double *,
                      double *, double *, double *, double *, double *,
                      double *, char *, int *, double *, double *, double *,
                      double *, double *, double *, float *,
                      int, int, int, int, int, int, int, int,
                      double *, double *,
                      atmos_data_struct *, layer_data_struct *,
                      layer_data_struct *, soil_con_struct *,
                      veg_var_struct *, veg_var_struct *);
int    snow_melt(double, double, double, double, double *, double, double *, double,
		 double, double, double, double, double, double, double,
                 double, double, double, double, double, double,
                 double *, double *, double *, double *, double *, double *,
                 double *, double *, double *, double *, double *, double *,
                 int, int, int, int, snow_data_struct *, soil_con_struct *);
double SnowPackEnergyBalance(double, va_list);
void   soil_carbon_balance(soil_con_struct *, energy_bal_struct *,
                           cell_data_struct *, veg_var_struct *);
double soil_conductivity(double, double, double, double, double, double, double, double);
// remove: void   soil_thermal_calc(soil_con_struct *, layer_data_struct *,
//			 energy_bal_struct, double *, double *, double *,
//			 int, int);
double soil_thermal_eqn(double, va_list);
double solve_atmos_energy_bal(double Tcanopy, ...);
double solve_atmos_moist_bal(double , ...);
double solve_canopy_energy_bal(double Tfoliage, ...);
// lake: solve_lake
double solve_snow(char, double, double, double, double, double, double,
                  double, double, double, double, double, double,
                  double *, double *, double *, double *, double *,
                  double *, double *, double *, double *, double *,
                  double *, double *, double *, double *, double *, double *,
                  double *, double *, double *, double *, double *, double *,
                  double *, double *, double *, double *, double *, double *,
                  float *, int, int, int, int, int, int, int, int, int, int *,
                  double *, double *,
                  dmy_struct *, atmos_data_struct *, energy_bal_struct *,
                  layer_data_struct *, layer_data_struct *,
                  snow_data_struct *, soil_con_struct *,
                  veg_var_struct *, veg_var_struct *);
// remove: double solve_snow_ground_flux(double Tsurf, ...);
double solve_surf_energy_bal(double Tsurf, ...);
int    solve_T_profile(double *, double *, char *, int *, double *, double *,double *,
		       double *, double, double *, double *, double *,
		       double *, double *, double *, double *, double, double *,
		       int, int *, int, int, int, int);
int   solve_T_profile_implicit(double *, double *, char *, int *, double *, double *, double *,
			       double *, double, double *, double *, double *,
			       double *, double *, double *, double *, double, int, int *,
			       int, int, int, int,
			       double *, double *, double *, double *, double *, double *, double *);
// lake: specheat
double StabilityCorrection(double, double, double, double, double, double);
// CalcBlowingSnow: sub_with_height
int    surface_fluxes(char, double, double, double, double,
		      double, double, double *, double *, double **,
                      double *, double *, double *, double *,
                      double *, double *, double *, double *, double *,
		      float *, int, int, int, int, int,
                      int, int, int, int, atmos_data_struct *, dmy_struct *,
                      energy_bal_struct *, global_param_struct *,
                      cell_data_struct *, cell_data_struct *,
                      snow_data_struct *, soil_con_struct *,
                      veg_var_struct *, veg_var_struct *, float, float, float, double *);
double svp(double);
double svp_slope(double);
// lake: temp_area
// lake: tracer_mixer
void transpiration(layer_data_struct *, int, int, double, double, double,
		   double, double, double, double, double, double,
		   double *, double *, double *, double *, double *,
                   double *, float *, double *, double, double *,
                   double, double *, double *, double *, double *);
// CalcBlowingSnow: transport_with_height
// CalcBlowingSnow" trapzd
// remove: void tridag(double *,double *,double *,double *,double *,int);
// lake: tridiag
void tridiag(double *, double *, double *, double *, unsigned);
// driver: int update_thermal_nodes(dist_prcp_struct *,
//			  int, int, int, soil_con_struct *, veg_con_struct *);
// driver: void usage(char *);
// driver: void   vicerror(char *);
double volumetric_heat_capacity(double,double,double,double);
// lake: water_balance
// lake: water_energy_balance
// lake: water_under_ice
void wrap_compute_zwt(soil_con_struct *, cell_data_struct *);
// driver: void write_data(out_data_file_struct *, out_data_struct *, dmy_struct *, int);
// driver: void write_dist_prcp(dist_prcp_struct *);
// driver: void write_forcing_file(atmos_data_struct *, int, out_data_file_struct *, out_data_struct *);
// driver: void write_header(out_data_file_struct *, out_data_struct *, dmy_struct *, global_param_struct);
// driver: void write_layer(layer_data_struct *, int, int,
//                 double *, double *);
// driver: void write_model_state(dist_prcp_struct *, global_param_struct *, int,
//		       int, filep_struct *, soil_con_struct *, char *,
//		       int *, lake_con_struct);
// driver: void write_vegvar(veg_var_struct *, int);
// driver: void zero_output_list(out_data_struct *);

