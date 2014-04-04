/* header file for vic_run routines */

#ifndef VIC_RUN_H
#define VIC_RUN_H

void advect_carbon_storage(double, double, lake_var_struct *, cell_data_struct *);
void advect_snow_storage(double, double, double, snow_data_struct *);
void advect_soil_veg_storage(double, double, double, double *, soil_con_struct *, veg_con_struct *, cell_data_struct *, veg_var_struct *, lake_con_struct);
double advected_sensible_heat(double, double, double, double, double);
void alblake(double, double, double *, double *, float *, float *, double, double, 
         int, int *, double, double, char *, int, double);
double arno_evap(layer_data_struct *, layer_data_struct *, double, double,
		 double, double, double, double, double, double, double, double,
		 double, double *);
unsigned char average_moisture_for_storm(double *, double *, double, double);
double calc_atmos_energy_bal(double, double, double, double, double, double,
                             double, double, double, double, double, double,
                             double, double, double, double,
                             double *, double *, double *, double *,
                             double *, double *, double *, double *, char *, int *);
double calc_density(double);
double calc_energy_balance_error(int, double, double, double, double, double);
int    calc_layer_average_thermal_props(energy_bal_struct *, layer_data_struct *,
          layer_data_struct *, layer_data_struct *,
          soil_con_struct *, int, int, double *);
void calc_Nscale_factors(char, double *, double, double, double, double,
                         dmy_struct, double *);
double calc_rainonly(double,double,double,double,double);
double calc_rc(double,double,float,double,double,double,double,char);
void calc_rc_ps(char, double, double, double, double *, double,
                double, double *, double, double, double *,
                double, double, double, double *, double *);
double calc_snow_coverage(int *, double, double, double, double, double,
                          double, double, double *, double *, double *,
                          double *, double *);
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
double calc_water_balance_error(int, double, double, double);
int    CalcAerodynamic(char, double, double, double, double, double,
	  	       double *, double *, double *, double *, double *);
double CalcBlowingSnow(double, double, int, double, double, double, double,
                       double, double, double, double, double, float,
                       float, double, int, int, float, double, double, double *);
double CalcIcePackEnergyBalance(double Tsurf, ...);
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
void colavg (double *, double *, double *, float, double *, int, double, double);
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
float dragcoeff(float, double, double);
void eddy (int, double, double * , double *, double *, double, int, double, double);
void energycalc(double *, double *, int, double, double,double *, double *, double *);
double error_calc_atmos_energy_bal(double Tcanopy, ...);
double error_calc_atmos_moist_bal(double , ...);
double error_calc_canopy_energy_bal(double Tsurf, ...);
double error_calc_surf_energy_bal(double Tsurf, ...);
double error_print_atmos_energy_bal(double, va_list);
double error_print_atmos_moist_bal(double, va_list);
double error_print_canopy_energy_bal(double, va_list);
double error_print_solve_T_profile(double, va_list);
double error_print_surf_energy_bal(double, va_list);
double error_solve_T_profile(double Tsurf, ...);
double ErrorIcePackEnergyBalance(double Tsurf, ...);
double ErrorPrintIcePackEnergyBalance(double, va_list);
double ErrorPrintSnowPackEnergyBalance(double, va_list);
double ErrorSnowPackEnergyBalance(double Tsurf, ...);
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
void   faparl(double *, double, double, double, double, double *, double *);
void   fda_heat_eqn(double *, double *, int, int, ...);
void   fdjac3(double *, double *, double *, double *, double *,
            void (*vecfunc)(double *, double *, int, int, ...),
            int);
void   find_0_degree_fronts(energy_bal_struct *, double *, double *, int);
layer_data_struct find_average_layer(layer_data_struct *, layer_data_struct *,
				     double, double);
int    full_energy(char, int, int, atmos_data_struct *, dist_prcp_struct *,
		   dmy_struct *, global_param_struct *, lake_con_struct *,
                   soil_con_struct *, veg_con_struct *);
double func_atmos_energy_bal(double, va_list);
double func_atmos_moist_bal(double, va_list);
double func_canopy_energy_bal(double, va_list);
double func_surf_energy_bal(double, va_list);
int get_depth(lake_con_struct, double, double *);
// CalcBlowingSnow: get_prob
int get_sarea(lake_con_struct, double, double *);
// CalcBlowingSnow: get_shear
// CalcBlowingSnow: get_thresh
int get_volume(lake_con_struct, double, double *);
double hiTinhib(double);
int initialize_lake(lake_var_struct *, lake_con_struct, soil_con_struct *, cell_data_struct *, double, int);
int ice_melt(double, double, double *, double, snow_data_struct *, lake_var_struct *, int, double, double, double, double, double, double, double, double, double, double, double, double, double, double, double *, double *, double *, double *, double *, double *, double *, double *, double *, double);
double IceEnergyBalance(double, va_list);
void iceform (double *,double *,double ,double,double *,int, int, double, double, double *, double *, double *, double *, double *, double);
void icerad(double,double ,double,double *, double *,double *);
int    initialize_new_storm(cell_data_struct ***, veg_var_struct ***,
			    int, int, int, double, double);
int lakeice(double *, double, double, double, double, int, 
        double, double, double *, double, double, int, dmy_struct, double *, double *, double, double);
void   latent_heat_from_snow(double, double, double, double, double,
                             double, double, double *, double *,
                             double *, double *, double *);
void latsens(double,double, double, double, double, double, double, double,
         double *, double *, double);
double linear_interp(double,double,double,double,double);
float lkdrag(float, double, double, double, double);
void   MassRelease(double *,double *,double *,double *);
double maximum_unfrozen_water(double, double, double, double);
double modify_Ksat(double);
double new_snow_density(double);
int    newt_raph(void (*vecfunc)(double *, double *, int, int, ...),
               double *, int);
double penman(double, double, double, double, double, double, double);
void photosynth(char, double, double, double, double, double, double,
                double, double, double, char *, double *, double *,
                double *, double *, double *);
// CalcBlowingSnow: polint
void   prepare_full_energy(int, int, int, dist_prcp_struct *,
			   soil_con_struct *, double *, double *);
int    put_data(dist_prcp_struct *, atmos_data_struct *,
		soil_con_struct *, veg_con_struct *,
                lake_con_struct *, out_data_file_struct *,
		out_data_struct *, save_data_struct *,
 	        dmy_struct *, int);
// CalcBlowingSnow: qromb
int    redistribute_during_storm(cell_data_struct ***, veg_var_struct ***,
				 int, int, int, double, double, double,
				 double *);
unsigned char redistribute_moisture_for_storm(double *, double *, double,
					      double, double);
void rescale_snow_energy_fluxes(double, double, snow_data_struct *, energy_bal_struct *);
void rescale_snow_storage(double, double, snow_data_struct *);
void rescale_soil_veg_fluxes(double, double, cell_data_struct *, veg_var_struct *);
void rhoinit(double *, double);
double root_brent(double, double, char *, double (*Function)(double, va_list), ...);
// CalcBlowingSnow: rtnewt
int    runoff(cell_data_struct *, cell_data_struct *,
              energy_bal_struct *, soil_con_struct *, double *,
              double *, double, int, int, int, int, int);
void set_node_parameters(double *, double *, double *, double *, double *, double *,
			 double *, double *, double *, double *, double *,
			 double *, double *, int, int, char);
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
double soil_thermal_eqn(double, va_list);
double solve_atmos_energy_bal(double Tcanopy, ...);
double solve_atmos_moist_bal(double , ...);
double solve_canopy_energy_bal(double Tfoliage, ...);
// lake: solve_lake
int solve_lake(double, double, double, double, double, double, double, double, 
        double, double, lake_var_struct *, lake_con_struct, 
        soil_con_struct, int, int, double, dmy_struct, double);
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
double specheat (double);
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
void temp_area(double, double, double, double *, double *, double *, double *, int, double *, int, double, double, double*, double *, double *);
void tracer_mixer(double *, int *, int, double*, int, double, double, double *);
void transpiration(layer_data_struct *, int, int, double, double, double,
		   double, double, double, double, double, double,
		   double *, double *, double *, double *, double *,
                   double *, float *, double *, double, double *,
                   double, double *, double *, double *, double *);
// CalcBlowingSnow: transport_with_height
// CalcBlowingSnow" trapzd
void tridia(int, double *, double *, double *, double *, double *);
void tridiag(double *, double *, double *, double *, unsigned);
double volumetric_heat_capacity(double,double,double,double);
int water_balance (lake_var_struct *, lake_con_struct, int, dist_prcp_struct *, int, int, int, double, soil_con_struct, veg_con_struct);
int water_energy_balance(int, double*, double*, int, int, double, double, double, double, double, double, double, double, double, double, double, double, double, double *, double *, double *, double*, double *, double *, double *, double, double *, double *, double *, double *, double *, double);
int water_under_ice(int, double,  double, double *, double *, double, int, double, double, double, double *, double *, double *, double *, int, double, double, double, double *);
void wrap_compute_zwt(soil_con_struct *, cell_data_struct *);

#endif
