#include <math.h>
#include <vicNl_def.h>

/*** SubRoutine Prototypes ***/

/** Data Read, Write and Handling Routines **/
atmos_data_struct *read_forcing_data(infiles_struct,int *,int);
void read_sawd(atmos_data_struct *, FILE *, int *, int, int);
void read_snowmodel(atmos_data_struct *, FILE *, int, int);
void read_atmosdata(atmos_data_struct *, FILE *, int *, int);
void read_rosemount(atmos_data_struct *, FILE *, int *, int);
void read_PILPSc(atmos_data_struct *, FILE *, int*, int);
void initialize_atmos(atmos_data_struct *, dmy_struct *, double, double, 
                      double, double, double, double *, int, int);
void write_atmosdata(atmos_data_struct *, int);

void initialize_global();

dmy_struct *make_dmy(global_param_struct);

soil_con_struct read_soilparam(FILE *);
void write_soilparam(soil_con_struct);
cell_data_struct **make_cell_data(int, int);
void initialize_soil(cell_data_struct **, soil_con_struct, int);
void write_layer(layer_data_struct *, int, int, double *);
soil_con_struct read_soilparam_arc(FILE *, char *, int *, int *, int);
double read_arcinfo_value(char *, double, double);
int read_arcinfo_info(char *, double **, double **, int **);
void read_snowband(FILE *, int, int, double, double **, double **, double **);
void read_PILPS2c(atmos_data_struct *, FILE *, int *, int);
void read_sawd_binary(atmos_data_struct *, FILE *, int *, int, int);

veg_lib_struct *read_veglib(FILE *, int *);
veg_con_struct *read_vegparam(FILE *, int, int);
void write_vegparam(veg_con_struct *);
veg_var_struct **make_veg_var(int);
void initialize_veg( veg_var_struct **, veg_con_struct *,
                     global_param_struct);
void write_vegvar(veg_var_struct, int);

/** IO and Data Storage Routines **/
void put_data(dist_prcp_struct *, atmos_data_struct *, veg_con_struct *, 
              outfiles_struct, double *, double *, dmy_struct *, int);
out_data_struct *make_out_data(int);
void write_data(out_data_struct *, outfiles_struct, dmy_struct *);
void close_files(infiles_struct, outfiles_struct, filenames_struct);
void compress_files(char string[]);
FILE *open_file(char string[], char type[]);
void open_debug();
filenames_struct cmd_proc(int argc, char *argv[]);
void usage(char *);
void check_files(infiles_struct *, filenames_struct);
global_param_struct get_global_param(filenames_struct *,FILE *);
filenames_struct make_in_and_outfiles(infiles_struct *, filenames_struct, 
                          soil_con_struct, outfiles_struct *);

/** Runoff Routines **/
void runoff(layer_data_struct *, layer_data_struct *, energy_bal_struct *, 
	    soil_con_struct, double *, double *, double *, double *, 
	    double *, double, int, int, int, int, int);
void vertical(layer_data_struct *, soil_con_struct, double, int);
void tridag(double *,double *,double *,double *,double *,int);

/** Evaporation Routines **/
double priestley(double, double);
double canopy_evap(layer_data_struct *, layer_data_struct *,
		   veg_var_struct *, veg_var_struct *, char, int, int, 
		   double, double *, double, double, double, double, 
		   double, double, double,  double, double, double, double, 
		   double *, double *, double *, double *);
void transpiration(layer_data_struct *, int, int, double, double, double, 
		   double, double, double, double, double, double, double, 
		   double, double *, double *, double *, double *, double *, 
		   double *);
double penman(double, double, double, double, double, double, double, double,
              double, double, float);
double arno_evap(layer_data_struct *, layer_data_struct *, double, double, 
		 double, double, double, double, double, double, double, 
		 double, double, double, double, double, double);

/** Radiation Routines **/
void rad_and_vpd(atmos_data_struct *, soil_con_struct, int, dmy_struct *,
		 double **, double **);
double estimate_vapor_pressure(dmy_struct *,double *,double *,double,
			       double,double,int);
double calc_trans(double, double);
double calc_netshort(double, int, double);
void calc_long_shortwave(double *, double *, double *, double,
                         double, double, double, double, double,
                         double, char, char, char);
double fltrad(double, double, double, int, double); 
double shrad(double,double,double,double,double,int,double);
double svp(double);
double svp_slope(double);
double net_out_longwave(double,double,double,double,double *);
double in_shortwave(float,int,double);

/* Spatially Distributed Precipitation Routines */
dist_prcp_struct make_dist_prcp(int, int *);
void initialize_new_storm(cell_data_struct ***, veg_var_struct ***,
			  int, int, int, double, double);
void redistribute_during_storm(cell_data_struct ***, veg_var_struct ***,
			       int, int, int, double, double);
unsigned char redistribute_moisture_for_storm(double *, double *, double, 
					      double);
unsigned char average_moisture_for_storm(double *, double *, double, double);
void write_dist_prcp(dist_prcp_struct);
void free_dist_prcp(dist_prcp_struct *, int);
double my_qromo(double,double,double,double,double,double,double,double,int,
                int*);
double my_midpnt(double,double,int,double,double,double,double,double,double,
                 int);
double funct_design(double,double,double,double,double,double,double,int);
void polint(double *,double *,int, double,double *,double *);
void nrerror(char *);
void vicerror(char *);
double *read_mu_prec(int,atmos_data_struct*,FILE*);

/** Frozen Soils Routines **/
void initialize_energy_bal(energy_bal_struct **, cell_data_struct ***,
                           soil_con_struct *, double, double *, int, int,
			   int, FILE *);
energy_bal_struct **make_energy_bal(int, int *);
double soil_conductivity(double,double,double,double);
double frozen_soil_conductivity(double,double,double,double,double);
double volumetric_heat_capacity(double,double,double);
void find_0_degree_fronts(energy_bal_struct *,layer_data_struct *,
                          double, double *, double *, int, int);
void redistribute_moisture(layer_data_struct *, double *, double *,
                           double *, double *, double *, int);
void find_sublayer_temperatures(layer_data_struct *, double *, double *,
                                double *, double, double, int, int);
void soil_thermal_calc(soil_con_struct, layer_data_struct *,
                       energy_bal_struct, double *, double *, double *,
                       int, int);
void distribute_soil_property(double *,double,double,double**,
                                 int, int, double *, double *);
double maximum_unfrozen_water(double, double, double, double);
layer_data_struct find_average_layer(layer_data_struct, layer_data_struct,
				     double, double);
double find_total_layer_moisture(layer_data_struct,double);
double find_total_layer_ice(layer_data_struct,double);
void setup_frozen_soil(soil_con_struct, layer_data_struct  *,
		       layer_data_struct  *, layer_data_struct *,
		       energy_bal_struct, int, int, int, double,
		       double *, double *, double *);
void finish_frozen_soil_calcs(energy_bal_struct *, layer_data_struct *,
			      layer_data_struct *, layer_data_struct *,
			      soil_con_struct, int, int, double, 
			      double *, double *, double *, double *, 
			      double *, double *, double *);
void solve_T_profile(double *, double *, double *, double *,double *,
		     double *, double, double *, double, double *,
		     double *, double *, double *, double *, int);
double error_solve_T_profile(double Tsurf, ...);
double error_print_solve_T_profile(double, va_list);
double soil_thermal_eqn(double, va_list);
double linear_interp(double,double,double,double,double);
double modify_Ksat(double);
void correct_precip(double *, double *, double);


/** Snow Routines **/
snow_data_struct **make_snow_data(int);
void initialize_snow(snow_data_struct **,int,FILE *);
double snow_albedo(double, double, double, double, int);
double snow_density(int, double, double, double, double, double, double);
double calc_snow_ground_flux(int, int, int, int, double, double, double, 
			     double, double, double *, double *, double *, 
			     double *, double *, energy_bal_struct *, 
			     snow_data_struct *, layer_data_struct *,
                             layer_data_struct *, soil_con_struct);
double func_snow_ground_flux(double, va_list);
double solve_snow_ground_flux(double Tsurf, ...);
double error_calc_snow_ground_flux(double Tsurf, ...);
double error_print_snow_ground_flux(double, va_list);
double f(double, double, double, double, double, double, double, double,
         double, double, int, double *, double, double, double, double *,
         double *, double *, double *, double *, double *);
double solve_snow(snow_data_struct *, layer_data_struct *, 
		  layer_data_struct *, veg_var_struct *, veg_var_struct *, 
		  dmy_struct *, energy_bal_struct *, soil_con_struct, 
		  char, char *, char *, int, int, int, int, int, int, int, 
		  int, double, double , 
		  double, double, double, double, double, double, double,
		  double, double, double, double, double, double, double,
		  double, double, double, double, double *, double *,
		  double *, double *, double *, double *, double *,
		  double *, double *, double *, double *, double *,
		  double *, double *, double *, double *);
double root_brent(double, double, double (*Function)(double, va_list), ...);
void snow_melt(soil_con_struct, int, int, double, double, double,
               snow_data_struct *, double, double, double, double, double, 
	       double, double, double, double, double, double, double, 
	       double, double, double, double *, double *, double *,
	       double *, double *, double *, double *, double *, double *);
double CalcSnowPackEnergyBalance(double Tsurf, ...);
double SnowPackEnergyBalance(double, va_list);
double ErrorSnowPackEnergyBalance(double Tsurf, ...);
double ErrorPrintSnowPackEnergyBalance(double, va_list);
double WaterDensity(double);
double StabilityCorrection(double, double, double, double, double, double);
void snow_intercept(double, double, double, double, double, double, double,
		    double, double, double, double, double, double, double,
                    double *, double *, double *, double *,
                    double *, double *, double *, double *, int, int, int);
void MassRelease(double *,double *,double *,double *);

/** Full Energy Balance Routines **/
void dist_prec(atmos_data_struct *,dist_prcp_struct *,soil_con_struct,
               veg_con_struct *,dmy_struct *,global_param_struct,
               outfiles_struct,int,int,char,char,double *,double *);
void full_energy(int,atmos_data_struct *,soil_con_struct,
                 veg_con_struct *, dist_prcp_struct *,
                 dmy_struct *,global_param_struct,int,char,
		 double *,double *,double *,double *,double *);
void store_moisture_for_debug(int,int,double *,cell_data_struct ***,
			      veg_var_struct ***,snow_data_struct **,
			      soil_con_struct);
void prepare_full_energy(int, int, int, dist_prcp_struct *,
			 soil_con_struct, double *, double *);
double calc_surf_energy_bal(char, int, int, int, int, int, int, int, int, 
			    double, double, double, double, double, 
			    double, double, double, double, double, 
			    double, double, double, double, double, 
			    double *, double *, double *, double *, 
			    atmos_data_struct *, veg_var_struct *, 
			    veg_var_struct *, energy_bal_struct *, 
			    layer_data_struct *, layer_data_struct *, 
			    soil_con_struct, dmy_struct);
double resolve_surf_energy_bal(double, char, int, int, int, int, int, int, int,
			       int, double, double, double, double, 
			       double, double, double, double, double, 
			       double, double, double, double, double, 
			       double, double *, double *, double *, 
			       double *, atmos_data_struct *, 
			       veg_var_struct *, veg_var_struct *, 
			       energy_bal_struct *, layer_data_struct *, 
			       layer_data_struct *, soil_con_struct, 
			       dmy_struct);
double func_surf_energy_bal(double, va_list);
double solve_surf_energy_bal(double Tsurf, ...);
double error_calc_surf_energy_bal(double Tsurf, ...);
double error_print_surf_energy_bal(double, va_list);
double estimate_T1(double, double, double, double, double, double, double, 
		   double, double, double, double);
void write_debug(atmos_data_struct, soil_con_struct, cell_data_struct *,
                 energy_bal_struct *, snow_data_struct *, veg_var_struct *,
                 dmy_struct, global_param_struct,
                 double, double, int, int, int, int, int, char);
void CalcAerodynamic(char, int, int, double, double, double, double,
		     double *, double *, double *, double, double *,
		     double *);
double *calc_aero_resist(char,char,int,int,double,double,double,double *,
                        double *,double *,double *);
double func_aero_resist(double,double,double,double,double);
double calc_veg_displacement(double);
double calc_veg_height(double);
double calc_veg_roughness(double);
/**double calc_air_temperature(double *, double *, int);**/
double HourlyT(int,int *,double *,int *,double *,double *);
double hermint(double, int, double *, double *, double *, double *, double *);
void hermite(int, double *, double *, double *, double *, double *);
double calc_rainonly(double,double,double,double);
void store_max_min_temp(atmos_data_struct *, double *, int *,
                        double *, int *, int, int);
void calc_water_balance_error(int, double, double, double);
void calc_energy_balance_error(int, double, double, double, double, double);
