/******************************************************************************
 * @section DESCRIPTION
 *
 * Header file for vic_run routines
 *
 * @section LICENSE
 *
 * The Variable Infiltration Capacity (VIC) macroscale hydrological model
 * Copyright (C) 2016 The Computational Hydrology Group, Department of Civil
 * and Environmental Engineering, University of Washington.
 *
 * The VIC model is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 *****************************************************************************/

#ifndef VIC_RUN_H
#define VIC_RUN_H

#include <vic_def.h>

void advect_carbon_storage(double, double, lake_var_struct *,
                           cell_data_struct *);
void advect_snow_storage(double, double, double, snow_data_struct *);
void advect_soil_veg_storage(double, double, double, double *,
                             soil_con_struct *, veg_con_struct *,
                             cell_data_struct *, veg_var_struct *,
                             lake_con_struct);
double advected_sensible_heat(double, double, double, double, double);
void alblake(double, double, double *, double *, double *, double *, double,
             double, double, unsigned int *, double, bool *, unsigned short int,
             double);
double arno_evap(layer_data_struct *, double, double, double, double, double,
                 double, double, double, double, double, double *);
bool assert_close_double(double x, double y, double rtol, double abs_tol);
bool assert_close_float(float x, float y, float rtol, float abs_tol);
double calc_atmos_energy_bal(double, double, double, double, double, double,
                             double, double, double, double, double, double,
                             double, double *, double *, double *, double *,
                             double *, double *, bool *, unsigned int*);
double calc_density(double);
void calc_gridcell_avg_albedo(double *, double, size_t, bool,
                              energy_bal_struct **, snow_data_struct **,
                              veg_con_struct *, soil_con_struct *);
double calc_latent_heat_of_sublimation(double temp);
double calc_latent_heat_of_vaporization(double temp);
int calc_layer_average_thermal_props(energy_bal_struct *, layer_data_struct *,
                                     soil_con_struct *, size_t, double *);
double calc_outgoing_longwave(double temp, double emis);
double calc_scale_height(double tair, double elevation);
double calc_sensible_heat(double atmos_density, double t1, double t0,
                          double Ra);
void calc_Nscale_factors(bool, double *, double, double, double *);
double calc_rainonly(double, double, double, double);
double calc_rc(double, double, double, double, double, double, double, char);
void calc_rc_ps(char, double, double, double, double *, double, double,
                double *, double, double, double *, double, double, double,
                double *, double *);
double calc_snow_coverage(bool *, double, double, double, double, double,
                          double, double, double *, double, double *, double *,
                          double *);
int calc_soil_thermal_fluxes(int, double *, double *, char *, unsigned int *,
                             double *, double *, double *, double *, double *,
                             double *, double *, double *, double *, double *,
                             double *, int, int, int);
double calc_surf_energy_bal(double, double, double, double, double, double,
                            double, double, double, double, double, double,
                            double, double, double, double, double, double,
                            double, double, double, double, double, double *,
                            double *, double *, double *, double *, double *,
                            double, double *, double *, double, double *,
                            double *, int, int, size_t, size_t, double, size_t,
                            unsigned short int, int, unsigned short int,
                            double *, double *, force_data_struct *,
                            dmy_struct *, energy_bal_struct *,
                            layer_data_struct *, snow_data_struct *,
                            soil_con_struct *, veg_var_struct *);
double calc_veg_displacement(double);
double calc_veg_height(double);
double calc_veg_roughness(double);
int CalcAerodynamic(bool, double, double, double, double, double, double *,
                    double *, double *, double *, double *);
double CalcBlowingSnow(double, double, unsigned int, double, double, double,
                       double, double, double, double, double, double, double,
                       double, int, int, double, double, double, double *);
double CalcIcePackEnergyBalance(double Tsurf, ...);
double CalcSnowPackEnergyBalance(double Tsurf, ...);
double CalcSubFlux(double EactAir, double es, double Zrh, double AirDens,
                   double utshear, double ushear, double fe, double Tsnow,
                   double Tair, double U10, double Zo_salt, double F,
                   double *Transport);
void canopy_assimilation(char, double, double, double, double *, double, double,
                         double *, double, double, double *, double, char *,
                         double *, double *, double *, double *, double *,
                         double *, double *, double *, double *, double *);
double canopy_evap(layer_data_struct *, veg_var_struct *, bool,
                   unsigned short int, double *, double, double, double, double,
                   double, double, double, double, double *, double *, double *,
                   double *, double *, double *, double, double, double *);
void colavg(double *, double *, double *, double, double *, int, double,
            double);
double compute_coszen(double, double, double, unsigned short int, unsigned int);
void compute_derived_lake_dimensions(lake_var_struct *, lake_con_struct);
void compute_pot_evap(size_t, double, double, double, double, double, double,
                      double, double, double, double *, char, double, double,
                      double, double *);
void compute_runoff_and_asat(soil_con_struct *, double *, double, double *,
                             double *);
double calc_Q12(double, double, double, double, double);
void compute_soil_resp(int, double *, double, double, double *, double *,
                       double, double, double, double *, double *, double *);
void compute_soil_layer_thermal_properties(layer_data_struct *, double *,
                                           double *, double *, double *,
                                           double *, double *, double *,
                                           double *, size_t);
double compute_zwt(soil_con_struct *, int, double);
void correct_precip(double *, double, double, double, double);
double darkinhib(double);
int distribute_node_moisture_properties(double *, double *, double *, double *,
                                        double *, double *, double *, double *,
                                        double *, double *, double *, double *,
                                        double *, double *, double *, double *,
                                        double *, int, int, char);
void eddy(int, double, double *, double *, double, int, double, double);
void energycalc(double *, double *, int, double, double, double *, double *,
                double *);
double error_calc_atmos_energy_bal(double Tcanopy, ...);
double error_calc_atmos_moist_bal(double, ...);
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
int ErrorPrintSnowPackEnergyBalance(double, va_list);
int ErrorSnowPackEnergyBalance(double Tsurf, ...);
void estimate_frost_temperature_and_depth(double ***, double **, double *,
                                          double *, double *, double *, double,
                                          size_t, size_t);
int estimate_layer_ice_content(layer_data_struct *, double ***, double **,
                               double *, double *, double *, double *, double *,
                               size_t, size_t, char);
int estimate_layer_temperature(layer_data_struct *, double ***, double **,
                               double *, double *, size_t, size_t);
int estimate_layer_temperature_quick_flux(layer_data_struct *, double *, double,
                                          double, double, double);
int estimate_layer_ice_content_quick_flux(layer_data_struct *, double *,
                                          double *, double *, double *,
                                          double *, double, char);
double estimate_T1(double, double, double, double, double, double, double,
                   double, double, double);
void faparl(double *, double, double, double, double, double *, double *);
void fda_heat_eqn(double *, double *, int, int, ...);
void fdjac3(double *, double *, double *, double *, double *, void (*vecfunc)(
                double *, double *, int, int, ...), int);
void find_0_degree_fronts(energy_bal_struct *, double *, double *, int);
void free_2d_double(size_t *shape, double **array);
void free_3d_double(size_t *shape, double ***array);
double func_atmos_energy_bal(double, va_list);
double func_atmos_moist_bal(double, va_list);
double func_canopy_energy_bal(double, va_list);
double func_surf_energy_bal(double, va_list);
double (*funcd)(double z, double es, double Wind, double AirDens, double ZO,
                double EactAir, double F, double hsalt, double phi_r,
                double ushear,
                double Zrh);
int get_depth(lake_con_struct, double, double *);
double get_prob(double Tair, double Age, double SurfaceLiquidWater, double U10);
int get_sarea(lake_con_struct, double, double *);
void get_shear(double x, double *f, double *df, double Ur, double Zr);
double get_thresh(double Tair, double SurfaceLiquidWater, double Zo_salt);
int get_volume(lake_con_struct, double, double *);
double hiTinhib(double);
int ice_melt(double, double, double *, double, snow_data_struct *,
             lake_var_struct *, double, double, double, double, double, double,
             double, double, double, double, double, double, double, double,
             double, double *, double *, double *, double *, double *, double *,
             double *, double *, double *);
double IceEnergyBalance(double, va_list);
void iceform(double *, double *, double, double, double *, int, double, double,
             double, double *, double *, double *, double *, double);
void icerad(double, double, double, double *, double *, double *);
void initialize_lake(lake_var_struct *, lake_con_struct, soil_con_struct *,
                     cell_data_struct *, bool);
int lakeice(double, double, double, double, double, double *, double, double *,
            double *, double, double);
void latent_heat_from_snow(double, double, double, double, double, double,
                           double, double *, double *, double *, double *,
                           double *);
void latsens(double, double, double, double, double, double, double, double,
             double *, double *, double);
double linear_interp(double, double, double, double, double);
double lkdrag(double, double, double, double, double);
void malloc_2d_double(size_t *shape, double ***array);
void malloc_3d_double(size_t *shape, double ****array);
void MassRelease(double *, double *, double *, double *);
double maximum_unfrozen_water(double, double, double, double);
double new_snow_density(double);
int newt_raph(void (*vecfunc)(double *, double *, int, int,
                              ...), double *, int);
double penman(double, double, double, double, double, double, double);
void photosynth(char, double, double, double, double, double, double, double,
                double, double, char *, double *, double *, double *, double *,
                double *);
void polint(double xa[], double ya[], int n, double x, double *y, double *dy);
void prepare_full_energy(cell_data_struct *, energy_bal_struct *,
                         soil_con_struct *, double *, double *);
double qromb(
    double (*sub_with_height)(), double es, double Wind, double AirDens, double ZO, double EactAir, double F, double hsalt, double phi_r, double ushear, double Zrh, double a, double b);
void rescale_snow_energy_fluxes(double, double, snow_data_struct *,
                                energy_bal_struct *);
void rescale_snow_storage(double, double, snow_data_struct *);
void rescale_soil_veg_fluxes(double, double, cell_data_struct *,
                             veg_var_struct *);
void rhoinit(double *, double);
double root_brent(double, double, double (*Function)(double, va_list), ...);
double rtnewt(double x1, double x2, double xacc, double Ur, double Zr);
int runoff(cell_data_struct *, energy_bal_struct *, soil_con_struct *, double,
           double *, int);
void set_node_parameters(double *, double *, double *, double *, double *,
                         double *, double *, double *, double *, double *,
                         double *, int, int);
void shear_stress(double U10, double ZO, double *ushear, double *Zo_salt,
                  double utshear);
double snow_albedo(double, double, double, double, double, int, bool);
double snow_density(snow_data_struct *, double, double, double, double);
int snow_intercept(double, double, double, double, double, double, double,
                   double, double, double, double *, double *, double *,
                   double *, double *, double *, double *, double *, double *,
                   double *, double *, double *, double *, double *, double *,
                   double *, bool *, unsigned int *, double *, double *,
                   double *, double *, double *, double *, double *, int, int,
                   int, int, int, unsigned short int, double *, double *,
                   force_data_struct *, layer_data_struct *, soil_con_struct *,
                   veg_var_struct *);
int snow_melt(double, double, double, double, double *, double, double *,
              double, double, double, double, double, double, double, double,
              double, double, double, double, double, double *, double *,
              double *, double *, double *, double *, double *, double *,
              double *, double *, double *, double *, int, int, int,
              snow_data_struct *);
double SnowPackEnergyBalance(double, va_list);
void soil_carbon_balance(soil_con_struct *, energy_bal_struct *,
                         cell_data_struct *, veg_var_struct *);
double soil_conductivity(double, double, double, double, double, double, double,
                         double);
double soil_thermal_eqn(double, va_list);
double solve_atmos_energy_bal(double Tcanopy, ...);
double solve_atmos_moist_bal(double, ...);
double solve_canopy_energy_bal(double Tfoliage, ...);
int solve_lake(double, double, double, double, double, double, double, double,
               double, double, lake_var_struct *, soil_con_struct, double,
               double, dmy_struct, double);
double solve_snow(char, double, double, double, double, double, double, double,
                  double, double, double *, double *, double *, double *,
                  double *, double *, double *, double *, double *, double *,
                  double *, double *, double *, double *, double *, double *,
                  double *, double *, double *, double *, double *, double *,
                  double *, double *, double *, double *, double *, double *,
                  int, size_t, unsigned short int, unsigned short int, double,
                  size_t, int, int *, double *, double *, dmy_struct *,
                  force_data_struct *, energy_bal_struct *, layer_data_struct *,
                  snow_data_struct *, soil_con_struct *, veg_var_struct *);
double solve_surf_energy_bal(double Tsurf, ...);
int solve_T_profile(double *, double *, char *, unsigned int *, double *,
                    double *, double *, double *, double, double *, double *,
                    double *, double *, double *, double *, double *, double,
                    int, int *, int, int, int);
int solve_T_profile_implicit(double *, double *, char *, unsigned int *,
                             double *, double *, double *, double *, double,
                             double *, double *, double *, double *, double *,
                             double *, double *, double, int, int *, int, int,
                             double *, double *, double *, double *, double *,
                             double *, double *);
double specheat(double);
double StabilityCorrection(double, double, double, double, double, double);
double sub_with_height(double z, double es, double Wind, double AirDens,
                       double ZO, double EactAir, double F, double hsalt,
                       double phi_r, double ushear, double Zrh);
int surface_fluxes(bool, double, double, double, double, double *, double *,
                   double *, double *, double *, double *, double *, double *,
                   double *, double *, double *, double *, double *, size_t,
                   size_t, unsigned short int, double, unsigned short int,
                   unsigned short int, force_data_struct *, dmy_struct *,
                   energy_bal_struct *, global_param_struct *,
                   cell_data_struct *, snow_data_struct *, soil_con_struct *,
                   veg_var_struct *, double, double, double, double *);
double svp(double);
double svp_slope(double);
void temp_area(double, double, double, double *, double *, double *, double *,
               double, double *, int, double, double, double *, double *,
               double *);
void tracer_mixer(double *, int *, double *, int, double, double, double *);
void transpiration(layer_data_struct *, veg_var_struct *, unsigned short int,
                   double, double, double, double, double, double, double,
                   double, double *, double *, double *, double *, double *,
                   double *, double, double, double *);
double transport_with_height(double z, double es, double Wind, double AirDens,
                             double ZO, double EactAir, double F, double hsalt,
                             double phi_r, double ushear, double Zrh);
double trapzd(
    double (*funcd)(), double es, double Wind, double AirDens, double ZO, double EactAir, double F, double hsalt, double phi_r, double ushear, double Zrh, double a, double b, int n);
void tridia(int, double *, double *, double *, double *, double *);
void tridiag(double *, double *, double *, double *, unsigned int);
int vic_run(force_data_struct *, all_vars_struct *, dmy_struct *,
            global_param_struct *, lake_con_struct *, soil_con_struct *,
            veg_con_struct *, veg_lib_struct *);
double volumetric_heat_capacity(double, double, double, double);
int water_balance(lake_var_struct *, lake_con_struct, double, all_vars_struct *,
                  int, int, double, soil_con_struct, veg_con_struct);
int water_energy_balance(int, double *, double *, double, double, double,
                         double, double, double, double, double, double, double,
                         double, double, double, double *, double *, double *,
                         double *, double *, double *, double *, double,
                         double *, double *, double *, double *, double *,
                         double);
int water_under_ice(int, double, double, double *, double *, double, int,
                    double, double, double, double *, double *, double *,
                    double *, int, double, double, double, double *);
void wrap_compute_zwt(soil_con_struct *, cell_data_struct *);
void write_layer(layer_data_struct *, int, double *);
void write_vegvar(veg_var_struct *, int);

#endif
