/******************************************************************************
 * @section DESCRIPTION
 *
 * Header file for CESM-VIC routines
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

#ifndef VIC_DRIVER_CESM_H
#define VIC_DRIVER_CESM_H

#include <vic_driver_shared_image.h>

#define VIC_DRIVER "CESM"

#define MAXDIMS 10
// Hardcoded filenames (see bld/build_vic_namelist for more info)
#define GLOBALPARAM "vic.globalconfig.txt"
#define RPOINTER "rpointer.lnd"
#define SHR_CONST_SPVAL 1.0e30   /**< CESM missing value */

/******************************************************************************
 * @brief   CESM run types
 *****************************************************************************/
enum
{
    CESM_RUNTYPE_CLEANSTART,
    CESM_RUNTYPE_RESTART,
    CESM_RUNTYPE_BRANCH
};

/******************************************************************************
 * @brief   This structure stores clock information. See also ESMF_Clock.
 *          Order is important and any changes here must be echoed in
 *          vic_cesm_def_mod_f.F90
 *****************************************************************************/
typedef struct {
    int timestep;       // timestep in seconds
    short int current_year;  // current year
    short int current_month;  // current month
    short int current_day;  // current day
    int current_dayseconds;        // current dayseconds
    bool state_flag;      // state flag
    bool stop_flag;      // stop flag
    char calendar[MAXSTRING];      // calendar
} vic_clock;


/******************************************************************************
 * @brief   This structure stores clock information. See also ESMF_Clock.
 *          Order is important and any changes here must be echoed in
 *          vic_cesm_def_mod_f.F90
 *****************************************************************************/
typedef struct {
    char caseid[MAXSTRING];
    char casedesc[MAXSTRING];
    char starttype[MAXSTRING];
    char model_version[MAXSTRING];
    char hostname[MAXSTRING];
    char username[MAXSTRING];
} case_metadata;

/******************************************************************************
 * @brief   This structure is a c type container for the x2l fields.
 *          Order is important and any changes here must be echoed in
 *          vic_cesm_def_mod_f.F90
 *****************************************************************************/
typedef struct {
    double x2l_Sa_z;  /** bottom atm level height */
    double x2l_Sa_u;  /** bottom atm level zon wind */
    double x2l_Sa_v;  /** bottom atm level mer wind */
    double x2l_Sa_ptem;  /** bottom atm level pot temp */
    double x2l_Sa_shum;  /** bottom atm level spec hum */
    double x2l_Sa_pbot;  /** bottom atm level pressure */
    double x2l_Sa_tbot;  /** bottom atm level temp */
    double x2l_Faxa_lwdn;  /** downward lw heat flux */
    double x2l_Faxa_rainc;  /** prec: liquid "convective" */
    double x2l_Faxa_rainl;  /** prec: liquid "large scale" */
    double x2l_Faxa_snowc;  /** prec: frozen "convective" */
    double x2l_Faxa_snowl;  /** prec: frozen "large scale" */
    double x2l_Faxa_swndr;  /** sw: nir direct  downward */
    double x2l_Faxa_swvdr;  /** sw: vis direct  downward */
    double x2l_Faxa_swndf;  /** sw: nir diffuse downward */
    double x2l_Faxa_swvdf;  /** sw: vis diffuse downward */
    double x2l_Sa_co2prog;  /** bottom atm level prognostic co2 */
    double x2l_Sa_co2diag;  /** bottom atm level diagnostic co2 */
    double x2l_Faxa_bcphidry;  /** flux: Black Carbon hydrophilic dry deposition */
    double x2l_Faxa_bcphodry;  /** flux: Black Carbon hydrophobic dry deposition */
    double x2l_Faxa_bcphiwet;  /** flux: Black Carbon hydrophilic wet deposition */
    double x2l_Faxa_ocphidry;  /** flux: Organic Carbon hydrophilic dry deposition */
    double x2l_Faxa_ocphodry;  /** flux: Organic Carbon hydrophobic dry deposition */
    double x2l_Faxa_ocphiwet;  /** flux: Organic Carbon hydrophilic dry deposition */
    double x2l_Faxa_dstwet1;  /** flux: Size 1 dust -- wet deposition */
    double x2l_Faxa_dstwet2;  /** flux: Size 2 dust -- wet deposition */
    double x2l_Faxa_dstwet3;  /** flux: Size 3 dust -- wet deposition */
    double x2l_Faxa_dstwet4;  /** flux: Size 4 dust -- wet deposition */
    double x2l_Faxa_dstdry1;  /** flux: Size 1 dust -- dry deposition */
    double x2l_Faxa_dstdry2;  /** flux: Size 2 dust -- dry deposition */
    double x2l_Faxa_dstdry3;  /** flux: Size 3 dust -- dry deposition */
    double x2l_Faxa_dstdry4;  /** flux: Size 4 dust -- dry deposition */
    double x2l_Flrr_flood;  /** rtm->lnd rof (flood) flux */
    bool x2l_vars_set; /** x2l set flag */
} x2l_data_struct;

/******************************************************************************
 * @brief   This structure is a c type container for the l2x fields.
 *          Order is important and any changes here must be echoed in
 *          vic_cesm_def_mod_f.F90
 *****************************************************************************/
typedef struct {
    double l2x_Sl_t;  /**< temperature */
    double l2x_Sl_tref;  /**< 2m reference temperature */
    double l2x_Sl_qref;  /**< 2m reference specific humidity */
    double l2x_Sl_avsdr;  /**< albedo: direct , visible */
    double l2x_Sl_anidr;  /**< albedo: direct , near-ir */
    double l2x_Sl_avsdf;  /**< albedo: diffuse, visible */
    double l2x_Sl_anidf;  /**< albedo: diffuse, near-ir */
    double l2x_Sl_snowh;  /**< snow height */
    double l2x_Sl_u10;  /**< 10m wind */
    double l2x_Sl_ddvel;  /**< dry deposition velocities (optional) */
    double l2x_Sl_fv;  /**< friction velocity */
    double l2x_Sl_ram1;  /**< aerodynamical resistance */
    double l2x_Sl_logz0;  /**< log z0 */
    double l2x_Fall_taux;  /**< wind stress, zonal */
    double l2x_Fall_tauy;  /**< wind stress, meridional */
    double l2x_Fall_lat;  /**< latent heat flux */
    double l2x_Fall_sen;  /**< sensible heat flux */
    double l2x_Fall_lwup;  /**< upward longwave heat flux */
    double l2x_Fall_evap;  /**< evaporation water flux */
    double l2x_Fall_swnet;  /**< heat flux shortwave net */
    double l2x_Fall_fco2_lnd;  /**< co2 flux **For testing set to 0 */
    double l2x_Fall_flxdst1;  /**< dust flux size bin 1 */
    double l2x_Fall_flxdst2;  /**< dust flux size bin 2 */
    double l2x_Fall_flxdst3;  /**< dust flux size bin 3 */
    double l2x_Fall_flxdst4;  /**< dust flux size bin 4 */
    double l2x_Fall_flxvoc;  /**< MEGAN fluxes */
    double l2x_Flrl_rofliq;  /**< lnd->rtm input fluxes */
    double l2x_Flrl_rofice;  /**< lnd->rtm input fluxes */
    bool l2x_vars_set; /** l2x set flag */
} l2x_data_struct;

void advance_vic_time(void);
void assert_time_insync(vic_clock *vclock, dmy_struct *dmy);
void finalize_cesm_time(vic_clock *vclock);
void get_global_param(FILE *);
void initialize_cesm_time(void);
void initialize_l2x_data(void);
void initialize_vic_cesm_mpi(MPI_Fint *MPI_COMM_VIC_F);
void initialize_x2l_data(void);
void make_dummy_forcings(x2l_data_struct *x2l);
void print_case_metadata(case_metadata *cmeta);
void print_l2x_data(l2x_data_struct *l2x);
void print_vic_clock(vic_clock *vclock);
void print_x2l_data(x2l_data_struct *x2l);
void read_rpointer_file(char *fname);
unsigned short int start_type_from_char(char *start_str);
char *trimstr(char *str);
void validate_filenames(filenames_struct *filenames);
void validate_global_param(global_param_struct *global_param);
void validate_options(option_struct *options);
void vic_cesm_alloc(void);
int vic_cesm_init_mpi(int MPI_COMM_VIC_F);
int vic_cesm_init(vic_clock *vclock, case_metadata *cmeta);
int vic_cesm_final(vic_clock *vclock);
void vic_cesm_finalize(void);
int vic_cesm_run(vic_clock *vclock);
void vic_force(void);
void vic_cesm_put_data(void);
void vic_cesm_run_model(void);
void vic_cesm_start(vic_clock *vclock, case_metadata *cmeta);
void vic_initialize_albedo(void);
void vic_initialize_lwup(void);
void vic_initialize_temperature(void);
void vic_populate_model_state(char *runtype_str, dmy_struct *dmy_current);
void write_rpointer_file(char *fname);

#endif
