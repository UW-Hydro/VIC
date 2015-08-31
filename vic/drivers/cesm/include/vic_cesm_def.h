/******************************************************************************
 * @section DESCRIPTION
 *
 * Header file for defining datatypes and structures for the CESM driver
 *
 * @section LICENSE
 *
 * The Variable Infiltration Capacity (VIC) macroscale hydrological model
 * Copyright (C) 2014 The Land Surface Hydrology Group, Department of Civil
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

#ifndef VIC_CESM_DEF_H
#define VIC_CESM_DEF_H

#define VIC_DRIVER "CESM"

#define MAXDIMS 10

/******************************************************************************
 * @brief    Structure to store location information for individual grid cells.
 * @details  The global and local indices show the position of the grid cell
 *           within the global and local (processor) domains. When the model is
 *           run on a single processor, the glonal and local domains are
 *           identical. The model is run over a list of cells.
 *****************************************************************************/
typedef struct {
    double latitude; /**< latitude of grid cell center */
    double longitude; /**< longitude of grid cell center */
    double area; /**< area of grid cell */
    double frac; /**< fraction of grid cell that is active */
    size_t nveg; /**< number of vegetation type according to parameter file */
    size_t global_idx; /**< index of grid cell in global list of grid cells */
    size_t io_idx; /**< index of cell in 1-D I/O arrays */
    size_t local_idx; /**< index of grid cell in local list of grid cells */
} location_struct;


/******************************************************************************
 * @brief    Structure to store local and global domain information. If the
 *           model is run on a single processor, then the two are identical.
 *****************************************************************************/
typedef struct {
    size_t ncells; /**< number of active grid cell on global domain */
    size_t n_nx; /**< size of x-index; */
    size_t n_ny; /**< size of y-index */
    location_struct *locations; /**< locations structs for local domain */
} domain_struct;

/******************************************************************************
 * @brief    Structure for netcdf file information. Initially to store
 *           information for the output files (state and history)
 *****************************************************************************/
typedef struct {
    char fname[MAXSTRING + 1];
    char c_fillvalue;
    int i_fillvalue;
    double d_fillvalue;
    float f_fillvalue;
    int nc_id;
    int band_dimid;
    int front_dimid;
    int frost_dimid;
    int layer_dimid;
    int ni_dimid;
    int nj_dimid;
    int node_dimid;
    int root_zone_dimid;
    int time_dimid;
    int veg_dimid;
    size_t band_size;
    size_t front_size;
    size_t frost_size;
    size_t layer_size;
    size_t ni_size;
    size_t nj_size;
    size_t node_size;
    size_t root_zone_size;
    size_t time_size;
    size_t veg_size;
    bool open;
} nc_file_struct;

/******************************************************************************
 * @brief    Structure for netcdf variable information
 *****************************************************************************/
typedef struct {
    char nc_var_name[MAXSTRING]; /**< variable name */
    char nc_units[MAXSTRING]; /**< variable name */
    int nc_dimids[MAXDIMS]; /**< ids of dimensions */
    int nc_counts[MAXDIMS]; /**< size of dimid */
    int nc_type; /**< netcdf type */
    int nc_aggtype; /**< aggregation type as defined in vic_def.h */
    int nc_dims; /**< number of dimensions */
    int nc_write; /**< TRUE: write to file; FALSE: don't */
} nc_var_struct;

/******************************************************************************
 * @brief    Structure for mapping the vegetation types for each grid cell as
 *           stored in VIC's veg_con_struct to a regular array.
 *****************************************************************************/
typedef struct {
    size_t nv_types;  /**< total number of vegetation types */
                      /**< size of vidx and Cv arrays */
    size_t nv_active; /**< number of active vegetation types. Because of the */
                      /**< way that VIC defines nveg, this is nveg+1 */
                      /**< (for bare soil) or nveg+2 (if the treeline option */
                      /**< is active as well) */
    int *vidx;        /**< array of indices for active vegetation types */
    double *Cv;       /**< array of fractional coverage for nc_types */
} veg_con_map_struct;

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
    double l2x_Sl_soilw;  /**< volumetric soil water */
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
} l2x_data_struct;

#endif
