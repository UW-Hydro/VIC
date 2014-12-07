/******************************************************************************
 * @section DESCRIPTION
 *
 * Header file for vic_driver_shared routines
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

#ifndef VIC_DRIVER_SHARED_H
#define VIC_DRIVER_SHARED_H

#include <stdio.h>
#include <vic_def.h>
#include <vic_physical_constants.h>

#define VERSION "5.0 beta 2014-Dec-03"
#define SHORT_VERSION "5.0.beta"

/******************************************************************************
 * @brief    Stores forcing file input information.
 *****************************************************************************/
typedef struct {
    size_t N_ELEM; /**< number of elements per record; for LAI and ALBEDO,
                        1 element per veg tile; for others N_ELEM = 1; */
    char SIGNED;
    int SUPPLIED;
    double multiplier;
} force_type_struct;

/******************************************************************************
 * @brief    This structure records the parameters set by the forcing file
             input routines.  Those filled, are used to estimate the paramters
             needed for the model run in initialize_atmos.c.
 *****************************************************************************/
typedef struct {
    force_type_struct TYPE[N_FORCING_TYPES];
    unsigned short FORCE_DT[2];    /**< forcing file time step */
    int FORCE_ENDIAN[2];            /**< endian-ness of input file, used for
                                       DAILY_BINARY format */
    int FORCE_FORMAT[2];            /**< ASCII or BINARY */
    int FORCE_INDEX[2][N_FORCING_TYPES];
    unsigned N_TYPES[2];
} param_set_struct;

void calc_root_fractions(veg_con_struct *veg_con, soil_con_struct *soil_con);
void compute_treeline(atmos_data_struct *, dmy_struct *, double, double *,
                      bool *);
void cmd_proc(int argc, char **argv, char *globalfilename);
void compress_files(char string[]);
char*get_current_datetime(void);
void display_current_settings(int);
void free_all_vars(all_vars_struct *all_vars, int Nveg);
void free_dmy(dmy_struct **dmy);
void free_vegcon(veg_con_struct **veg_con);
double get_dist(double lat1, double long1, double lat2, double long2);
char*get_logname(const char *pathtofile);
void get_next_time_step(unsigned short *, unsigned short *, unsigned short *,
                        unsigned short *, unsigned short *, unsigned short);
void get_parameters(FILE *paramfile);
void initialize_forcing_files(void);
void initialize_filenames(void);
void initialize_fileps(void);
void initialize_global(void);
void initialize_log(void);
void initialize_options(void);
void initialize_parameters(void);
void initialize_snow(snow_data_struct **snow, size_t veg_num);
void initialize_soil(cell_data_struct **cell, soil_con_struct *soil_con,
                     size_t veg_num);
void initialize_veg(veg_var_struct **veg_var, size_t nveg);
all_vars_struct make_all_vars(size_t nveg);
cell_data_struct **make_cell_data(size_t veg_type_num);
dmy_struct *make_dmy(global_param_struct *global);
energy_bal_struct **make_energy_bal(size_t nveg);
snow_data_struct **make_snow_data(size_t nveg);
veg_var_struct **make_veg_var(size_t veg_type_num);
FILE  *open_file(char string[], char type[]);
void print_cell_data(cell_data_struct *cell, size_t nlayers, size_t nfrost,
                     size_t npet);
void print_dmy(dmy_struct *dmy);
void print_energy_bal(energy_bal_struct *eb, size_t nnodes, size_t nfronts);
void print_filenames(filenames_struct *fnames);
void print_filep(filep_struct *fp);
void print_force_type(force_type_struct *force_type);
void print_global_param(global_param_struct *gp);
void print_lake_con(lake_con_struct *lcon, size_t nlnodes);
void print_lake_var(lake_var_struct *lvar, size_t nlnodes, size_t nfronts,
                    size_t nlayers, size_t nnodes, size_t nfrost, size_t npet);
void print_layer_data(layer_data_struct *ldata, size_t nfrost);
void print_license(void);
void print_option(option_struct *option);
void print_out_data(out_data_struct *out, size_t nelem);
void print_out_data_file(out_data_file_struct *outf);
void print_param_set(param_set_struct *param_set);
void print_parameters(parameters_struct *param);
void print_save_data(save_data_struct *save);
void print_snow_data(snow_data_struct *snow);
void print_soil_con(soil_con_struct *scon, size_t nlayers, size_t nnodes,
                    size_t nfrost, size_t nbands, size_t nzwt);
void print_veg_con(veg_con_struct *vcon, size_t nroots, char blowing, char lake,
                   char carbon, size_t ncanopy);
void print_veg_lib(veg_lib_struct *vlib, char carbon);
void print_veg_var(veg_var_struct *vvar, size_t ncanopy);
void print_version(char *);
void print_usage(char *);
void setup_logging(void);
void soil_moisture_from_water_table(soil_con_struct *soil_con, size_t nlayers);
#endif
