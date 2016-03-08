/******************************************************************************
 * @section DESCRIPTION
 *
 * MPI support routines for VIC
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

#include <vic_driver_shared_image.h>

/******************************************************************************
 * @brief   Initialize MPI functionality
 *****************************************************************************/
void
initialize_mpi(void)
{
    extern MPI_Datatype mpi_global_struct_type;
    extern MPI_Datatype mpi_filenames_struct_type;
    extern MPI_Datatype mpi_location_struct_type;
    extern MPI_Datatype mpi_nc_file_struct_type;
    extern MPI_Datatype mpi_option_struct_type;
    extern MPI_Datatype mpi_param_struct_type;
    extern int          mpi_rank;
    extern int          mpi_size;
    int                 status;

    // get MPI mpi_rank and mpi_size
    status = MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    if (status != MPI_SUCCESS) {
        log_err("MPI error in initialize_mpi(): %d\n", status);
    }

    status = MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    if (status != MPI_SUCCESS) {
        log_err("MPI error in initialize_mpi(): %d\n", status);
    }

    // initialize MPI data structures
    create_MPI_global_struct_type(&mpi_global_struct_type);
    create_MPI_filenames_struct_type(&mpi_filenames_struct_type);
    create_MPI_location_struct_type(&mpi_location_struct_type);
    create_MPI_nc_file_struct_type(&mpi_nc_file_struct_type);
    create_MPI_option_struct_type(&mpi_option_struct_type);
    create_MPI_param_struct_type(&mpi_param_struct_type);
}

/******************************************************************************
 * @brief   Create an MPI_Datatype that represents the global_param_struct
 * @details This allows MPI operations in which the entire global_param_struct
 *          can be treated as an MPI_Datatype. NOTE: This function needs to be
 *          kept in-sync with the global_param_struct data type in vic_def.h.
 *
 * @param mpi_type MPI_Datatype that can be used in MPI operations
 *****************************************************************************/
void
create_MPI_global_struct_type(MPI_Datatype *mpi_type)
{
    int           nitems; // number of elements in struct
    int           status;
    int          *blocklengths;
    size_t        i;
    MPI_Aint     *offsets;
    MPI_Datatype *mpi_types;

    // nitems has to equal the number of elements in global_param_struct
    nitems = 33;
    blocklengths = malloc(nitems * sizeof(*blocklengths));
    if (blocklengths == NULL) {
        log_err("Memory allocation error in create_MPI_global_struct_type().")
    }

    offsets = malloc(nitems * sizeof(*offsets));
    if (offsets == NULL) {
        log_err("Memory allocation error in create_MPI_global_struct_type().")
    }

    mpi_types = malloc(nitems * sizeof(*mpi_types));
    if (mpi_types == NULL) {
        log_err("Memory allocation error in create_MPI_global_struct_type().")
    }

    // most of the elements in global_param_struct are not arrays. Use 1 as
    // the default block length and reset as needed
    for (i = 0; i < (size_t) nitems; i++) {
        blocklengths[i] = 1;
    }

    // reset i
    i = 0;

    // double wind_h;
    offsets[i] = offsetof(global_param_struct, wind_h);
    mpi_types[i++] = MPI_DOUBLE;

    // double resolution;
    offsets[i] = offsetof(global_param_struct, resolution);
    mpi_types[i++] = MPI_DOUBLE;

    // double dt;
    offsets[i] = offsetof(global_param_struct, dt);
    mpi_types[i++] = MPI_DOUBLE;

    // double snow_dt;
    offsets[i] = offsetof(global_param_struct, snow_dt);
    mpi_types[i++] = MPI_DOUBLE;

    // double runoff_dt;
    offsets[i] = offsetof(global_param_struct, runoff_dt);
    mpi_types[i++] = MPI_DOUBLE;

    // double atmos_dt;
    offsets[i] = offsetof(global_param_struct, atmos_dt);
    mpi_types[i++] = MPI_DOUBLE;

    // double out_dt;
    offsets[i] = offsetof(global_param_struct, out_dt);
    mpi_types[i++] = MPI_DOUBLE;

    // size_t model_steps_per_day;
    offsets[i] = offsetof(global_param_struct, model_steps_per_day);
    mpi_types[i++] = MPI_AINT;

    // size_t snow_steps_per_day;
    offsets[i] = offsetof(global_param_struct, snow_steps_per_day);
    mpi_types[i++] = MPI_AINT;

    // size_t runoff_steps_per_day;
    offsets[i] = offsetof(global_param_struct, runoff_steps_per_day);
    mpi_types[i++] = MPI_AINT;

    // size_t atmos_steps_per_day;
    offsets[i] = offsetof(global_param_struct, atmos_steps_per_day);
    mpi_types[i++] = MPI_AINT;

    // size_t output_steps_per_day;
    offsets[i] = offsetof(global_param_struct, output_steps_per_day);
    mpi_types[i++] = MPI_AINT;

    // unsigned short endday;
    offsets[i] = offsetof(global_param_struct, endday);
    mpi_types[i++] = MPI_UNSIGNED_SHORT;

    // unsigned short endmonth;
    offsets[i] = offsetof(global_param_struct, endmonth);
    mpi_types[i++] = MPI_UNSIGNED_SHORT;

    // unsigned short endyear;
    offsets[i] = offsetof(global_param_struct, endyear);
    mpi_types[i++] = MPI_UNSIGNED_SHORT;

    // unsigned short forceday[2];
    offsets[i] = offsetof(global_param_struct, forceday);
    blocklengths[i] = 2;
    mpi_types[i++] = MPI_UNSIGNED_SHORT;

    // unsigned int forcesec[2];
    offsets[i] = offsetof(global_param_struct, forcesec);
    blocklengths[i] = 2;
    mpi_types[i++] = MPI_UNSIGNED;

    // unsigned short forcemonth[2];
    offsets[i] = offsetof(global_param_struct, forcemonth);
    blocklengths[i] = 2;
    mpi_types[i++] = MPI_UNSIGNED_SHORT;

    // unsigned short forceoffset[2];
    offsets[i] = offsetof(global_param_struct, forceoffset);
    blocklengths[i] = 2;
    mpi_types[i++] = MPI_UNSIGNED_SHORT;

    // unsigned int forceskip[2];
    offsets[i] = offsetof(global_param_struct, forceskip);
    blocklengths[i] = 2;
    mpi_types[i++] = MPI_UNSIGNED;

    // unsigned short int forceyear[2];
    offsets[i] = offsetof(global_param_struct, forceyear);
    blocklengths[i] = 2;
    mpi_types[i++] = MPI_UNSIGNED_SHORT;

    // size_t nrecs;
    offsets[i] = offsetof(global_param_struct, nrecs);
    mpi_types[i++] = MPI_AINT;

    // unsigned short int skipyear;
    offsets[i] = offsetof(global_param_struct, skipyear);
    mpi_types[i++] = MPI_UNSIGNED_SHORT;

    // unsigned short int startday;
    offsets[i] = offsetof(global_param_struct, startday);
    mpi_types[i++] = MPI_UNSIGNED_SHORT;

    // unsigned int startsec;
    offsets[i] = offsetof(global_param_struct, startsec);
    mpi_types[i++] = MPI_UNSIGNED;

    // unsigned short int startmonth;
    offsets[i] = offsetof(global_param_struct, startmonth);
    mpi_types[i++] = MPI_UNSIGNED_SHORT;

    // unsigned short int startyear;
    offsets[i] = offsetof(global_param_struct, startyear);
    mpi_types[i++] = MPI_UNSIGNED_SHORT;

    // unsigned short int stateday;
    offsets[i] = offsetof(global_param_struct, stateday);
    mpi_types[i++] = MPI_UNSIGNED_SHORT;

    // unsigned short int statemonth;
    offsets[i] = offsetof(global_param_struct, statemonth);
    mpi_types[i++] = MPI_UNSIGNED_SHORT;

    // unsigned short int stateyear;
    offsets[i] = offsetof(global_param_struct, stateyear);
    mpi_types[i++] = MPI_UNSIGNED_SHORT;

    // unsigned short int calendar;
    offsets[i] = offsetof(global_param_struct, calendar);
    mpi_types[i++] = MPI_UNSIGNED_SHORT;

    // unsigned short int time_units;
    offsets[i] = offsetof(global_param_struct, time_units);
    mpi_types[i++] = MPI_UNSIGNED_SHORT;

    // double time_origin_num;
    offsets[i] = offsetof(global_param_struct, time_origin_num);
    mpi_types[i++] = MPI_DOUBLE;

    // make sure that the we have the right number of elements
    if (i != (size_t) nitems) {
        log_err("Miscount in create_MPI_global_struct_type(): "
                "%zd not equal to %d\n", i, nitems);
    }

    status = MPI_Type_create_struct(nitems, blocklengths, offsets, mpi_types,
                                    mpi_type);
    if (status != MPI_SUCCESS) {
        log_err("MPI error in create_MPI_global_struct_type(): %d\n", status);
    }

    status = MPI_Type_commit(mpi_type);
    if (status != MPI_SUCCESS) {
        log_err("MPI error in create_MPI_global_struct_type(): %d\n", status);
    }

    // cleanup
    free(blocklengths);
    free(offsets);
    free(mpi_types);
}

/******************************************************************************
 * @brief   Create an MPI_Datatype that represents the filenames_struct
 * @details This allows MPI operations in which the entire filenames_struct
 *          can be treated as an MPI_Datatype. NOTE: This function needs to be
 *          kept in-sync with the global_param_struct data type in vic_def.h.
 *
 * @param mpi_type MPI_Datatype that can be used in MPI operations
 *****************************************************************************/
void
create_MPI_filenames_struct_type(MPI_Datatype *mpi_type)
{
    int           nitems; // number of elements in struct
    int           status;
    int          *blocklengths;
    size_t        i;
    MPI_Aint     *offsets;
    MPI_Datatype *mpi_types;

    // nitems has to equal the number of elements in filenames_struct
    nitems = 14;
    blocklengths = malloc(nitems * sizeof(*blocklengths));
    if (blocklengths == NULL) {
        log_err("Memory allocation error in create_MPI_filenames_struct_type().")
    }

    offsets = malloc(nitems * sizeof(*offsets));
    if (offsets == NULL) {
        log_err("Memory allocation error in create_MPI_filenames_struct_type().")
    }

    mpi_types = malloc(nitems * sizeof(*mpi_types));
    if (mpi_types == NULL) {
        log_err("Memory allocation error in create_MPI_filenames_struct_type().")
    }

    // most of the elements in filenames_struct are character arrays. Use
    // MAXSTRING as the default block length and reset as needed
    for (i = 0; i < (size_t) nitems; i++) {
        blocklengths[i] = MAXSTRING;
    }

    // reset i
    i = 0;

    // char forcing[2][MAXSTRING];
    offsets[i] = offsetof(filenames_struct, forcing);
    blocklengths[i] *= 2;
    mpi_types[i++] = MPI_CHAR;

    // char f_path_pfx[2][MAXSTRING];
    offsets[i] = offsetof(filenames_struct, f_path_pfx);
    blocklengths[i] *= 2;
    mpi_types[i++] = MPI_CHAR;

    // char global[MAXSTRING];
    offsets[i] = offsetof(filenames_struct, global);
    mpi_types[i++] = MPI_CHAR;

    // char domain[MAXSTRING];
    offsets[i] = offsetof(filenames_struct, domain);
    mpi_types[i++] = MPI_CHAR;

    // char constants[MAXSTRING];
    offsets[i] = offsetof(filenames_struct, constants);
    mpi_types[i++] = MPI_CHAR;

    // char init_state[MAXSTRING];
    offsets[i] = offsetof(filenames_struct, init_state);
    mpi_types[i++] = MPI_CHAR;

    // char lakeparam[MAXSTRING];
    offsets[i] = offsetof(filenames_struct, lakeparam);
    mpi_types[i++] = MPI_CHAR;

    // char result_dir[MAXSTRING];
    offsets[i] = offsetof(filenames_struct, result_dir);
    mpi_types[i++] = MPI_CHAR;

    // char snowband[MAXSTRING];
    offsets[i] = offsetof(filenames_struct, snowband);
    mpi_types[i++] = MPI_CHAR;

    // char soil[MAXSTRING];
    offsets[i] = offsetof(filenames_struct, soil);
    mpi_types[i++] = MPI_CHAR;

    // char statefile[MAXSTRING];
    offsets[i] = offsetof(filenames_struct, statefile);
    mpi_types[i++] = MPI_CHAR;

    // char veg[MAXSTRING];
    offsets[i] = offsetof(filenames_struct, veg);
    mpi_types[i++] = MPI_CHAR;

    // char veglib[MAXSTRING];
    offsets[i] = offsetof(filenames_struct, veglib);
    mpi_types[i++] = MPI_CHAR;

    // char log_path[MAXSTRING];
    offsets[i] = offsetof(filenames_struct, log_path);
    mpi_types[i++] = MPI_CHAR;


    // make sure that the we have the right number of elements
    if (i != (size_t) nitems) {
        log_err("Miscount in create_MPI_filenames_struct_type(): "
                "%zd not equal to %d\n", i, nitems);
    }

    status = MPI_Type_create_struct(nitems, blocklengths, offsets, mpi_types,
                                    mpi_type);
    if (status != MPI_SUCCESS) {
        log_err("MPI error in create_MPI_filenames_struct_type(): %d\n",
                status);
    }

    status = MPI_Type_commit(mpi_type);
    if (status != MPI_SUCCESS) {
        log_err("MPI error in create_MPI_filenames_struct_type(): %d\n",
                status);
    }

    // cleanup
    free(blocklengths);
    free(offsets);
    free(mpi_types);
}

/******************************************************************************
 * @brief   Create an MPI_Datatype that represents the location_struct
 * @details This allows MPI operations in which the entire location_struct
 *          can be treated as an MPI_Datatype. NOTE: This function needs to be
 *          kept in-sync with the location_struct data type in
 *          vic_image_driver.h.
 *
 * @param mpi_type MPI_Datatype that can be used in MPI operations
 *****************************************************************************/
void
create_MPI_location_struct_type(MPI_Datatype *mpi_type)
{
    int           nitems; // number of elements in struct
    int           status;
    int          *blocklengths;
    size_t        i;
    MPI_Aint     *offsets;
    MPI_Datatype *mpi_types;

    // nitems has to equal the number of elements in global_param_struct
    nitems = 9;
    blocklengths = malloc(nitems * sizeof(*blocklengths));
    if (blocklengths == NULL) {
        log_err("Memory allocation error in create_MPI_location_struct_type().")
    }

    offsets = malloc(nitems * sizeof(*offsets));
    if (offsets == NULL) {
        log_err("Memory allocation error in create_MPI_location_struct_type().")
    }

    mpi_types = malloc(nitems * sizeof(*mpi_types));
    if (mpi_types == NULL) {
        log_err("Memory allocation error in create_MPI_location_struct_type().")
    }

    // none of the elements in location_struct are arrays.
    for (i = 0; i < (size_t) nitems; i++) {
        blocklengths[i] = 1;
    }

    // reset i
    i = 0;

    // size_t run;
    offsets[i] = offsetof(location_struct, run);
    mpi_types[i++] = MPI_C_BOOL;

    // double latitude;
    offsets[i] = offsetof(location_struct, latitude);
    mpi_types[i++] = MPI_DOUBLE;

    // double longitude;
    offsets[i] = offsetof(location_struct, longitude);
    mpi_types[i++] = MPI_DOUBLE;

    // double area;
    offsets[i] = offsetof(location_struct, area);
    mpi_types[i++] = MPI_DOUBLE;

    // double frac;
    offsets[i] = offsetof(location_struct, frac);
    mpi_types[i++] = MPI_DOUBLE;

    // size_t nveg;
    offsets[i] = offsetof(location_struct, nveg);
    mpi_types[i++] = MPI_AINT;

    // size_t global_idx;
    offsets[i] = offsetof(location_struct, global_idx);
    mpi_types[i++] = MPI_AINT;

    // size_t io_idx;
    offsets[i] = offsetof(location_struct, io_idx);
    mpi_types[i++] = MPI_AINT;

    // size_t local_idx;
    offsets[i] = offsetof(location_struct, local_idx);
    mpi_types[i++] = MPI_AINT;

    // make sure that the we have the right number of elements
    if (i != (size_t) nitems) {
        log_err("Miscount in create_MPI_location_struct_type(): "
                "%zd not equal to %d\n", i, nitems);
    }

    status = MPI_Type_create_struct(nitems, blocklengths, offsets, mpi_types,
                                    mpi_type);
    if (status != MPI_SUCCESS) {
        log_err("MPI error in create_MPI_location_struct_type(): %d\n", status);
    }

    status = MPI_Type_commit(mpi_type);
    if (status != MPI_SUCCESS) {
        log_err("MPI error in create_MPI_location_struct_type(): %d\n", status);
    }

    // cleanup
    free(blocklengths);
    free(offsets);
    free(mpi_types);
}

/******************************************************************************
 * @brief   Create an MPI_Datatype that represents the nc_file_struct
 * @details This allows MPI operations in which the entire nc_file_struct can
 *          be treated as an MPI_Datatype. NOTE: This function needs to be kept
 *          in-sync with the nc_file_struct data type in vic_driver_image.h.
 *
 * @param mpi_type MPI_Datatype that can be used in MPI operations
 *****************************************************************************/
void
create_MPI_nc_file_struct_type(MPI_Datatype *mpi_type)
{
    int           nitems; // number of elements in struct
    int           status;
    int          *blocklengths;
    size_t        i;
    MPI_Aint     *offsets;
    MPI_Datatype *mpi_types;

    // nitems has to equal the number of elements in nc_file_struct
    nitems = 27;
    blocklengths = malloc(nitems * sizeof(*blocklengths));
    if (blocklengths == NULL) {
        log_err("Memory allocation error in create_MPI_nc_file_struct_type().")
    }

    offsets = malloc(nitems * sizeof(*offsets));
    if (offsets == NULL) {
        log_err("Memory allocation error in create_MPI_nc_file_struct_type().")
    }

    mpi_types = malloc(nitems * sizeof(*mpi_types));
    if (mpi_types == NULL) {
        log_err("Memory allocation error in create_MPI_nc_file_struct_type().")
    }

    // only the first element in nc_file_struct is an array
    i = 0;
    blocklengths[i] = MAXSTRING + 1;
    for (i = 1; i < (size_t) nitems; i++) {
        blocklengths[i] = 1;
    }

    // reset i
    i = 0;

    // char fname[MAXSTRING + 1];
    offsets[i] = offsetof(nc_file_struct, fname);
    mpi_types[i++] = MPI_CHAR;

    // char c_fillvalue;
    offsets[i] = offsetof(nc_file_struct, c_fillvalue);
    mpi_types[i++] = MPI_CHAR;

    // int i_fillvalue;
    offsets[i] = offsetof(nc_file_struct, i_fillvalue);
    mpi_types[i++] = MPI_INT;

    // double d_fillvalue;
    offsets[i] = offsetof(nc_file_struct, d_fillvalue);
    mpi_types[i++] = MPI_DOUBLE;

    // float f_fillvalue;
    offsets[i] = offsetof(nc_file_struct, f_fillvalue);
    mpi_types[i++] = MPI_FLOAT;

    // int nc_id;
    offsets[i] = offsetof(nc_file_struct, nc_id);
    mpi_types[i++] = MPI_INT;

    // int band_dimid;
    offsets[i] = offsetof(nc_file_struct, band_dimid);
    mpi_types[i++] = MPI_INT;

    // int front_dimid;
    offsets[i] = offsetof(nc_file_struct, front_dimid);
    mpi_types[i++] = MPI_INT;

    // int frost_dimid;
    offsets[i] = offsetof(nc_file_struct, frost_dimid);
    mpi_types[i++] = MPI_INT;

    // int layer_dimid;
    offsets[i] = offsetof(nc_file_struct, layer_dimid);
    mpi_types[i++] = MPI_INT;

    // int ni_dimid;
    offsets[i] = offsetof(nc_file_struct, ni_dimid);
    mpi_types[i++] = MPI_INT;

    // int nj_dimid;
    offsets[i] = offsetof(nc_file_struct, nj_dimid);
    mpi_types[i++] = MPI_INT;

    // int node_dimid;
    offsets[i] = offsetof(nc_file_struct, node_dimid);
    mpi_types[i++] = MPI_INT;

    // int root_zone_dimid;
    offsets[i] = offsetof(nc_file_struct, root_zone_dimid);
    mpi_types[i++] = MPI_INT;

    // int time_dimid;
    offsets[i] = offsetof(nc_file_struct, time_dimid);
    mpi_types[i++] = MPI_INT;

    // int veg_dimid;
    offsets[i] = offsetof(nc_file_struct, veg_dimid);
    mpi_types[i++] = MPI_INT;

    // size_t band_size;
    offsets[i] = offsetof(nc_file_struct, band_size);
    mpi_types[i++] = MPI_AINT;

    // size_t front_size;
    offsets[i] = offsetof(nc_file_struct, front_size);
    mpi_types[i++] = MPI_AINT;

    // size_t frost_size;
    offsets[i] = offsetof(nc_file_struct, frost_size);
    mpi_types[i++] = MPI_AINT;

    // size_t layer_size;
    offsets[i] = offsetof(nc_file_struct, layer_size);
    mpi_types[i++] = MPI_AINT;

    // size_t ni_size;
    offsets[i] = offsetof(nc_file_struct, ni_size);
    mpi_types[i++] = MPI_AINT;

    // size_t nj_size;
    offsets[i] = offsetof(nc_file_struct, nj_size);
    mpi_types[i++] = MPI_AINT;

    // size_t node_size;
    offsets[i] = offsetof(nc_file_struct, node_size);
    mpi_types[i++] = MPI_AINT;

    // size_t root_zone_size;
    offsets[i] = offsetof(nc_file_struct, root_zone_size);
    mpi_types[i++] = MPI_AINT;

    // size_t time_size;
    offsets[i] = offsetof(nc_file_struct, time_size);
    mpi_types[i++] = MPI_AINT;

    // size_t veg_size;
    offsets[i] = offsetof(nc_file_struct, veg_size);
    mpi_types[i++] = MPI_AINT;

    // bool open;
    offsets[i] = offsetof(nc_file_struct, open);
    mpi_types[i++] = MPI_C_BOOL;

    // make sure that the we have the right number of elements
    if (i != (size_t) nitems) {
        log_err("Miscount in create_MPI_nc_file_struct_type(): "
                "%zd not equal to %d\n", i, nitems);
    }

    status = MPI_Type_create_struct(nitems, blocklengths, offsets, mpi_types,
                                    mpi_type);
    if (status != MPI_SUCCESS) {
        log_err("MPI error in create_MPI_nc_file_struct_type(): %d\n", status);
    }
    status = MPI_Type_commit(mpi_type);
    if (status != MPI_SUCCESS) {
        log_err("MPI error in create_MPI_nc_file_struct_type(): %d\n", status);
    }

    // cleanup
    free(blocklengths);
    free(offsets);
    free(mpi_types);
}

/******************************************************************************
 * @brief   Create an MPI_Datatype that represents the option_struct
 * @details This allows MPI operations in which the entire option_struct can
 *          be treated as an MPI_Datatype. NOTE: This function needs to be kept
 *          in-sync with the option_struct data type in vic_def.h.
 *
 * @param mpi_type MPI_Datatype that can be used in MPI operations
 *****************************************************************************/
void
create_MPI_option_struct_type(MPI_Datatype *mpi_type)
{
    int           nitems; // number of elements in struct
    int           status;
    int          *blocklengths;
    size_t        i;
    MPI_Aint     *offsets;
    MPI_Datatype *mpi_types;

    // nitems has to equal the number of elements in option_struct
    nitems = 60;
    blocklengths = malloc(nitems * sizeof(*blocklengths));
    if (blocklengths == NULL) {
        log_err("Memory allocation error in create_MPI_option_struct_type().")
    }

    offsets = malloc(nitems * sizeof(*offsets));
    if (offsets == NULL) {
        log_err("Memory allocation error in create_MPI_option_struct_type().")
    }

    mpi_types = malloc(nitems * sizeof(*mpi_types));
    if (mpi_types == NULL) {
        log_err("Memory allocation error in create_MPI_option_struct_type().")
    }

    // none of the elements in option_struct are arrays
    for (i = 0; i < (size_t) nitems; i++) {
        blocklengths[i] = 1;
    }

    // reset i
    i = 0;

    // short AboveTreelineVeg;
    offsets[i] = offsetof(option_struct, AboveTreelineVeg);
    mpi_types[i++] = MPI_SHORT;

    // unsigned short AERO_RESIST_CANSNOW;
    offsets[i] = offsetof(option_struct, AERO_RESIST_CANSNOW);
    mpi_types[i++] = MPI_UNSIGNED_SHORT;

    // bool BLOWING;
    offsets[i] = offsetof(option_struct, BLOWING);
    mpi_types[i++] = MPI_C_BOOL;

    // bool BLOWING_VAR_THRESHOLD;
    offsets[i] = offsetof(option_struct, BLOWING_VAR_THRESHOLD);
    mpi_types[i++] = MPI_C_BOOL;

    // bool BLOWING_CALC_PROB;
    offsets[i] = offsetof(option_struct, BLOWING_CALC_PROB);
    mpi_types[i++] = MPI_C_BOOL;

    // bool BLOWING_SIMPLE;
    offsets[i] = offsetof(option_struct, BLOWING_SIMPLE);
    mpi_types[i++] = MPI_C_BOOL;

    // bool BLOWING_FETCH;
    offsets[i] = offsetof(option_struct, BLOWING_FETCH);
    mpi_types[i++] = MPI_C_BOOL;

    // bool BLOWING_SPATIAL_WIND;
    offsets[i] = offsetof(option_struct, BLOWING_SPATIAL_WIND);
    mpi_types[i++] = MPI_C_BOOL;

    // bool CARBON;
    offsets[i] = offsetof(option_struct, CARBON);
    mpi_types[i++] = MPI_C_BOOL;

    // bool CLOSE_ENERGY;
    offsets[i] = offsetof(option_struct, CLOSE_ENERGY);
    mpi_types[i++] = MPI_C_BOOL;

    // bool COMPUTE_TREELINE;
    offsets[i] = offsetof(option_struct, COMPUTE_TREELINE);
    mpi_types[i++] = MPI_C_BOOL;

    // bool CONTINUEONERROR;
    offsets[i] = offsetof(option_struct, CONTINUEONERROR);
    mpi_types[i++] = MPI_C_BOOL;

    // bool CORRPREC;
    offsets[i] = offsetof(option_struct, CORRPREC);
    mpi_types[i++] = MPI_C_BOOL;

    // bool EQUAL_AREA;
    offsets[i] = offsetof(option_struct, EQUAL_AREA);
    mpi_types[i++] = MPI_C_BOOL;

    // bool EXP_TRANS;
    offsets[i] = offsetof(option_struct, EXP_TRANS);
    mpi_types[i++] = MPI_C_BOOL;

    // bool FROZEN_SOIL;
    offsets[i] = offsetof(option_struct, FROZEN_SOIL);
    mpi_types[i++] = MPI_C_BOOL;

    // bool FULL_ENERGY;
    offsets[i] = offsetof(option_struct, FULL_ENERGY);
    mpi_types[i++] = MPI_C_BOOL;

    // unsigned short GRND_FLUX_TYPE;
    offsets[i] = offsetof(option_struct, GRND_FLUX_TYPE);
    mpi_types[i++] = MPI_UNSIGNED_SHORT;

    // bool IMPLICIT;
    offsets[i] = offsetof(option_struct, IMPLICIT);
    mpi_types[i++] = MPI_C_BOOL;

    // bool JULY_TAVG_SUPPLIED;
    offsets[i] = offsetof(option_struct, JULY_TAVG_SUPPLIED);
    mpi_types[i++] = MPI_C_BOOL;

    // bool LAKES;
    offsets[i] = offsetof(option_struct, LAKES);
    mpi_types[i++] = MPI_C_BOOL;

    // size_t Ncanopy;
    offsets[i] = offsetof(option_struct, Ncanopy);
    mpi_types[i++] = MPI_AINT; // note there is no MPI_SIZE_T equivalent

    // size_t Nfrost;
    offsets[i] = offsetof(option_struct, Nfrost);
    mpi_types[i++] = MPI_AINT; // note there is no MPI_SIZE_T equivalent

    // size_t Nlakenode;
    offsets[i] = offsetof(option_struct, Nlakenode);
    mpi_types[i++] = MPI_AINT; // note there is no MPI_SIZE_T equivalent

    // size_t Nlayer;
    offsets[i] = offsetof(option_struct, Nlayer);
    mpi_types[i++] = MPI_AINT; // note there is no MPI_SIZE_T equivalent

    // size_t Nnode;
    offsets[i] = offsetof(option_struct, Nnode);
    mpi_types[i++] = MPI_AINT; // note there is no MPI_SIZE_T equivalent

    // bool NOFLUX;
    offsets[i] = offsetof(option_struct, NOFLUX);
    mpi_types[i++] = MPI_C_BOOL;

    // size_t NVEGTYPES;
    offsets[i] = offsetof(option_struct, NVEGTYPES);
    mpi_types[i++] = MPI_AINT; // note there is no MPI_SIZE_T equivalent

    // unsigned short RC_MODE;
    offsets[i] = offsetof(option_struct, RC_MODE);
    mpi_types[i++] = MPI_UNSIGNED_SHORT;

    // size_t ROOT_ZONES;
    offsets[i] = offsetof(option_struct, ROOT_ZONES);
    mpi_types[i++] = MPI_AINT; // note there is no MPI_SIZE_T equivalent

    // bool QUICK_FLUX;
    offsets[i] = offsetof(option_struct, QUICK_FLUX);
    mpi_types[i++] = MPI_C_BOOL;

    // bool QUICK_SOLVE;
    offsets[i] = offsetof(option_struct, QUICK_SOLVE);
    mpi_types[i++] = MPI_C_BOOL;

    // bool SHARE_LAYER_MOIST;
    offsets[i] = offsetof(option_struct, SHARE_LAYER_MOIST);
    mpi_types[i++] = MPI_C_BOOL;

    // unsigned short SNOW_DENSITY;
    offsets[i] = offsetof(option_struct, SNOW_DENSITY);
    mpi_types[i++] = MPI_UNSIGNED_SHORT;

    // size_t SNOW_BAND;
    offsets[i] = offsetof(option_struct, SNOW_BAND);
    mpi_types[i++] = MPI_AINT; // note there is no MPI_SIZE_T equivalent

    // bool SPATIAL_FROST;
    offsets[i] = offsetof(option_struct, SPATIAL_FROST);
    mpi_types[i++] = MPI_C_BOOL;

    // bool SPATIAL_SNOW;
    offsets[i] = offsetof(option_struct, SPATIAL_SNOW);
    mpi_types[i++] = MPI_C_BOOL;

    // bool TFALLBACK;
    offsets[i] = offsetof(option_struct, TFALLBACK);
    mpi_types[i++] = MPI_C_BOOL;

    // bool BASEFLOW;
    offsets[i] = offsetof(option_struct, BASEFLOW);
    mpi_types[i++] = MPI_C_BOOL;

    // unsigned short GRID_DECIMAL;
    offsets[i] = offsetof(option_struct, GRID_DECIMAL);
    mpi_types[i++] = MPI_UNSIGNED_SHORT;

    // bool VEGLIB_PHOTO;
    offsets[i] = offsetof(option_struct, VEGLIB_PHOTO);
    mpi_types[i++] = MPI_C_BOOL;

    // bool VEGLIB_FCAN;
    offsets[i] = offsetof(option_struct, VEGLIB_FCAN);
    mpi_types[i++] = MPI_C_BOOL;

    // bool VEGPARAM_ALB;
    offsets[i] = offsetof(option_struct, VEGPARAM_ALB);
    mpi_types[i++] = MPI_C_BOOL;

    // bool VEGPARAM_LAI;
    offsets[i] = offsetof(option_struct, VEGPARAM_LAI);
    mpi_types[i++] = MPI_C_BOOL;

    // bool VEGPARAM_FCAN;
    offsets[i] = offsetof(option_struct, VEGPARAM_FCAN);
    mpi_types[i++] = MPI_C_BOOL;

    // unsigned short ALB_SRC;
    offsets[i] = offsetof(option_struct, ALB_SRC);
    mpi_types[i++] = MPI_UNSIGNED_SHORT;

    // unsigned short LAI_SRC;
    offsets[i] = offsetof(option_struct, LAI_SRC);
    mpi_types[i++] = MPI_UNSIGNED_SHORT;

    // unsigned short FCAN_SRC;
    offsets[i] = offsetof(option_struct, FCAN_SRC);
    mpi_types[i++] = MPI_UNSIGNED_SHORT;

    // bool LAKE_PROFILE;
    offsets[i] = offsetof(option_struct, LAKE_PROFILE);
    mpi_types[i++] = MPI_C_BOOL;

    // bool ORGANIC_FRACT;
    offsets[i] = offsetof(option_struct, ORGANIC_FRACT);
    mpi_types[i++] = MPI_C_BOOL;

    // bool BINARY_STATE_FILE;
    offsets[i] = offsetof(option_struct, BINARY_STATE_FILE);
    mpi_types[i++] = MPI_C_BOOL;

    // bool INIT_STATE;
    offsets[i] = offsetof(option_struct, INIT_STATE);
    mpi_types[i++] = MPI_C_BOOL;

    // bool SAVE_STATE;
    offsets[i] = offsetof(option_struct, SAVE_STATE);
    mpi_types[i++] = MPI_C_BOOL;

    // bool ALMA_OUTPUT;
    offsets[i] = offsetof(option_struct, ALMA_OUTPUT);
    mpi_types[i++] = MPI_C_BOOL;

    // bool BINARY_OUTPUT;
    offsets[i] = offsetof(option_struct, BINARY_OUTPUT);
    mpi_types[i++] = MPI_C_BOOL;

    // bool COMPRESS;
    offsets[i] = offsetof(option_struct, COMPRESS);
    mpi_types[i++] = MPI_C_BOOL;

    // bool MOISTFRACT;
    offsets[i] = offsetof(option_struct, MOISTFRACT);
    mpi_types[i++] = MPI_C_BOOL;

    // size_t Noutfiles;
    offsets[i] = offsetof(option_struct, Noutfiles);
    mpi_types[i++] = MPI_AINT; // note there is no MPI_SIZE_T equivalent

    // bool PRT_HEADER;
    offsets[i] = offsetof(option_struct, PRT_HEADER);
    mpi_types[i++] = MPI_C_BOOL;

    // bool PRT_SNOW_BAND;
    offsets[i] = offsetof(option_struct, PRT_SNOW_BAND);
    mpi_types[i++] = MPI_C_BOOL;

    // make sure that the we have the right number of elements
    if (i != (size_t) nitems) {
        log_err("Miscount in create_MPI_option_struct_type(): "
                "%zd not equal to %d\n", i, nitems);
    }

    status = MPI_Type_create_struct(nitems, blocklengths, offsets, mpi_types,
                                    mpi_type);
    if (status != MPI_SUCCESS) {
        log_err("MPI error in create_MPI_option_struct_type(): %d\n", status);
    }
    status = MPI_Type_commit(mpi_type);
    if (status != MPI_SUCCESS) {
        log_err("MPI error in create_MPI_option_struct_type(): %d\n", status);
    }

    // cleanup
    free(blocklengths);
    free(offsets);
    free(mpi_types);
}

/******************************************************************************
 * @brief   Create an MPI_Datatype that represents the parameters_struct
 * @details This allows MPI operations in which the entire parameters_struct
 *          can be treated as an MPI_Datatype. NOTE: This function needs to be
 *          kept in-sync with the parameters_struct data type in vic_def.h.
 *
 * @param mpi_type MPI_Datatype that can be used in MPI operations
 *****************************************************************************/
void
create_MPI_param_struct_type(MPI_Datatype *mpi_type)
{
    int           nitems; // number of elements in struct
    int           status;
    int          *blocklengths;
    size_t        i;
    MPI_Aint     *offsets;
    MPI_Datatype *mpi_types;

    // nitems has to equal the number of elements in parameters_struct
    nitems = 139;
    blocklengths = malloc(nitems * sizeof(*blocklengths));
    if (blocklengths == NULL) {
        log_err("Memory allocation error in create_MPI_param_struct_type().")
    }

    offsets = malloc(nitems * sizeof(*offsets));
    if (offsets == NULL) {
        log_err("Memory allocation error in create_MPI_param_struct_type().")
    }

    mpi_types = malloc(nitems * sizeof(*mpi_types));
    if (mpi_types == NULL) {
        log_err("Memory allocation error in create_MPI_param_struct_type().")
    }

    // none of the elements in parameters_struct are arrays
    for (i = 0; i < (size_t) nitems; i++) {
        blocklengths[i] = 1;
    }

    // reset i
    i = 0;

    // double LAPSE_RATE;
    offsets[i] = offsetof(parameters_struct, LAPSE_RATE);
    mpi_types[i++] = MPI_DOUBLE;

    // double GAUGE_HEIGHT;
    offsets[i] = offsetof(parameters_struct, GAUGE_HEIGHT);
    mpi_types[i++] = MPI_DOUBLE;

    // double HUGE_RESIST;
    offsets[i] = offsetof(parameters_struct, HUGE_RESIST);
    mpi_types[i++] = MPI_DOUBLE;

    // double ALBEDO_BARE_SOIL;
    offsets[i] = offsetof(parameters_struct, ALBEDO_BARE_SOIL);
    mpi_types[i++] = MPI_DOUBLE;

    // double ALBEDO_H20_SURF;
    offsets[i] = offsetof(parameters_struct, ALBEDO_H20_SURF);
    mpi_types[i++] = MPI_DOUBLE;

    // double EMISS_GRND;
    offsets[i] = offsetof(parameters_struct, EMISS_GRND);
    mpi_types[i++] = MPI_DOUBLE;

    // double EMISS_VEG
    offsets[i] = offsetof(parameters_struct, EMISS_VEG);
    mpi_types[i++] = MPI_DOUBLE;

    // double EMISS_ICE
    offsets[i] = offsetof(parameters_struct, EMISS_ICE);
    mpi_types[i++] = MPI_DOUBLE;

    // double EMISS_SNOW
    offsets[i] = offsetof(parameters_struct, EMISS_SNOW);
    mpi_types[i++] = MPI_DOUBLE;

    // double EMISS_H2O
    offsets[i] = offsetof(parameters_struct, EMISS_H2O);
    mpi_types[i++] = MPI_DOUBLE;

    // double SOIL_RESID_MOIST
    offsets[i] = offsetof(parameters_struct, SOIL_RESID_MOIST);
    mpi_types[i++] = MPI_DOUBLE;

    // double SOIL_SLAB_MOIST_FRACT
    offsets[i] = offsetof(parameters_struct, SOIL_SLAB_MOIST_FRACT);
    mpi_types[i++] = MPI_DOUBLE;

    // double VEG_LAI_SNOW_MULTIPLIER
    offsets[i] = offsetof(parameters_struct, VEG_LAI_SNOW_MULTIPLIER);
    mpi_types[i++] = MPI_DOUBLE;

    // double VEG_MIN_INTERCEPTION_STORAGE
    offsets[i] = offsetof(parameters_struct, VEG_MIN_INTERCEPTION_STORAGE);
    mpi_types[i++] = MPI_DOUBLE;

    // double VEG_LAI_WATER_FACTOR
    offsets[i] = offsetof(parameters_struct, VEG_LAI_WATER_FACTOR);
    mpi_types[i++] = MPI_DOUBLE;

    // double CANOPY_CLOSURE
    offsets[i] = offsetof(parameters_struct, CANOPY_CLOSURE);
    mpi_types[i++] = MPI_DOUBLE;

    // double CANOPY_RSMAX
    offsets[i] = offsetof(parameters_struct, CANOPY_RSMAX);
    mpi_types[i++] = MPI_DOUBLE;

    // double CANOPY_VPDMINFACTOR
    offsets[i] = offsetof(parameters_struct, CANOPY_VPDMINFACTOR);
    mpi_types[i++] = MPI_DOUBLE;

    // double LAKE_TMELT
    offsets[i] = offsetof(parameters_struct, LAKE_TMELT);
    mpi_types[i++] = MPI_DOUBLE;

    // double LAKE_MAX_SURFACE
    offsets[i] = offsetof(parameters_struct, LAKE_MAX_SURFACE);
    mpi_types[i++] = MPI_DOUBLE;

    // double LAKE_BETA
    offsets[i] = offsetof(parameters_struct, LAKE_BETA);
    mpi_types[i++] = MPI_DOUBLE;

    // double LAKE_FRACMIN
    offsets[i] = offsetof(parameters_struct, LAKE_FRACMIN);
    mpi_types[i++] = MPI_DOUBLE;

    // double LAKE_FRACLIM
    offsets[i] = offsetof(parameters_struct, LAKE_FRACLIM);
    mpi_types[i++] = MPI_DOUBLE;

    // double LAKE_DM
    offsets[i] = offsetof(parameters_struct, LAKE_DM);
    mpi_types[i++] = MPI_DOUBLE;

    // double LAKE_SNOWCRIT
    offsets[i] = offsetof(parameters_struct, LAKE_SNOWCRIT);
    mpi_types[i++] = MPI_DOUBLE;

    // double LAKE_ZWATER
    offsets[i] = offsetof(parameters_struct, LAKE_ZWATER);
    mpi_types[i++] = MPI_DOUBLE;

    // double LAKE_ZSNOW
    offsets[i] = offsetof(parameters_struct, LAKE_ZSNOW);
    mpi_types[i++] = MPI_DOUBLE;

    // double LAKE_RHOSNOW
    offsets[i] = offsetof(parameters_struct, LAKE_RHOSNOW);
    mpi_types[i++] = MPI_DOUBLE;

    // double LAKE_CONDI
    offsets[i] = offsetof(parameters_struct, LAKE_CONDI);
    mpi_types[i++] = MPI_DOUBLE;

    // double LAKE_CONDS
    offsets[i] = offsetof(parameters_struct, LAKE_CONDS);
    mpi_types[i++] = MPI_DOUBLE;

    // double LAKE_LAMISW
    offsets[i] = offsetof(parameters_struct, LAKE_LAMISW);
    mpi_types[i++] = MPI_DOUBLE;

    // double LAKE_LAMILW
    offsets[i] = offsetof(parameters_struct, LAKE_LAMILW);
    mpi_types[i++] = MPI_DOUBLE;

    // double LAKE_LAMSSW
    offsets[i] = offsetof(parameters_struct, LAKE_LAMSSW);
    mpi_types[i++] = MPI_DOUBLE;

    // double LAKE_LAMSLW
    offsets[i] = offsetof(parameters_struct, LAKE_LAMSLW);
    mpi_types[i++] = MPI_DOUBLE;

    // double LAKE_LAMWSW
    offsets[i] = offsetof(parameters_struct, LAKE_LAMWSW);
    mpi_types[i++] = MPI_DOUBLE;

    // double LAKE_LAMWLW
    offsets[i] = offsetof(parameters_struct, LAKE_LAMWLW);
    mpi_types[i++] = MPI_DOUBLE;

    // double LAKE_A1
    offsets[i] = offsetof(parameters_struct, LAKE_A1);
    mpi_types[i++] = MPI_DOUBLE;

    // double LAKE_A2
    offsets[i] = offsetof(parameters_struct, LAKE_A2);
    mpi_types[i++] = MPI_DOUBLE;

    // double LAKE_QWTAU
    offsets[i] = offsetof(parameters_struct, LAKE_QWTAU);
    mpi_types[i++] = MPI_DOUBLE;

    // int LAKE_MAX_ITER
    offsets[i] = offsetof(parameters_struct, LAKE_MAX_ITER);
    mpi_types[i++] = MPI_INT;

    // double SVP_A
    offsets[i] = offsetof(parameters_struct, SVP_A);
    mpi_types[i++] = MPI_DOUBLE;

    // double SVP_B
    offsets[i] = offsetof(parameters_struct, SVP_B);
    mpi_types[i++] = MPI_DOUBLE;

    // double SVP_C
    offsets[i] = offsetof(parameters_struct, SVP_C);
    mpi_types[i++] = MPI_DOUBLE;

    // double PHOTO_OMEGA
    offsets[i] = offsetof(parameters_struct, PHOTO_OMEGA);
    mpi_types[i++] = MPI_DOUBLE;

    // double PHOTO_LAIMAX
    offsets[i] = offsetof(parameters_struct, PHOTO_LAIMAX);
    mpi_types[i++] = MPI_DOUBLE;

    // double PHOTO_LAILIMIT
    offsets[i] = offsetof(parameters_struct, PHOTO_LAILIMIT);
    mpi_types[i++] = MPI_DOUBLE;

    // double PHOTO_LAIMIN
    offsets[i] = offsetof(parameters_struct, PHOTO_LAIMIN);
    mpi_types[i++] = MPI_DOUBLE;

    // double PHOTO_EPAR
    offsets[i] = offsetof(parameters_struct, PHOTO_EPAR);
    mpi_types[i++] = MPI_DOUBLE;

    // double PHOTO_FCMAX
    offsets[i] = offsetof(parameters_struct, PHOTO_FCMAX);
    mpi_types[i++] = MPI_DOUBLE;

    // double PHOTO_FCMIN
    offsets[i] = offsetof(parameters_struct, PHOTO_FCMIN);
    mpi_types[i++] = MPI_DOUBLE;

    // double PHOTO_ZENITHMIN
    offsets[i] = offsetof(parameters_struct, PHOTO_ZENITHMIN);
    mpi_types[i++] = MPI_DOUBLE;

    // double PHOTO_ZENITHMINPAR
    offsets[i] = offsetof(parameters_struct, PHOTO_ZENITHMINPAR);
    mpi_types[i++] = MPI_DOUBLE;

    // double PHOTO_ALBSOIPARMIN
    offsets[i] = offsetof(parameters_struct, PHOTO_ALBSOIPARMIN);
    mpi_types[i++] = MPI_DOUBLE;

    // double PHOTO_MINMAXETRANS
    offsets[i] = offsetof(parameters_struct, PHOTO_MINMAXETRANS);
    mpi_types[i++] = MPI_DOUBLE;

    // double PHOTO_MINSTOMCOND
    offsets[i] = offsetof(parameters_struct, PHOTO_MINSTOMCOND);
    mpi_types[i++] = MPI_DOUBLE;

    // double PHOTO_FCI1C3
    offsets[i] = offsetof(parameters_struct, PHOTO_FCI1C3);
    mpi_types[i++] = MPI_DOUBLE;

    // double PHOTO_FCI1C4
    offsets[i] = offsetof(parameters_struct, PHOTO_FCI1C4);
    mpi_types[i++] = MPI_DOUBLE;

    // double PHOTO_OX
    offsets[i] = offsetof(parameters_struct, PHOTO_OX);
    mpi_types[i++] = MPI_DOUBLE;

    // double PHOTO_KCPHOTO_KC
    offsets[i] = offsetof(parameters_struct, PHOTO_KC);
    mpi_types[i++] = MPI_DOUBLE;

    // double PHOTO_KO
    offsets[i] = offsetof(parameters_struct, PHOTO_KO);
    mpi_types[i++] = MPI_DOUBLE;

    // double PHOTO_EC
    offsets[i] = offsetof(parameters_struct, PHOTO_EC);
    mpi_types[i++] = MPI_DOUBLE;

    // double PHOTO_EO
    offsets[i] = offsetof(parameters_struct, PHOTO_EO);
    mpi_types[i++] = MPI_DOUBLE;

    // double PHOTO_EV
    offsets[i] = offsetof(parameters_struct, PHOTO_EV);
    mpi_types[i++] = MPI_DOUBLE;

    // double PHOTO_ER
    offsets[i] = offsetof(parameters_struct, PHOTO_ER);
    mpi_types[i++] = MPI_DOUBLE;

    // double PHOTO_ALC3
    offsets[i] = offsetof(parameters_struct, PHOTO_ALC3);
    mpi_types[i++] = MPI_DOUBLE;

    // double PHOTO_FRDC3
    offsets[i] = offsetof(parameters_struct, PHOTO_FRDC3);
    mpi_types[i++] = MPI_DOUBLE;

    // double PHOTO_EK
    offsets[i] = offsetof(parameters_struct, PHOTO_EK);
    mpi_types[i++] = MPI_DOUBLE;

    // double PHOTO_ALC4
    offsets[i] = offsetof(parameters_struct, PHOTO_ALC4);
    mpi_types[i++] = MPI_DOUBLE;

    // double PHOTO_FRDC4
    offsets[i] = offsetof(parameters_struct, PHOTO_FRDC4);
    mpi_types[i++] = MPI_DOUBLE;

    // double PHOTO_THETA
    offsets[i] = offsetof(parameters_struct, PHOTO_THETA);
    mpi_types[i++] = MPI_DOUBLE;

    // double PHOTO_FRLEAF
    offsets[i] = offsetof(parameters_struct, PHOTO_FRLEAF);
    mpi_types[i++] = MPI_DOUBLE;

    // double PHOTO_FRGROWTH
    offsets[i] = offsetof(parameters_struct, PHOTO_FRGROWTH);
    mpi_types[i++] = MPI_DOUBLE;

    // double SRESP_E0_LT
    offsets[i] = offsetof(parameters_struct, SRESP_E0_LT);
    mpi_types[i++] = MPI_DOUBLE;

    // double SRESP_T0_LT
    offsets[i] = offsetof(parameters_struct, SRESP_T0_LT);
    mpi_types[i++] = MPI_DOUBLE;

    // double SRESP_WMINFM
    offsets[i] = offsetof(parameters_struct, SRESP_WMINFM);
    mpi_types[i++] = MPI_DOUBLE;

    // double SRESP_WMAXFM
    offsets[i] = offsetof(parameters_struct, SRESP_WMAXFM);
    mpi_types[i++] = MPI_DOUBLE;

    // double SRESP_WOPTFM
    offsets[i] = offsetof(parameters_struct, SRESP_WOPTFM);
    mpi_types[i++] = MPI_DOUBLE;

    // double SRESP_RHSAT
    offsets[i] = offsetof(parameters_struct, SRESP_RHSAT);
    mpi_types[i++] = MPI_DOUBLE;

    // double SRESP_RFACTOR
    offsets[i] = offsetof(parameters_struct, SRESP_RFACTOR);
    mpi_types[i++] = MPI_DOUBLE;

    // double SRESP_TAULITTER
    offsets[i] = offsetof(parameters_struct, SRESP_TAULITTER);
    mpi_types[i++] = MPI_DOUBLE;

    // double SRESP_TAUINTER
    offsets[i] = offsetof(parameters_struct, SRESP_TAUINTER);
    mpi_types[i++] = MPI_DOUBLE;

    // double SRESP_TAUSLOW
    offsets[i] = offsetof(parameters_struct, SRESP_TAUSLOW);
    mpi_types[i++] = MPI_DOUBLE;

    // double SRESP_FAIR
    offsets[i] = offsetof(parameters_struct, SRESP_FAIR);
    mpi_types[i++] = MPI_DOUBLE;

    // double SRESP_FINTER
    offsets[i] = offsetof(parameters_struct, SRESP_FINTER);
    mpi_types[i++] = MPI_DOUBLE;

    // double SNOW_MAX_SURFACE_SWE
    offsets[i] = offsetof(parameters_struct, SNOW_MAX_SURFACE_SWE);
    mpi_types[i++] = MPI_DOUBLE;

    // double SNOW_LIQUID_WATER_CAPACITY
    offsets[i] = offsetof(parameters_struct, SNOW_LIQUID_WATER_CAPACITY);
    mpi_types[i++] = MPI_DOUBLE;

    // double SNOW_NEW_SNOW_DENSITY
    offsets[i] = offsetof(parameters_struct, SNOW_NEW_SNOW_DENSITY);
    mpi_types[i++] = MPI_DOUBLE;

    // double SNOW_DENS_DMLIMIT
    offsets[i] = offsetof(parameters_struct, SNOW_DENS_DMLIMIT);
    mpi_types[i++] = MPI_DOUBLE;

    // double SNOW_DENS_MAX_CHANGE
    offsets[i] = offsetof(parameters_struct, SNOW_DENS_MAX_CHANGE);
    mpi_types[i++] = MPI_DOUBLE;

    // double SNOW_DENS_ETA0
    offsets[i] = offsetof(parameters_struct, SNOW_DENS_ETA0);
    mpi_types[i++] = MPI_DOUBLE;

    // double SNOW_DENS_C1
    offsets[i] = offsetof(parameters_struct, SNOW_DENS_C1);
    mpi_types[i++] = MPI_DOUBLE;

    // double SNOW_DENS_C2
    offsets[i] = offsetof(parameters_struct, SNOW_DENS_C2);
    mpi_types[i++] = MPI_DOUBLE;

    // double SNOW_DENS_C5
    offsets[i] = offsetof(parameters_struct, SNOW_DENS_C5);
    mpi_types[i++] = MPI_DOUBLE;

    // double SNOW_DENS_C6
    offsets[i] = offsetof(parameters_struct, SNOW_DENS_C6);
    mpi_types[i++] = MPI_DOUBLE;

    // double SNOW_DENS_F
    offsets[i] = offsetof(parameters_struct, SNOW_DENS_F);
    mpi_types[i++] = MPI_DOUBLE;

    // double SNOW_MIN_SWQ_EB_THRES
    offsets[i] = offsetof(parameters_struct, SNOW_MIN_SWQ_EB_THRES);
    mpi_types[i++] = MPI_DOUBLE;

    // double SNOW_A1
    offsets[i] = offsetof(parameters_struct, SNOW_A1);
    mpi_types[i++] = MPI_DOUBLE;

    // double SNOW_A2
    offsets[i] = offsetof(parameters_struct, SNOW_A2);
    mpi_types[i++] = MPI_DOUBLE;

    // double SNOW_L1
    offsets[i] = offsetof(parameters_struct, SNOW_L1);
    mpi_types[i++] = MPI_DOUBLE;

    // double SNOW_L2
    offsets[i] = offsetof(parameters_struct, SNOW_L2);
    mpi_types[i++] = MPI_DOUBLE;

    // double SNOW_NEW_SNOW_ALB
    offsets[i] = offsetof(parameters_struct, SNOW_NEW_SNOW_ALB);
    mpi_types[i++] = MPI_DOUBLE;

    // double SNOW_ALB_ACCUM_A
    offsets[i] = offsetof(parameters_struct, SNOW_ALB_ACCUM_A);
    mpi_types[i++] = MPI_DOUBLE;

    // double SNOW_ALB_ACCUM_B
    offsets[i] = offsetof(parameters_struct, SNOW_ALB_ACCUM_B);
    mpi_types[i++] = MPI_DOUBLE;

    // double SNOW_ALB_THAW_A
    offsets[i] = offsetof(parameters_struct, SNOW_ALB_THAW_A);
    mpi_types[i++] = MPI_DOUBLE;

    // double SNOW_ALB_THAW_B
    offsets[i] = offsetof(parameters_struct, SNOW_ALB_THAW_B);
    mpi_types[i++] = MPI_DOUBLE;

    // double SNOW_TRACESNOW
    offsets[i] = offsetof(parameters_struct, SNOW_TRACESNOW);
    mpi_types[i++] = MPI_DOUBLE;

    // double SNOW_CONDUCT
    offsets[i] = offsetof(parameters_struct, SNOW_CONDUCT);
    mpi_types[i++] = MPI_DOUBLE;

    // double SNOW_MAX_SNOW_TEMP
    offsets[i] = offsetof(parameters_struct, SNOW_MAX_SNOW_TEMP);
    mpi_types[i++] = MPI_DOUBLE;

    // double SNOW_MIN_RAIN_TEMP
    offsets[i] = offsetof(parameters_struct, SNOW_MIN_RAIN_TEMP);
    mpi_types[i++] = MPI_DOUBLE;

    // double BLOWING_KA
    offsets[i] = offsetof(parameters_struct, BLOWING_KA);
    mpi_types[i++] = MPI_DOUBLE;

    // double BLOWING_CSALT
    offsets[i] = offsetof(parameters_struct, BLOWING_CSALT);
    mpi_types[i++] = MPI_DOUBLE;

    // double BLOWING_UTHRESH
    offsets[i] = offsetof(parameters_struct, BLOWING_UTHRESH);
    mpi_types[i++] = MPI_DOUBLE;

    // double BLOWING_KIN_VIS
    offsets[i] = offsetof(parameters_struct, BLOWING_KIN_VIS);
    mpi_types[i++] = MPI_DOUBLE;

    // int BLOWING_MAX_ITER
    offsets[i] = offsetof(parameters_struct, BLOWING_MAX_ITER);
    mpi_types[i++] = MPI_INT;

    // int BLOWING_K
    offsets[i] = offsetof(parameters_struct, BLOWING_K);
    mpi_types[i++] = MPI_INT;

    // double BLOWING_SETTLING
    offsets[i] = offsetof(parameters_struct, BLOWING_SETTLING);
    mpi_types[i++] = MPI_DOUBLE;

    // int BLOWING_NUMINCS
    offsets[i] = offsetof(parameters_struct, BLOWING_NUMINCS);
    mpi_types[i++] = MPI_INT;

    // double TREELINE_TEMPERATURE
    offsets[i] = offsetof(parameters_struct, TREELINE_TEMPERATURE);
    mpi_types[i++] = MPI_DOUBLE;

    // double SNOW_DT
    offsets[i] = offsetof(parameters_struct, SNOW_DT);
    mpi_types[i++] = MPI_DOUBLE;

    // double SURF_DT
    offsets[i] = offsetof(parameters_struct, SURF_DT);
    mpi_types[i++] = MPI_DOUBLE;

    // double SOIL_DT
    offsets[i] = offsetof(parameters_struct, SOIL_DT);
    mpi_types[i++] = MPI_DOUBLE;

    // double CANOPY_DT
    offsets[i] = offsetof(parameters_struct, CANOPY_DT);
    mpi_types[i++] = MPI_DOUBLE;

    // double CANOPY_VP
    offsets[i] = offsetof(parameters_struct, CANOPY_VP);
    mpi_types[i++] = MPI_DOUBLE;

    // double TOL_GRND
    offsets[i] = offsetof(parameters_struct, TOL_GRND);
    mpi_types[i++] = MPI_DOUBLE;

    // double TOL_OVER
    offsets[i] = offsetof(parameters_struct, TOL_OVER);
    mpi_types[i++] = MPI_DOUBLE;

    // int FROZEN_MAXITER
    offsets[i] = offsetof(parameters_struct, FROZEN_MAXITER);
    mpi_types[i++] = MPI_INT;

    // int NEWT_RAPH_MAXTRIAL
    offsets[i] = offsetof(parameters_struct, NEWT_RAPH_MAXTRIAL);
    mpi_types[i++] = MPI_INT;

    // double NEWT_RAPH_TOLX
    offsets[i] = offsetof(parameters_struct, NEWT_RAPH_TOLX);
    mpi_types[i++] = MPI_DOUBLE;

    // double NEWT_RAPH_TOLF
    offsets[i] = offsetof(parameters_struct, NEWT_RAPH_TOLF);
    mpi_types[i++] = MPI_DOUBLE;

    // double NEWT_RAPH_R_MAX
    offsets[i] = offsetof(parameters_struct, NEWT_RAPH_R_MAX);
    mpi_types[i++] = MPI_DOUBLE;

    // double NEWT_RAPH_R_MIN
    offsets[i] = offsetof(parameters_struct, NEWT_RAPH_R_MIN);
    mpi_types[i++] = MPI_DOUBLE;

    // double NEWT_RAPH_RELAX1
    offsets[i] = offsetof(parameters_struct, NEWT_RAPH_RELAX1);
    mpi_types[i++] = MPI_DOUBLE;

    // double NEWT_RAPH_RELAX2
    offsets[i] = offsetof(parameters_struct, NEWT_RAPH_RELAX2);
    mpi_types[i++] = MPI_DOUBLE;

    // double NEWT_RAPH_RELAX3
    offsets[i] = offsetof(parameters_struct, NEWT_RAPH_RELAX3);
    mpi_types[i++] = MPI_DOUBLE;

    // double NEWT_RAPH_EPS2
    offsets[i] = offsetof(parameters_struct, NEWT_RAPH_EPS2);
    mpi_types[i++] = MPI_DOUBLE;

    // int ROOT_BRENT_MAXTRIES
    offsets[i] = offsetof(parameters_struct, ROOT_BRENT_MAXTRIES);
    mpi_types[i++] = MPI_INT;

    // int ROOT_BRENT_MAXITER
    offsets[i] = offsetof(parameters_struct, ROOT_BRENT_MAXITER);
    mpi_types[i++] = MPI_INT;

    // double ROOT_BRENT_TSTEP
    offsets[i] = offsetof(parameters_struct, ROOT_BRENT_TSTEP);
    mpi_types[i++] = MPI_DOUBLE;

    // double ROOT_BRENT_T
    offsets[i] = offsetof(parameters_struct, ROOT_BRENT_T);
    mpi_types[i++] = MPI_DOUBLE;

    // make sure that the we have the right number of elements
    if (i != (size_t) nitems) {
        log_err("Miscount in create_MPI_param_struct_type(): "
                "%zd not equal to %d\n", i, nitems);
    }

    status = MPI_Type_create_struct(nitems, blocklengths, offsets, mpi_types,
                                    mpi_type);
    if (status != MPI_SUCCESS) {
        log_err("MPI error in create_MPI_param_struct_type(): %d\n", status);
    }
    status = MPI_Type_commit(mpi_type);
    if (status != MPI_SUCCESS) {
        log_err("MPI error in create_MPI_param_struct_type(): %d\n", status);
    }

    // cleanup
    free(blocklengths);
    free(offsets);
    free(mpi_types);
}

/******************************************************************************
 * @brief   Type-agnostic mapping function
 * @details Reorders the elements in 'from' to 'to' according to the ordering
 *          specified in 'map'.
 *          Note that this function can also be used for filtering, i.e. you
 *          can use a smaller number of elements in 'map' and 'to' than in
 *          'from' to get only a subset of the elements.
 *
 *          to[to_map[i]] = from[from_map[i]]
 *
 * @param size size of the datatype of 'from' and 'to', e.g. sizeof(int)
 * @param n number of elements in 'from_map' and 'to_map'
 * @param from_map array of length n with 'from' indices, if from_map == NULL,
 *        then the 'from' indices are sequential
 * @param to_map array of length n with 'to' indices, if to_map == NULL, then
 *        the 'to' indices are sequential
 * @param from array of with entries of size 'size' (unchanged)
 * @param to array of with entries of size 'size' (changed)
 *****************************************************************************/
void
map(size_t  size,
    size_t  n,
    size_t *from_map,
    size_t *to_map,
    void   *from,
    void   *to)
{
    size_t i;

    if (to_map == NULL && from_map == NULL) {
        for (i = 0; i < n; i++) {
            // type-agnostic version of to[i] = from[i];
            memcpy(to + i * size, from + i * size, size);
        }
    }
    if (to_map == NULL) {
        for (i = 0; i < n; i++) {
            // type-agnostic version of to[i] = from[from_map[i]];
            memcpy(to + i * size, from + from_map[i] * size, size);
        }
    }
    else if (from_map == NULL) {
        for (i = 0; i < n; i++) {
            // type-agnostic version of to[to_map[i]] = from[i];
            memcpy(to + to_map[i] * size, from + i * size, size);
        }
    }
    else {
        for (i = 0; i < n; i++) {
            // type-agnostic version of to[to_map[i]] = from[from_map[i]];
            memcpy(to + to_map[i] * size, from + from_map[i] * size, size);
        }
    }
}

/******************************************************************************
 * @brief   Decompose the domain for MPI operations
 * @details This function sets up the arrays needed to scatter and gather
 *          data from and to the master process to the individual mpi
 *          processes.
 *
 * @param ncells total number of cells
 * @param mpi_size number of mpi processes
 * @param mpi_map_local_array_sizes address of integer array with number of
 *        cells assigned to each node (MPI_Scatterv:sendcounts and
 *        MPI_Gatherv:recvcounts)
 * @param mpi_map_global_array_offsets address of integer array with offsets
 *        for sending and receiving data (MPI_Scatterv:displs and
 *        MPI_Gatherv:displs)
 * @param mpi_map_mapping_array address of size_t array with indices to prepare
 *        an array on the master process for MPI_Scatterv or map back after
 *        MPI_Gatherv
 *****************************************************************************/
void
mpi_map_decomp_domain(size_t   ncells,
                      size_t   mpi_size,
                      int    **mpi_map_local_array_sizes,
                      int    **mpi_map_global_array_offsets,
                      size_t **mpi_map_mapping_array)
{
    size_t i;
    size_t j;
    size_t k;
    size_t n;

    *mpi_map_local_array_sizes = calloc(mpi_size,
                                        sizeof(*(*mpi_map_local_array_sizes)));
    *mpi_map_global_array_offsets = calloc(mpi_size,
                                           sizeof(*(*
                                                    mpi_map_global_array_offsets)));
    *mpi_map_mapping_array = calloc(ncells, sizeof(*(*mpi_map_mapping_array)));

    // determine number of cells per node
    for (n = ncells, i = 0; n > 0; n--, i++) {
        if (i >= mpi_size) {
            i = 0;
        }
        (*mpi_map_local_array_sizes)[i] += 1;
    }

    // determine offsets to use for MPI_Scatterv and MPI_Gatherv
    for (i = 1; i < mpi_size; i++) {
        for (j = 0; j < i; j++) {
            (*mpi_map_global_array_offsets)[i] +=
                (*mpi_map_local_array_sizes)[j];
        }
    }

    // set mapping array
    for (i = 0, k = 0; i < (size_t) mpi_size; i++) {
        for (j = 0; j < (size_t) (*mpi_map_local_array_sizes)[i]; j++) {
            (*mpi_map_mapping_array)[k++] = (size_t) (i + j * mpi_size);
        }
    }
}

/******************************************************************************
 * @brief   Gather and write double precision NetCDF field
 * @details Values are gathered to the master node and then written from the
 *          master node
 *****************************************************************************/
void
gather_put_nc_field_double(char   *nc_name,
                           bool   *open,
                           int    *nc_id,
                           double  fillval,
                           int    *dimids,
                           int     ndims,
                           char   *var_name,
                           size_t *start,
                           size_t *count,
                           double *var)
{
    int                  status;
    double              *dvar = NULL;
    double              *dvar_gathered = NULL;
    double              *dvar_remapped = NULL;
    size_t               grid_size;
    size_t               i;
    extern domain_struct global_domain;
    extern domain_struct local_domain;
    extern int           mpi_rank;
    extern int          *mpi_map_global_array_offsets;
    extern int          *mpi_map_local_array_sizes;
    extern size_t       *filter_active_cells;
    extern size_t       *mpi_map_mapping_array;

    if (mpi_rank == 0) {
        grid_size = global_domain.n_nx * global_domain.n_ny;
        dvar = malloc(grid_size * sizeof(*dvar));
        if (dvar == NULL) {
            log_err("Memory allocation error in gather_put_nc_field_double().");
        }
        for (i = 0; i < grid_size; i++) {
            dvar[i] = fillval;
        }

        dvar_gathered =
            malloc(global_domain.ncells_active * sizeof(*dvar_gathered));
        if (dvar_gathered == NULL) {
            log_err("Memory allocation error in gather_put_nc_field_double().");
        }

        dvar_remapped =
            malloc(global_domain.ncells_active * sizeof(*dvar_remapped));
        if (dvar_remapped == NULL) {
            log_err("Memory allocation error in gather_put_nc_field_double().");
        }
    }
    // Gather the results from the nodes, result for the local node is in the
    // array *var (which is a function argument)
    status = MPI_Gatherv(var, local_domain.ncells_active, MPI_DOUBLE,
                         dvar_gathered, mpi_map_local_array_sizes,
                         mpi_map_global_array_offsets, MPI_DOUBLE,
                         0, MPI_COMM_WORLD);
    if (status != MPI_SUCCESS) {
        fprintf(stderr, "MPI error in main(): %d\n", status);
        exit(EXIT_FAILURE);
    }
    if (mpi_rank == 0) {
        // remap the array
        map(sizeof(double), global_domain.ncells_active, NULL,
            mpi_map_mapping_array,
            dvar_gathered, dvar_remapped);
        // expand to full grid size
        map(sizeof(double), global_domain.ncells_active, NULL,
            filter_active_cells,
            dvar_remapped, dvar);
        // write to file
        put_nc_field_double(nc_name, open, nc_id, fillval, dimids, ndims,
                            var_name, start, count, dvar);
        // cleanup
        free(dvar);
        free(dvar_gathered);
        free(dvar_remapped);
    }
}

/******************************************************************************
 * @brief   Gather and write integer NetCDF field
 * @details Values are gathered to the master node and then written from the
 *          master node
 *****************************************************************************/
void
gather_put_nc_field_int(char   *nc_name,
                        bool   *open,
                        int    *nc_id,
                        int     fillval,
                        int    *dimids,
                        int     ndims,
                        char   *var_name,
                        size_t *start,
                        size_t *count,
                        int    *var)
{
    int                  status;
    int                 *ivar = NULL;
    int                 *ivar_gathered = NULL;
    int                 *ivar_remapped = NULL;
    size_t               grid_size;
    size_t               i;
    extern domain_struct global_domain;
    extern domain_struct local_domain;
    extern int           mpi_rank;
    extern int          *mpi_map_global_array_offsets;
    extern int          *mpi_map_local_array_sizes;
    extern size_t       *filter_active_cells;
    extern size_t       *mpi_map_mapping_array;

    if (mpi_rank == 0) {
        grid_size = global_domain.n_nx * global_domain.n_ny;
        ivar = malloc(grid_size * sizeof(*ivar));
        if (ivar == NULL) {
            log_err("Memory allocation error in gather_put_nc_field_int().");
        }
        for (i = 0; i < grid_size; i++) {
            ivar[i] = fillval;
        }

        ivar_gathered =
            malloc(global_domain.ncells_active * sizeof(*ivar_gathered));
        if (ivar_gathered == NULL) {
            log_err("Memory allocation error in gather_put_nc_field_int().");
        }

        ivar_remapped =
            malloc(global_domain.ncells_active * sizeof(*ivar_remapped));
        if (ivar_remapped == NULL) {
            log_err("Memory allocation error in gather_put_nc_field_int().");
        }
    }
    // Gather the results from the nodes, result for the local node is in the
    // array *var (which is a function argument)
    status = MPI_Gatherv(var, local_domain.ncells_active, MPI_INT,
                         ivar_gathered, mpi_map_local_array_sizes,
                         mpi_map_global_array_offsets, MPI_INT,
                         0, MPI_COMM_WORLD);
    if (status != MPI_SUCCESS) {
        fprintf(stderr, "MPI error in main(): %d\n", status);
        exit(EXIT_FAILURE);
    }

    if (mpi_rank == 0) {
        // remap the array
        map(sizeof(int), global_domain.ncells_active, NULL,
            mpi_map_mapping_array,
            ivar_gathered, ivar_remapped);
        // expand to full grid size
        map(sizeof(int), global_domain.ncells_active, NULL, filter_active_cells,
            ivar_remapped, ivar);
        // write to file
        put_nc_field_int(nc_name, open, nc_id, fillval, dimids, ndims,
                         var_name, start, count, ivar);
        // cleanup
        free(ivar);
        free(ivar_gathered);
        free(ivar_remapped);
    }
}

/******************************************************************************
 * @brief   Read double precision NetCDF field from file and scatter
 * @details Read happens on the master node and is then scattered to the local
 *          nodes
 *****************************************************************************/
void
get_scatter_nc_field_double(char   *nc_name,
                            char   *var_name,
                            size_t *start,
                            size_t *count,
                            double *var)
{
    int                  status;
    double              *dvar = NULL;
    double              *dvar_filtered = NULL;
    double              *dvar_mapped = NULL;
    extern domain_struct global_domain;
    extern domain_struct local_domain;
    extern int           mpi_rank;
    extern int          *mpi_map_global_array_offsets;
    extern int          *mpi_map_local_array_sizes;
    extern size_t       *filter_active_cells;
    extern size_t       *mpi_map_mapping_array;

    if (mpi_rank == 0) {
        dvar = malloc(global_domain.ncells_total * sizeof(*dvar));
        if (dvar == NULL) {
            log_err("Memory allocation error in get_scatter_nc_field_double().");
        }
        dvar_filtered =
            malloc(global_domain.ncells_active * sizeof(*dvar_filtered));
        if (dvar_filtered == NULL) {
            log_err("Memory allocation error in get_scatter_nc_field_double().");
        }
        dvar_mapped =
            malloc(global_domain.ncells_active * sizeof(*dvar_mapped));
        if (dvar_mapped == NULL) {
            log_err("Memory allocation error in get_scatter_nc_field_double().");
        }
        get_nc_field_double(nc_name, var_name, start, count, dvar);
        // filter the active cells only
        map(sizeof(double), global_domain.ncells_active, filter_active_cells,
            NULL,
            dvar, dvar_filtered);
        // map to prepare for MPI_Scatterv
        map(sizeof(double), global_domain.ncells_active, mpi_map_mapping_array,
            NULL,
            dvar_filtered, dvar_mapped);
        free(dvar);
        free(dvar_filtered);
    }

    // Scatter the results to the nodes, result for the local node is in the
    // array *var (which is a function argument)
    status = MPI_Scatterv(dvar_mapped, mpi_map_local_array_sizes,
                          mpi_map_global_array_offsets, MPI_DOUBLE,
                          var, local_domain.ncells_active, MPI_DOUBLE,
                          0, MPI_COMM_WORLD);
    if (status != MPI_SUCCESS) {
        fprintf(stderr, "MPI error in main(): %d\n", status);
        exit(EXIT_FAILURE);
    }

    if (mpi_rank == 0) {
        free(dvar_mapped);
    }
}

/******************************************************************************
 * @brief   Read single precision NetCDF field from file and scatter
 * @details Read happens on the master node and is then scattered to the local
 *          nodes
 *****************************************************************************/
void
get_scatter_nc_field_float(char   *nc_name,
                           char   *var_name,
                           size_t *start,
                           size_t *count,
                           float  *var)
{
    int                  status;
    float               *fvar = NULL;
    float               *fvar_filtered = NULL;
    float               *fvar_mapped = NULL;
    extern domain_struct global_domain;
    extern domain_struct local_domain;
    extern int           mpi_rank;
    extern int          *mpi_map_global_array_offsets;
    extern int          *mpi_map_local_array_sizes;
    extern size_t       *filter_active_cells;
    extern size_t       *mpi_map_mapping_array;

    if (mpi_rank == 0) {
        fvar = malloc(global_domain.ncells_total * sizeof(*fvar));
        if (fvar == NULL) {
            log_err("Memory allocation error in get_scatter_nc_field_float().");
        }
        fvar_filtered =
            malloc(global_domain.ncells_active * sizeof(*fvar_filtered));
        if (fvar_filtered == NULL) {
            log_err("Memory allocation error in get_scatter_nc_field_float().");
        }
        fvar_mapped =
            malloc(global_domain.ncells_active * sizeof(*fvar_mapped));
        if (fvar_mapped == NULL) {
            log_err("Memory allocation error in get_scatter_nc_field_float().");
        }
        get_nc_field_float(nc_name, var_name, start, count, fvar);
        // filter the active cells only
        map(sizeof(float), global_domain.ncells_active, filter_active_cells,
            NULL,
            fvar, fvar_filtered);
        // map to prepare for MPI_Scatterv
        map(sizeof(float), global_domain.ncells_active, mpi_map_mapping_array,
            NULL,
            fvar_filtered, fvar_mapped);
        free(fvar);
        free(fvar_filtered);
    }

    // Scatter the results to the nodes, result for the local node is in the
    // array *var (which is a function argument)
    status = MPI_Scatterv(fvar_mapped, mpi_map_local_array_sizes,
                          mpi_map_global_array_offsets, MPI_FLOAT,
                          var, local_domain.ncells_active, MPI_FLOAT,
                          0, MPI_COMM_WORLD);
    if (status != MPI_SUCCESS) {
        fprintf(stderr, "MPI error in main(): %d\n", status);
        exit(EXIT_FAILURE);
    }

    if (mpi_rank == 0) {
        free(fvar_mapped);
    }
}

/******************************************************************************
 * @brief   Read integer NetCDF field from file and scatter
 * @details Read happens on the master node and is then scattered to the local
 *          nodes
 *****************************************************************************/
void
get_scatter_nc_field_int(char   *nc_name,
                         char   *var_name,
                         size_t *start,
                         size_t *count,
                         int    *var)
{
    int                  status;
    int                 *ivar = NULL;
    int                 *ivar_filtered = NULL;
    int                 *ivar_mapped = NULL;
    extern domain_struct global_domain;
    extern domain_struct local_domain;
    extern int           mpi_rank;
    extern int          *mpi_map_global_array_offsets;
    extern int          *mpi_map_local_array_sizes;
    extern size_t       *filter_active_cells;
    extern size_t       *mpi_map_mapping_array;

    if (mpi_rank == 0) {
        ivar = malloc(global_domain.ncells_total * sizeof(*ivar));
        if (ivar == NULL) {
            log_err("Memory allocation error in get_scatter_nc_field_int().");
        }
        ivar_filtered =
            malloc(global_domain.ncells_active * sizeof(*ivar_filtered));
        if (ivar_filtered == NULL) {
            log_err("Memory allocation error in get_scatter_nc_field_int().");
        }
        ivar_mapped =
            malloc(global_domain.ncells_active * sizeof(*ivar_mapped));
        if (ivar_mapped == NULL) {
            log_err("Memory allocation error in get_scatter_nc_field_int().");
        }
        get_nc_field_int(nc_name, var_name, start, count, ivar);
        // filter the active cells only
        map(sizeof(int), global_domain.ncells_active, filter_active_cells, NULL,
            ivar, ivar_filtered);
        // map to prepare for MPI_Scatterv
        map(sizeof(int), global_domain.ncells_active, mpi_map_mapping_array,
            NULL,
            ivar_filtered, ivar_mapped);
        free(ivar);
        free(ivar_filtered);
    }

    // Scatter the results to the nodes, result for the local node is in the
    // array *var (which is a function argument)
    status = MPI_Scatterv(ivar_mapped, mpi_map_local_array_sizes,
                          mpi_map_global_array_offsets, MPI_INT,
                          var, local_domain.ncells_active, MPI_INT,
                          0, MPI_COMM_WORLD);
    if (status != MPI_SUCCESS) {
        fprintf(stderr, "MPI error in main(): %d\n", status);
        exit(EXIT_FAILURE);
    }

    if (mpi_rank == 0) {
        free(ivar_mapped);
    }
}

#ifdef VIC_MPI_SUPPORT_TEST

#include <vic_driver_shared.h>
#include <assert.h>

// size_t              NF, NR;
// size_t              current;
size_t *filter_active_cells = NULL;
size_t *mpi_map_mapping_array = NULL;
// all_vars_struct    *all_vars = NULL;
// atmos_data_struct  *atmos = NULL;
// dmy_struct         *dmy = NULL;
filenames_struct filenames;
filep_struct     filep;
domain_struct    global_domain;
domain_struct    local_domain;
// global_param_struct global_param;
// lake_con_struct     lake_con;
MPI_Datatype     mpi_global_struct_type;
MPI_Datatype     mpi_location_struct_type;
MPI_Datatype     mpi_nc_file_struct_type;
MPI_Datatype     mpi_option_struct_type;
MPI_Datatype     mpi_param_struct_type;
int             *mpi_map_local_array_sizes = NULL;
int             *mpi_map_global_array_offsets = NULL;
int              mpi_rank;
int              mpi_size;
// nc_file_struct      nc_hist_file;
// nc_var_struct       nc_vars[N_OUTVAR_TYPES];
// option_struct       options;
// parameters_struct   param;
// out_data_struct   **out_data;
// save_data_struct   *save_data;
// param_set_struct    param_set;
// soil_con_struct    *soil_con = NULL;
// veg_con_map_struct *veg_con_map = NULL;
// veg_con_struct    **veg_con = NULL;
// veg_hist_struct   **veg_hist = NULL;
// veg_lib_struct    **veg_lib = NULL;
FILE *LOG_DEST;

/******************************************************************************
 * @brief   Test routine for VIC MPI support functions
 * @details Note: need to define VIC_MPI_SUPPORT_TEST to compile.
 *          For example (from drivers/image/src):
 *          mpicc-mpich-mp -Wall -Wextra -o vic_mpi_support vic_mpi_support.c
 *          get_nc_field.c put_nc_field.c ../../shared/src/vic_log.c
 *          ../../shared/src/open_file.c -I ../include/ -I ../../shared/include
 *          -I ../../../vic_run/include/ -I/opt/local/include
 *          -L/opt/local/lib -lnetcdf -lmpi -DVIC_MPI_SUPPORT_TEST
 *****************************************************************************/
int
main(int    argc,
     char **argv)
{
    int                 status;
    global_param_struct global;
    location_struct     location;
    nc_file_struct      ncfile;
    option_struct       option;
    parameters_struct   param;

    // Initialize Log Destination
    strcpy(filenames.log_path, "MISSING");
    initialize_log();

    status = MPI_Init(&argc, &argv);
    if (status != MPI_SUCCESS) {
        log_err("MPI error in main(): %d\n", status);
    }
    status = MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    if (status != MPI_SUCCESS) {
        log_err("MPI error in main(): %d\n", status);
    }
    status = MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    if (status != MPI_SUCCESS) {
        log_err("MPI error in main(): %d\n", status);
    }

    create_MPI_global_struct_type(&mpi_global_struct_type);
    create_MPI_location_struct_type(&mpi_location_struct_type);
    create_MPI_nc_file_struct_type(&mpi_nc_file_struct_type);
    create_MPI_option_struct_type(&mpi_option_struct_type);
    create_MPI_param_struct_type(&mpi_param_struct_type);

    if (mpi_rank == 0) {
        // populate the test structure on the master
        // make sure to test the last element of the structure, because any
        // problem with alignment would show there
        global.forceskip[0] = 4321;
        global.forceskip[1] = 8765;
        // last element of global
        global.time_origin_num = -12345.6789;

        location.latitude = 47.620607;
        location.longitude = -122.349281;
        // last element of location
        location.local_idx = 12345678;

        sprintf(ncfile.fname, "Space Needle, Seattle, WA");
        // last element of ncfile
        ncfile.open = true;

        option.CARBON = false;
        option.CLOSE_ENERGY = true;
        option.COMPUTE_TREELINE = true;
        option.CONTINUEONERROR = true;
        option.Nfrost = 1234;
        option.Nlakenode = 5678;
        option.PRT_HEADER = false;
        // last element of option
        option.PRT_SNOW_BAND = true;

        param.ROOT_BRENT_MAXTRIES = 6543;
        param.ROOT_BRENT_MAXITER = 1010;
        param.TOL_GRND = 0.001;
        // last element of param
        param.ROOT_BRENT_T = -98765432.;
    }

    // broadcast to the slaves
    status = MPI_Bcast(&global, 1, mpi_global_struct_type,
                       0, MPI_COMM_WORLD);
    if (status != MPI_SUCCESS) {
        log_err("MPI error in main(): %d\n", status);
    }

    status = MPI_Bcast(&location, 1, mpi_location_struct_type,
                       0, MPI_COMM_WORLD);
    if (status != MPI_SUCCESS) {
        log_err("MPI error in main(): %d\n", status);
    }

    status = MPI_Bcast(&ncfile, 1, mpi_nc_file_struct_type,
                       0, MPI_COMM_WORLD);
    if (status != MPI_SUCCESS) {
        log_err("MPI error in main(): %d\n", status);
    }

    status = MPI_Bcast(&option, 1, mpi_option_struct_type,
                       0, MPI_COMM_WORLD);
    if (status != MPI_SUCCESS) {
        log_err("MPI error in main(): %d\n", status);
    }

    status = MPI_Bcast(&param, 1, mpi_param_struct_type,
                       0, MPI_COMM_WORLD);
    if (status != MPI_SUCCESS) {
        log_err("MPI error in main(): %d\n", status);
    }

    // assert the values on all processes are the same as the original
    // assignment
    printf("%d: global.forceskip == %d\n", mpi_rank, global.forceskip[0]);
    assert(global.forceskip[0] == 4321);
    printf("%d: global.forceskip == %d\n", mpi_rank, global.forceskip[1]);
    assert(global.forceskip[1] == 8765);
    printf("%d: global.time_origin_num == %f\n", mpi_rank,
           global.time_origin_num);
    assert(global.time_origin_num == -12345.6789);
    printf("%d: location.latitude == %f\n", mpi_rank, location.latitude);
    assert(location.latitude == 47.620607);
    printf("%d: location.longitude == %f\n", mpi_rank, location.longitude);
    assert(location.longitude == -122.349281);
    printf("%d: location.local_idx == %zd\n", mpi_rank, location.local_idx);
    assert(location.local_idx == 12345678);
    printf("%d: ncfile.fname == %s\n", mpi_rank, ncfile.fname);
    assert(strcmp(ncfile.fname, "Space Needle, Seattle, WA") == 0);
    printf("%d: ncfile.open == %d\n", mpi_rank, ncfile.open);
    assert(ncfile.open == true);
    printf("%d: option.CARBON == %d\n", mpi_rank, option.CARBON);
    assert(option.CARBON == false);
    printf("%d: option.CLOSE_ENERGY == %d\n", mpi_rank, option.CLOSE_ENERGY);
    assert(option.CLOSE_ENERGY == true);
    printf("%d: option.COMPUTE_TREELINE == %d\n", mpi_rank,
           option.COMPUTE_TREELINE);
    assert(option.COMPUTE_TREELINE == true);
    printf("%d: option.CONTINUEONERROR == %d\n", mpi_rank,
           option.CONTINUEONERROR);
    assert(option.CONTINUEONERROR == true);
    printf("%d: option.Nfrost == %zd\n", mpi_rank, option.Nfrost);
    assert(option.Nfrost == 1234);
    printf("%d: option.Nlakenode == %zd\n", mpi_rank, option.Nlakenode);
    assert(option.Nlakenode == 5678);
    printf("%d: option.PRT_HEADER == %d\n", mpi_rank, option.PRT_HEADER);
    assert(option.PRT_HEADER == false);
    printf("%d: option.PRT_SNOW_BAND == %d\n", mpi_rank, option.PRT_SNOW_BAND);
    assert(option.PRT_SNOW_BAND == true);
    printf("%d: param.ROOT_BRENT_MAXTRIES == %d\n",
           mpi_rank, param.ROOT_BRENT_MAXTRIES);
    assert(param.ROOT_BRENT_MAXTRIES == 6543);
    printf("%d: param.ROOT_BRENT_MAXITER == %d\n",
           mpi_rank, param.ROOT_BRENT_MAXITER);
    assert(param.ROOT_BRENT_MAXITER == 1010);
    printf("%d: param.TOL_GRND == %f\n", mpi_rank, param.TOL_GRND);
    assert(param.TOL_GRND == 0.001);
    printf("%d: param.ROOT_BRENT_T == %f\n", mpi_rank, param.ROOT_BRENT_T);
    assert(param.ROOT_BRENT_T == -98765432.);

    status = MPI_Finalize();
    if (status != MPI_SUCCESS) {
        log_err("MPI error in main(): %d\n", status);
    }

    return EXIT_SUCCESS;
}

#endif
