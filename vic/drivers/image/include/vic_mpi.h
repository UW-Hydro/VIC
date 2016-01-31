/******************************************************************************
 * @section DESCRIPTION
 *
 * Header file for vic_mpi routines
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

#ifndef VIC_MPI_H
#define VIC_MPI_H

#include <stddef.h> // used for offsetof() macro
#include <stdbool.h>  // used for bool types
#include <mpi.h>

void create_MPI_filenames_struct_type(MPI_Datatype *mpi_type);
void create_MPI_global_struct_type(MPI_Datatype *mpi_type);
void create_MPI_location_struct_type(MPI_Datatype *mpi_type);
void create_MPI_nc_file_struct_type(MPI_Datatype *mpi_type);
void create_MPI_option_struct_type(MPI_Datatype *mpi_type);
void create_MPI_param_struct_type(MPI_Datatype *mpi_type);
void gather_put_nc_field_double(char *nc_name, bool *open, int *nc_id,
                                double fillval, int *dimids, int ndims,
                                char *var_name, size_t *start, size_t *count,
                                double *var);
void gather_put_nc_field_int(char *nc_name, bool *open, int *nc_id, int fillval,
                             int *dimids, int ndims, char *var_name,
                             size_t *start, size_t *count, int *var);
void get_scatter_nc_field_double(char *nc_name, char *var_name, size_t *start,
                                 size_t *count, double *var);
void get_scatter_nc_field_float(char *nc_name, char *var_name, size_t *start,
                                size_t *count, float *var);
void get_scatter_nc_field_int(char *nc_name, char *var_name, size_t *start,
                              size_t *count, int *var);
void initialize_mpi(void);
void map(size_t size, size_t n, size_t *from_map, size_t *to_map, void *from,
         void *to);
void mpi_map_decomp_domain(size_t ncells, size_t mpi_size,
                           int **mpi_map_local_array_sizes,
                           int **mpi_map_global_array_offsets,
                           size_t **mpi_map_mapping_array);

#endif
