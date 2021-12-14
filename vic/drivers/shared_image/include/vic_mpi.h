/******************************************************************************
 * @section DESCRIPTION
 *
 * Header file for vic_mpi routines
 *****************************************************************************/

#ifndef VIC_MPI_H
#define VIC_MPI_H

#include <vic_def.h>
#include <mpi.h>
#ifdef _OPENMP
    #include <omp.h>
#else
    #define omp_get_max_threads() 1
#endif

#define VIC_MPI_ROOT 0

/******************************************************************************
 * @brief   This structure stores netCDF file name and corresponding nc_id
 *****************************************************************************/
typedef struct {
    char nc_filename[MAXSTRING];
    int nc_id;
} nameid_struct;

void create_MPI_filenames_struct_type(MPI_Datatype *mpi_type);
void create_MPI_global_struct_type(MPI_Datatype *mpi_type);
void create_MPI_location_struct_type(MPI_Datatype *mpi_type);
void create_MPI_alarm_struct_type(MPI_Datatype *mpi_type);
void create_MPI_option_struct_type(MPI_Datatype *mpi_type);
void create_MPI_param_struct_type(MPI_Datatype *mpi_type);
void gather_field_double(double fillval, double *dvar, double *var);
void gather_put_nc_field_double(int nc_id, int var_id, double fillval,
                                size_t *start, size_t *count, double *var);
void gather_put_nc_field_float(int nc_id, int var_id, float fillval,
                               size_t *start, size_t *count, float *var);
void gather_put_nc_field_int(int nc_id, int var_id, int fillval, size_t *start,
                             size_t *count, int *var);
void gather_put_nc_field_short(int nc_id, int var_id, short int fillval,
                               size_t *start, size_t *count, short int *var);
void gather_put_nc_field_schar(int nc_id, int var_id, char fillval,
                               size_t *start, size_t *count, char *var);
void scatter_field_double(double *dvar, double *var);
void get_scatter_nc_field_double(nameid_struct *nc_nameid, char *var_name,
                                 size_t *start, size_t *count, double *var);
void get_scatter_nc_field_float(nameid_struct *nc_nameid, char *var_name,
                                size_t *start, size_t *count, float *var);
void get_scatter_nc_field_int(nameid_struct *nc_nameid, char *var_name,
                              size_t *start, size_t *count, int *var);
void initialize_mpi(void);
void map(size_t size, size_t n, size_t *from_map, size_t *to_map, void *from,
         void *to);
void mpi_map_decomp_domain(size_t ncells, size_t mpi_size,
                           int **mpi_map_local_array_sizes,
                           int **mpi_map_global_array_offsets,
                           size_t **mpi_map_mapping_array);
void print_mpi_error_str(int error_code);

#endif
