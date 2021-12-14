/******************************************************************************
 * @section DESCRIPTION
 *
 * Functions to support reading from netCDF field.
 *****************************************************************************/

#include <vic_driver_shared_image.h>

/******************************************************************************
 * @brief    Read double precision netCDF field from file.
 *****************************************************************************/
int
get_nc_field_double(nameid_struct *nc_nameid,
                    char          *var_name,
                    size_t        *start,
                    size_t        *count,
                    double        *var)
{
    int status;
    int var_id;

    /* get NetCDF variable */
    status = nc_inq_varid(nc_nameid->nc_id, var_name, &var_id);
    check_nc_status(status, "Error getting variable id for %s in %s", var_name,
                    nc_nameid->nc_filename);

    status = nc_get_vara_double(nc_nameid->nc_id, var_id, start, count, var);
    check_nc_status(status, "Error getting values for %s in %s", var_name,
                    nc_nameid->nc_filename);

    return status;
}

/******************************************************************************
 * @brief    Read single precision netCDF field from file.
 *****************************************************************************/
int
get_nc_field_float(nameid_struct *nc_nameid,
                   char          *var_name,
                   size_t        *start,
                   size_t        *count,
                   float         *var)
{
    int status;
    int var_id;

    /* get NetCDF variable */
    status = nc_inq_varid(nc_nameid->nc_id, var_name, &var_id);
    check_nc_status(status, "Error getting variable id for %s in %s", var_name,
                    nc_nameid->nc_filename);

    status = nc_get_vara_float(nc_nameid->nc_id, var_id, start, count, var);
    check_nc_status(status, "Error getting values for %s in %s", var_name,
                    nc_nameid->nc_filename);

    return status;
}

/******************************************************************************
 * @brief    Read integer netCDF field from file.
 *****************************************************************************/
int
get_nc_field_int(nameid_struct *nc_nameid,
                 char          *var_name,
                 size_t        *start,
                 size_t        *count,
                 int           *var)
{
    int status;
    int var_id;

    /* get NetCDF variable */
    status = nc_inq_varid(nc_nameid->nc_id, var_name, &var_id);
    check_nc_status(status, "Error getting variable id for %s in %s", var_name,
                    nc_nameid->nc_filename);

    status = nc_get_vara_int(nc_nameid->nc_id, var_id, start, count, var);
    check_nc_status(status, "Error getting values for %s in %s", var_name,
                    nc_nameid->nc_filename);

    return status;
}
