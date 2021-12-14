/******************************************************************************
 * @section DESCRIPTION
 *
 * Get netCDF dimension.
 *****************************************************************************/

#include <vic_driver_shared_image.h>

/******************************************************************************
 * @brief    Get netCDF dimension.
 *****************************************************************************/
int
get_nc_varndimensions(nameid_struct *nc_nameid,
                      char          *var_name)
{
    int var_id;
    int ndims;
    int status;

    // get variable id
    status = nc_inq_varid(nc_nameid->nc_id, var_name, &var_id);
    check_nc_status(status, "Error getting variable id %s in %s", var_name,
                    nc_nameid->nc_filename);

    // get number of dimensions
    status = nc_inq_varndims(nc_nameid->nc_id, var_id, &ndims);
    check_nc_status(status,
                    "Error getting number of dimensions for var %s in %s",
                    var_name,
                    nc_nameid->nc_filename);

    return ndims;
}
