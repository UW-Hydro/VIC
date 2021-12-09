/******************************************************************************
 * @section DESCRIPTION
 *
 * Get netCDF variable type.
 *****************************************************************************/

#include <vic_driver_shared_image.h>

/******************************************************************************
 * @brief    Get netCDF variable type.
 *****************************************************************************/
int
get_nc_var_type(nameid_struct *nc_nameid,
                char          *var_name)
{
    int var_id;
    int status;
    int xtypep;

    // get variable id
    status = nc_inq_varid(nc_nameid->nc_id, var_name, &var_id);
    check_nc_status(status, "Error getting variable id %s in %s", var_name,
                    nc_nameid->nc_filename);

    // get type ID
    status = nc_inq_var(nc_nameid->nc_id, var_id, NULL, &xtypep, NULL, NULL,
                        NULL);
    check_nc_status(status, "Error getting variable type %s in %s", var_name,
                    nc_nameid->nc_filename);

    return(xtypep);
}
