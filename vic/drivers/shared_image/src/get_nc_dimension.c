/******************************************************************************
 * @section DESCRIPTION
 *
 * Get netCDF dimension.
 *****************************************************************************/

#include <vic_driver_shared_image.h>

/******************************************************************************
 * @brief    Get netCDF dimension.
 *****************************************************************************/
size_t
get_nc_dimension(nameid_struct *nc_nameid,
                 char          *dim_name)
{
    int    dim_id;
    size_t dim_size;
    int    status;

    // get dimension id
    status = nc_inq_dimid(nc_nameid->nc_id, dim_name, &dim_id);
    check_nc_status(status, "Error getting dimension id %s in %s", dim_name,
                    nc_nameid->nc_filename);

    // get dimension size
    status = nc_inq_dimlen(nc_nameid->nc_id, dim_id, &dim_size);
    check_nc_status(status, "Error getting dimension size for dim %s in %s",
                    dim_name,
                    nc_nameid->nc_filename);

    return dim_size;
}
