/******************************************************************************
 * @section DESCRIPTION
 *
 * Put attribute to netCDF file
 *****************************************************************************/

#include <vic_driver_shared_image.h>

/******************************************************************************
 * @brief    Put text attribute to netCDF.
 *****************************************************************************/
void
put_nc_attr(int         nc_id,
            int         var_id,
            const char *name,
            const char *value)
{
    int status;

    status = nc_put_att_text(nc_id, var_id, name, strlen(value), value);
    check_nc_status(status, "Error adding %s attribute in ncid %d", name,
                    nc_id);
}
