#include <vic_def.h>
#include <vic_run.h>
#include <vic_driver_image.h>

#define ERR(e) {fprintf(stderr, "\nError(get_nc_field): %s\n", \
                        nc_strerror(e)); }

// TBD: note that in this implementation the assumption is that all NetCDF files
// are organized in the same way. This is a bad assumption and should be made
// more robust by actually checking the mapping from the lat lon to the cell
// indices

int
get_nc_field_double(char   *nc_name,
                    char   *var_name,
                    size_t *start,
                    size_t *count,
                    double *var)
{
    int    nc_id;
    int    status;
    int    var_id;

    size_t i;

    // open the netcdf file
    status = nc_open(nc_name, NC_NOWRITE, &nc_id);
    if (status != NC_NOERR) {
        ERR(status);
    }

    /* get NetCDF variable */
    status = nc_inq_varid(nc_id, var_name, &var_id);
    if (status != NC_NOERR) {
        ERR(status);
    }

    status = nc_get_vara_double(nc_id, var_id, start, count, var);
    if (status != NC_NOERR) {
        ERR(status);
    }

    // close the netcdf file
    status = nc_close(nc_id);
    if (status != NC_NOERR) {
        ERR(status);
    }

    return status;
}

int
get_nc_field_int(char   *nc_name,
                 char   *var_name,
                 size_t *start,
                 size_t *count,
                 int    *var)
{
    int    nc_id;
    int    status;
    int    var_id;

    size_t i;

    // open the netcdf file
    status = nc_open(nc_name, NC_NOWRITE, &nc_id);
    if (status != NC_NOERR) {
        ERR(status);
    }

    /* get NetCDF variable */
    status = nc_inq_varid(nc_id, var_name, &var_id);
    if (status != NC_NOERR) {
        ERR(status);
    }

    status = nc_get_vara_int(nc_id, var_id, start, count, var);
    if (status != NC_NOERR) {
        ERR(status);
    }

    // close the netcdf file
    status = nc_close(nc_id);
    if (status != NC_NOERR) {
        ERR(status);
    }

    return status;
}
