#include <vic_def.h>
#include <vic_run.h>
#include <vic_driver_image.h>

#define ERR(e) {fprintf(stderr, "\nError(get_nc_dimension): %s\n", \
                        nc_strerror(e)); }

size_t
get_nc_dimension(char *nc_name,
                 char *dim_name)
{
    int    nc_id;
    int    dim_id;
    size_t dim_size;
    int    status;

    printf("%s: %s\n", nc_name, dim_name);
    fflush(stdout);

    // open the netcdf file
    status = nc_open(nc_name, NC_NOWRITE, &nc_id);
    if (status != NC_NOERR) {
        ERR(status);
    }

    // get dimension id
    status = nc_inq_dimid(nc_id, dim_name, &dim_id);
    if (status != NC_NOERR) {
        ERR(status);
    }

    // get dimension size
    status = nc_inq_dimlen(nc_id, dim_id, &dim_size);
    if (status != NC_NOERR) {
        ERR(status);
    }

    // close the netcdf file
    status = nc_close(nc_id);
    if (status != NC_NOERR) {
        ERR(status);
    }

    return dim_size;
}
