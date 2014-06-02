#include <vic_def.h>
#include <vic_run.h>
#include <vic_driver_image.h>

#define ERR(e) {fprintf(stderr, "\nError(put_nc_field): %s\n", \
                        nc_strerror(e)); }

// TBD: note that in this implementation the assumption is that all NetCDF files
// are organized in the same way. This is a bad assumption and should be made
// more robust by actually checking the mapping from the lat lon to the cell
// indices

int
put_nc_field_double(char   *nc_name,
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
    int old_fill_mode;
    int status;
    int var_id;

    if (!(*open)) {
        // open the netcdf file
        status = nc_open(nc_name, 
                         NC_WRITE|NC_NOCLOBBER|NC_NETCDF4|NC_CLASSIC_MODEL,
                         nc_id);
        if (status != NC_NOERR) {
            ERR(status);
        } 
        *open = true;

        // set the NC_FILL attribute
        status = nc_set_fill(*nc_id, NC_FILL, &old_fill_mode);
            if (status != NC_NOERR) {
            ERR(status);
        }
    }

    /* get NetCDF variable */
    status = nc_inq_varid(*nc_id, var_name, &var_id);
    if (status == NC_ENOTVAR) {
        // enter define mode
        status = nc_redef(*nc_id);
        if (status != NC_NOERR) {
            ERR(status);
        }
        // define the variable
        status = nc_def_var(*nc_id, var_name, NC_DOUBLE, ndims, dimids, 
                            &var_id);
        if (status != NC_NOERR) {
            ERR(status);
        }
        // set the fill value attribute
        status = nc_put_att_double(*nc_id, var_id, "_FillValue", NC_DOUBLE, 1,
                                   &fillval);
        if (status != NC_NOERR) {
            fprintf(stderr, "nc_att_double %s\n", var_name);
            ERR(status);
        }
        // leave define mode        
        status = nc_enddef(*nc_id);
        if (status != NC_NOERR) {
            ERR(status);
        }
    } else if (status != NC_NOERR) {
        ERR(status);
    }

    status = nc_put_vara_double(*nc_id, var_id, start, count, var);
    if (status != NC_NOERR) {
        ERR(status);
    }

    // keep the file open
    
    return status;
}
