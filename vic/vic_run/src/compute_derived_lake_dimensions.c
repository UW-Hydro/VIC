/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine computes those lake variables that are completely dependent
 * on lake depth and basin dimensions.
 *****************************************************************************/

#include <vic_driver_shared_all.h>

/******************************************************************************
 * @brief    This routine computes those lake variables that are completely
 *           dependent on lake depth and basin dimensions.
 *****************************************************************************/
void
compute_derived_lake_dimensions(lake_var_struct *lake,
                                lake_con_struct  lake_con)
{
    extern parameters_struct param;

    int                      k;
    int                      status;
    double                   depth;
    double                   tmp_volume;

    /* number and thicknesses of lake layers */
    if (lake->ldepth > param.LAKE_MAX_SURFACE && lake->ldepth < 2 *
        param.LAKE_MAX_SURFACE) {
        /* Not quite enough for two full layers. */
        lake->surfdz = lake->ldepth / 2.;
        lake->dz = lake->ldepth / 2.;
        lake->activenod = 2;
    }
    else if (lake->ldepth >= 2 * param.LAKE_MAX_SURFACE) {
        /* More than two layers. */
        lake->surfdz = param.LAKE_MAX_SURFACE;
        lake->activenod = (int) (lake->ldepth / param.LAKE_MAX_SURFACE);
        if (lake->activenod > MAX_LAKE_NODES) {
            lake->activenod = MAX_LAKE_NODES;
        }
        lake->dz = (lake->ldepth - lake->surfdz) /
                   ((double) (lake->activenod - 1));
    }
    else if (lake->ldepth > DBL_EPSILON) {
        lake->surfdz = lake->ldepth;
        lake->dz = 0.0;
        lake->activenod = 1;
    }
    else {
        lake->surfdz = 0.0;
        lake->dz = 0.0;
        lake->activenod = 0;
        lake->ldepth = 0.0;
    }

    // lake_con.basin equals the surface area at specific depths as input by
    // the user in the lake parameter file or calculated in read_lakeparam(),
    // lake->surface equals the area at the top of each dynamic solution layer

    for (k = 0; k <= lake->activenod; k++) {
        if (k == 0) {
            depth = lake->ldepth;
        }
        else {
            depth = lake->dz * (lake->activenod - k);
        }
        status = get_sarea(lake_con, depth, &(lake->surface[k]));
        if (status < 0) {
            log_err("record = %d, depth = %f, "
                    "sarea = %e", 0, depth, lake->surface[k]);
        }
    }

    lake->sarea = lake->surface[0];
    status = get_volume(lake_con, lake->ldepth, &tmp_volume);
    if (status < 0) {
        log_err("record = %d, depth = %f, "
                "volume = %e", 0, depth, tmp_volume);
    }
    else if (status > 0) {
        log_err("lake depth exceeds maximum; "
                "setting to maximum; record = %d", 0);
    }
    lake->volume = tmp_volume + lake->ice_water_eq;
}
