/******************************************************************************
* @section DESCRIPTION
*
* This subroutine computes dependent lake parameters from the specified
* lake parameters.
******************************************************************************/

#include <vic_driver_shared_all.h>

/******************************************************************************
* @brief        This subroutine computes dependent lake parameters from the
*               specified lake parameters.
******************************************************************************/
void
compute_lake_params(lake_con_struct *lake_con,
                    soil_con_struct  soil_con)
{
    extern parameters_struct param;
    extern option_struct     options;

    size_t                   i;
    double                   tempdz;
    double                   radius;
    double                   x;
    int                      ErrFlag;

    // miscellaneous lake parameters
    lake_con->bpercent = lake_con->rpercent;
    lake_con->maxdepth = lake_con->z[0];
    lake_con->basin[0] = lake_con->Cl[0] * soil_con.cell_area;

    if (!options.LAKE_PROFILE) {
        // generate lake depth-area relationship
        tempdz = (lake_con->maxdepth) / ((double) lake_con->numnod);
        radius = sqrt(lake_con->basin[0] / CONST_PI);

        for (i = 1; i <= lake_con->numnod; i++) {
            lake_con->z[i] = (lake_con->numnod - i) * tempdz;
            if (lake_con->z[i] < 0.0) {
                lake_con->z[i] = 0.0;
            }
            x =
                pow(lake_con->z[i] / lake_con->maxdepth,
                    param.LAKE_BETA) * radius;
            lake_con->basin[i] = CONST_PI * x * x;
            lake_con->Cl[i] = lake_con->basin[i] / soil_con.cell_area;
        }
    }
    else {
        // final point in depth-area relationship
        lake_con->z[lake_con->numnod] = 0;
        lake_con->Cl[lake_con->numnod] = 0;

        // depth-area relationship specified (for area fractions)
        // compute basin node surface areas
        for (i = 1; i <= lake_con->numnod; i++) {
            lake_con->basin[i] = lake_con->Cl[i] * soil_con.cell_area;
        }
    }

    // compute max volume
    lake_con->maxvolume = 0.0;
    for (i = 1; i <= lake_con->numnod; i++) {
        lake_con->maxvolume += (lake_con->basin[i] + lake_con->basin[i - 1]) *
                               (lake_con->z[i - 1] - lake_con->z[i]) / 2.;
    }

    // compute volume corresponding to mindepth
    ErrFlag = get_volume(*lake_con, lake_con->mindepth, &(lake_con->minvolume));
    if (ErrFlag == ERROR) {
        log_err("Error calculating depth: depth %f volume %f",
                lake_con->mindepth, lake_con->minvolume);
    }
}
