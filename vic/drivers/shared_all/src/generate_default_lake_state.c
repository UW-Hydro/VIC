/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine initializes the lake model state (energy balance, water balance,
 * snow components) to default values.
 *****************************************************************************/

#include <vic_driver_shared_all.h>

/******************************************************************************
 * @brief    Initialize the lake model state (energy balance, water balance,
 *           and snow components) to default values.
 *****************************************************************************/
void
generate_default_lake_state(lake_var_struct *lake,
                            soil_con_struct *soil_con,
                            lake_con_struct  lake_con)
{
    extern option_struct options;

    size_t               k;

    /************************************************************************
       Initialize lake state variables
       TBD: currently setting depth to depth_in from parameter file, but
            in future we should initialize to mindepth as default and
            eliminate depth_in (require user to use a state file if they
            want control over initial depth)
    ************************************************************************/
    if (options.LAKES) {
        lake->ldepth = lake_con.depth_in;
        for (k = 0; k < lake->activenod; k++) {
            // lake model requires FULL_ENERGY set to true
            lake->temp[k] = soil_con->avg_temp;
        }
    }
}
