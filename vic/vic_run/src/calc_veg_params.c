/******************************************************************************
 * @section DESCRIPTION
 *
 * Calculate vegetation parameters.
 *****************************************************************************/

#include <vic_run.h>

/******************************************************************************
 * @brief    This subroutine estimates the displacement height of vegetation
 *           with a given average height based on equations on page 4.12 of the
 *           Handbook of Hydrology.
 *****************************************************************************/
double
calc_veg_displacement(double height)
{
    extern parameters_struct param;

    double                   value;

    value = param.VEG_RATIO_DH_HEIGHT * height;

    return (value);
}

/******************************************************************************
 * @brief    This subroutine backs the vegetation height out of the given
 *           displacement using the reverse procedure from cal_veg_displacement.
 *****************************************************************************/
double
calc_veg_height(double displacement)
{
    extern parameters_struct param;

    double                   value;

    value = displacement / param.VEG_RATIO_DH_HEIGHT;

    return (value);
}

/******************************************************************************
 * @brief    This subroutine estimates the roughness height of vegetation with
 *           a given average height based on equations on page 4.12 of the
 *           Handbook of Hydrology.
 *****************************************************************************/
double
calc_veg_roughness(double height)
{
    extern parameters_struct param;

    double                   value;

    value = param.VEG_RATIO_RL_HEIGHT * height;

    return (value);
}
