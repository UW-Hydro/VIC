/******************************************************************************
 * @section DESCRIPTION
 *
 * This set of functions calculate basic physical quantities that are
 * frequently calculated throughout the model
 *****************************************************************************/

#include <vic_run.h>

/******************************************************************************
 * @brief    Compute the latent heat of sublimation.
 * @note     Eq. 3.19, Bras 1990
 *****************************************************************************/
double
calc_latent_heat_of_sublimation(double temp)
{
    double ls;

    ls = (677. - 0.07 * temp) * JOULES_PER_CAL * GRAMS_PER_KG;

    return(ls);
}

/******************************************************************************
 * @brief    This subroutine computes the latent heat of vaporization.
 * @note     Eq. 4.2.1 in Handbook of Hydrology.
 *****************************************************************************/
double
calc_latent_heat_of_vaporization(double temp)
{
    double lv;

    lv = CONST_LATVAP - 2361. * temp;

    return(lv);
}

/******************************************************************************
 * @brief    Compute outgoing longwave using the Stefan-Boltzman Law.
 *****************************************************************************/
double
calc_outgoing_longwave(double temp,
                       double emiss)
{
    double lwout;

    lwout = emiss * CONST_STEBOL * temp * temp * temp * temp;

    return(lwout);
}

/******************************************************************************
 * @brief    Calculate scale height based on average temperature in the column
 *****************************************************************************/
double
calc_scale_height(double tair,
                  double elevation)
{
    extern parameters_struct param;

    double                   h;

    h = CONST_RDAIR / CONST_G *
        ((tair + CONST_TKFRZ) + 0.5 * elevation * param.LAPSE_RATE);

    return(h);
}

/******************************************************************************
 * @brief    Compute the sensible heat flux.
 *****************************************************************************/
double
calc_sensible_heat(double atmos_density,
                   double t1,
                   double t0,
                   double Ra)
{
    double sensible;

    sensible = CONST_CPMAIR * atmos_density * (t1 - t0) / Ra;

    return(sensible);
}
