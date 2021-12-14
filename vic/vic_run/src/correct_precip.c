/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine corrects preciptation measurements for gauge catch
 * deficiencies.
 *
 * Correction values read from Bras Figure 4.16.
 *
 * NOTE: Should locate better reference with fitted equations, or at least data
 * with which to fit an equation.
 *****************************************************************************/

#include <vic_run.h>

/******************************************************************************
 * @brief    Correct preciptation measurements for gauge catch deficiencies.
 *****************************************************************************/
void
correct_precip(double *gauge_correction,
               double  wind,
               double  wind_h,
               double  roughness,
               double  snow_roughness)
{
    extern parameters_struct param;

    double                   gauge_wind;

    gauge_wind = wind * (log((param.GAUGE_HEIGHT + roughness) / roughness) /
                         log(wind_h / roughness));

    gauge_correction[RAIN] = 100. / exp(4.606 - 0.041 *
                                        pow(gauge_wind, 0.69));

    gauge_wind = wind * (log((param.GAUGE_HEIGHT + snow_roughness) /
                             snow_roughness) /
                         log(wind_h / snow_roughness));

    gauge_correction[SNOW] = 100. / exp(4.606 - 0.036 *
                                        pow(gauge_wind, 1.75));
}
