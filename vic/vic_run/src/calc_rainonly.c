/******************************************************************************
 * @section DESCRIPTION
 *
 * Determines from the air temperature what fraction of incoming precipitation
 * is frozen and unfrozen (snow and rain).
 *****************************************************************************/

#include <vic_run.h>

/******************************************************************************
 * @brief    Determines from the air temperature what fraction of incoming
 *           precipitation is frozen and unfrozen (snow and rain).
 *****************************************************************************/
double
calc_rainonly(double air_temp,
              double prec,
              double MAX_SNOW_TEMP,
              double MIN_RAIN_TEMP)
{
    double rainonly;

    rainonly = 0.;
    if (MAX_SNOW_TEMP <= MIN_RAIN_TEMP) {
        log_err("MAX_SNOW_TEMP must be greater then MIN_RAIN_TEMP");
    }
    if (air_temp < MAX_SNOW_TEMP && air_temp > MIN_RAIN_TEMP) {
        rainonly = (air_temp - MIN_RAIN_TEMP) /
                   (MAX_SNOW_TEMP - MIN_RAIN_TEMP) * prec;
    }
    else if (air_temp >= MAX_SNOW_TEMP) {
        rainonly = prec;
    }

    return(rainonly);
}
