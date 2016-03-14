/******************************************************************************
 * @section DESCRIPTION
 *
 * Iteratively solve the atmospheric energy balance equation to estimate the
 * canopy air temperature.
 *
 * Basic concept for this was taken from:
 *    Sellers et al., J. Clim., v.9, April 1996, pp. 676-705.
 *    Dickinsen, BATS manual, NCAR Tech. Note (NCAR/TN-387+STR), August 1993.
 *
 * @section LICENSE
 *
 * The Variable Infiltration Capacity (VIC) macroscale hydrological model
 * Copyright (C) 2014 The Land Surface Hydrology Group, Department of Civil
 * and Environmental Engineering, University of Washington.
 *
 * The VIC model is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 *****************************************************************************/

#include <vic_run.h>

/******************************************************************************
 * @brief    Iteratively solve the atmospheric energy balance equation to
 *           estimate the canopy air temperature.
 *****************************************************************************/
double
calc_atmos_energy_bal(double    InOverSensible,
                      double    InUnderSensible,
                      double    LatentHeatOver,
                      double    LatentHeatUnder,
                      double    LatentHeatSubOver,
                      double    LatentHeatSubUnder,
                      double    NetLongOver,
                      double    NetLongUnder,
                      double    NetShortOver,
                      double    NetShortUnder,
                      double    Ra,
                      double    Tair,
                      double    atmos_density,
                      double   *Error,
                      double   *LatentHeat,
                      double   *LatentHeatSub,
                      double   *NetLongAtmos,
                      double   *NetShortAtmos,
                      double   *SensibleHeat,
                      bool     *Tcanopy_fbflag,
                      unsigned *Tcanopy_fbcount)
{
    extern option_struct     options;
    extern parameters_struct param;

    double                   F; // canopy closure fraction, not currently used by VIC
    double                   InSensible;
    double                   NetRadiation;
    double                   T_lower;
    double                   T_upper;
    double                   Tcanopy;
    char                     ErrorString[MAXSTRING];

    F = 1;

    // compute incoming sensible heat
    InSensible = InOverSensible + InUnderSensible;
    (*SensibleHeat) = InOverSensible + InUnderSensible;

    // compute net radiation
    (*NetLongAtmos) = (F * NetLongOver + (1. - F) * NetLongUnder);

    (*NetShortAtmos) = (NetShortOver + NetShortUnder);

    NetRadiation = (*NetShortAtmos + *NetLongAtmos);

    // compute total latent heat flux
    (*LatentHeat) = (LatentHeatOver + LatentHeatUnder);

    (*LatentHeatSub) = (LatentHeatSubOver + LatentHeatSubUnder);

    /******************************
       Find Canopy Air Temperature
    ******************************/

    if (options.CLOSE_ENERGY) {
        /* initialize Tcanopy_fbflag */
        *Tcanopy_fbflag = 0;

        /* set initial bounds for root brent **/
        T_lower = (Tair) - param.CANOPY_DT;
        T_upper = (Tair) + param.CANOPY_DT;

        // iterate for canopy air temperature
        Tcanopy = root_brent(T_lower, T_upper, ErrorString,
                             func_atmos_energy_bal, Ra, Tair, atmos_density,
                             InSensible, SensibleHeat);

        if (Tcanopy <= -998) {
            if (options.TFALLBACK) {
                Tcanopy = Tair;
                *Tcanopy_fbflag = 1;
                (*Tcanopy_fbcount)++;
            }
            else {
                // handle error flag from root brent
                (*Error) = error_calc_atmos_energy_bal(Tcanopy, (*LatentHeat) +
                                                       (*LatentHeatSub),
                                                       NetRadiation, Ra, Tair,
                                                       atmos_density,
                                                       InSensible,
                                                       SensibleHeat,
                                                       ErrorString);
                return (ERROR);
            }
        }
    }
    else {
        Tcanopy = Tair;
    }

    // compute variables based on final temperature
    (*Error) = solve_atmos_energy_bal(Tcanopy, Ra, Tair, atmos_density,
                                      InSensible, SensibleHeat);
    return(Tcanopy);
}

/******************************************************************************
 * @brief    Dummy function to allow calling func_atmos_energy_bal() directly.
 *****************************************************************************/
double
solve_atmos_energy_bal(double Tcanopy,
                       ...)
{
    va_list ap;

    double  error;

    va_start(ap, Tcanopy);
    error = func_atmos_energy_bal(Tcanopy, ap);
    va_end(ap);

    return error;
}

/******************************************************************************
 * @brief    Dummy function to allow calling error_calc_atmos_energy_bal()
 *           directly.
 *****************************************************************************/
double
error_calc_atmos_energy_bal(double Tcanopy,
                            ...)
{
    va_list ap;

    double  error;

    va_start(ap, Tcanopy);
    error = error_print_atmos_energy_bal(Tcanopy, ap);
    va_end(ap);

    return error;
}

/******************************************************************************
 * @brief    Print atmos energy balance terms.
 *****************************************************************************/
double
error_print_atmos_energy_bal(double  Tcanopy,
                             va_list ap)
{
    double  LatentHeat;
    double  NetRadiation;
    double  Ra;
    double  Tair;
    double  atmos_density;
    double  InSensible;

    double *SensibleHeat;
    char   *ErrorString;

    // extract variables from va_arg
    LatentHeat = (double)  va_arg(ap, double);
    NetRadiation = (double)  va_arg(ap, double);
    Ra = (double)  va_arg(ap, double);
    Tair = (double)  va_arg(ap, double);
    atmos_density = (double)  va_arg(ap, double);
    InSensible = (double)  va_arg(ap, double);

    SensibleHeat = (double *)va_arg(ap, double *);
    ErrorString = (char *)va_arg(ap, char *);

    // print variable values
    fprintf(LOG_DEST, "%s", ErrorString);
    fprintf(LOG_DEST, "ERROR: calc_atmos_energy_bal failed to converge to a "
            "solution in root_brent.  Variable values will be dumped "
            "to the screen, check for invalid values.\n");
    fprintf(LOG_DEST, "Tcanopy = %f\n", Tcanopy);
    fprintf(LOG_DEST, "LatentHeat = %f\n", LatentHeat);
    fprintf(LOG_DEST, "NetRadiation = %f\n", NetRadiation);
    fprintf(LOG_DEST, "Ra = %f\n", Ra);
    fprintf(LOG_DEST, "Tair = %f\n", Tair);
    fprintf(LOG_DEST, "atmos_density = %f\n", atmos_density);
    fprintf(LOG_DEST, "InSensible = %f\n", InSensible);

    fprintf(LOG_DEST, "*SensibleHeat = %f\n", *SensibleHeat);

    fprintf(LOG_DEST, "Finished writing calc_atmos_energy_bal variables.\n"
            "Try increasing CANOPY_DT to get model to complete cell.\n"
            "Then check output for instabilities.\n");

    return(ERROR);
}

/******************************************************************************
 * @brief   Dummy function to allow calling func_atmos_moist_bal() directly.
 *****************************************************************************/
double
solve_atmos_moist_bal(double VPcanopy,
                      ...)
{
    va_list ap;

    double  error;

    va_start(ap, VPcanopy);
    error = func_atmos_moist_bal(VPcanopy, ap);
    va_end(ap);

    return error;
}

/******************************************************************************
 * @brief    Dummy function to allow calling error_print_atmos_moist_bal()
 *           directly.
 *****************************************************************************/
double
error_calc_atmos_moist_bal(double VPcanopy,
                           ...)
{
    va_list ap;

    double  error;

    va_start(ap, VPcanopy);
    error = error_print_atmos_moist_bal(VPcanopy, ap);
    va_end(ap);

    return error;
}

/******************************************************************************
 * @brief    Print atmos moist energy balance terms.
 *****************************************************************************/
double
error_print_atmos_moist_bal(double  VPcanopy,
                            va_list ap)
{
    double  InLatent;
    double  Lv;
    double  Ra;
    double  atmos_density;
    double  gamma;
    double  vp;
    double *AtmosLatent;
    char   *ErrorString;

    // extract variables from va_arg
    InLatent = (double)  va_arg(ap, double);
    Lv = (double)  va_arg(ap, double);
    Ra = (double)  va_arg(ap, double);
    atmos_density = (double)  va_arg(ap, double);
    gamma = (double)  va_arg(ap, double);
    vp = (double)  va_arg(ap, double);
    AtmosLatent = (double *)va_arg(ap, double *);
    ErrorString = (char *)  va_arg(ap, char *);

    // print variable values
    fprintf(LOG_DEST, "%s", ErrorString);
    fprintf(LOG_DEST, "VPcanopy = %f\n", VPcanopy);
    fprintf(LOG_DEST, "InLatent = %f\n", InLatent);
    fprintf(LOG_DEST, "Lv = %f\n", Lv);
    fprintf(LOG_DEST, "Ra = %f\n", Ra);
    fprintf(LOG_DEST, "atmos_density = %f\n", atmos_density);
    fprintf(LOG_DEST, "gamma = %f\n", gamma);
    fprintf(LOG_DEST, "vp = %f\n", vp);
    fprintf(LOG_DEST, "AtmosLatent = %f\n", *AtmosLatent);

    log_err("Finished writing calc_atmos_moist_bal variables.\nTry increasing "
            "CANOPY_VP to get model to complete cell.\nThen check output for "
            "instabilities.");

    return(0.0);
}
