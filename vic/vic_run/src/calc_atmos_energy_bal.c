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
 * Copyright (C) 2016 The Computational Hydrology Group, Department of Civil
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
        Tcanopy = root_brent(T_lower, T_upper,
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
                                                       SensibleHeat);
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

    // extract variables from va_arg
    LatentHeat = (double)  va_arg(ap, double);
    NetRadiation = (double)  va_arg(ap, double);
    Ra = (double)  va_arg(ap, double);
    Tair = (double)  va_arg(ap, double);
    atmos_density = (double)  va_arg(ap, double);
    InSensible = (double)  va_arg(ap, double);

    SensibleHeat = (double *)va_arg(ap, double *);

    // print variable values
    log_warn("Failure to converge to a solution in root_brent.\n"
             "Check for invalid values.\n"
             "Tcanopy = %f\n"
             "LatentHeat = %f\n"
             "NetRadiation = %f\n"
             "Ra = %f\n"
             "Tair = %f\n"
             "atmos_density = %f\n"
             "InSensible = %f\n"
             "*SensibleHeat = %f\n"
             "Try increasing CANOPY_DT to get model to complete cell.\n"
             "Then check output for instabilities.",
             Tcanopy, LatentHeat, NetRadiation, Ra, Tair, atmos_density,
             InSensible, *SensibleHeat);

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

    // extract variables from va_arg
    InLatent = (double)  va_arg(ap, double);
    Lv = (double)  va_arg(ap, double);
    Ra = (double)  va_arg(ap, double);
    atmos_density = (double)  va_arg(ap, double);
    gamma = (double)  va_arg(ap, double);
    vp = (double)  va_arg(ap, double);
    AtmosLatent = (double *)va_arg(ap, double *);

    // print variable values
    log_err("VPcanopy = %f\n"
            "InLatent = %f\n"
            "Lv = %f\n"
            "Ra = %f\n"
            "atmos_density = %f\n"
            "gamma = %f\n"
            "vp = %f\n"
            "AtmosLatent = %f\n"
            "Try increasing CANOPY_VP to get model to complete cell.\n"
            "Then check output for instabilities.",
            VPcanopy, InLatent, Lv, Ra, atmos_density, gamma, vp,
            *AtmosLatent);

    return(0.0);
}
