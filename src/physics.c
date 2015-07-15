#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

double
calc_latent_heat_of_vaporization(double tair)
{
/**********************************************************************
   calc_latent_heat_of_vaporization      Joe Hamman        September 2, 2014

   This subroutine computes the latent heat of vaporization.  The simple 
   calculation is performed in a number of places within VIC and has been 
   moved to a function to improve consitancy and reusability.
   
   Eq. 4.2.1 in Handbook of Hydrology, assume Ts is Tair
**********************************************************************/

    double lv;

    lv = 2501000. - 2361. * tair;

    return(lv);
}

double
calc_latent_heat_of_sublimation(double temp)
{
/**********************************************************************
   calc_latent_heat_of_vaporization      Joe Hamman        September 2, 2014

   This subroutine computes the latent heat of sublimation.  The simple 
   calculation is performed in a number of places within VIC and has been 
   moved to a function to improve consitancy and reusability.
   
   Eq. 3.19, Bras 1990
**********************************************************************/

    double ls;

    ls = (677. - 0.07 * temp)  * JOULESPCAL * GRAMSPKG;

    return(ls);
}

double
calc_outgoing_longwave(double temp, double emiss)
{
/**********************************************************************
   calc_outgoing_longwave      Joe Hamman        September 2, 2014

   This subroutine computes outgoing longwave using the Stefan-Boltzman Law.
   The simple calculation is performed in a number of places within VIC and 
   has been moved to a function to improve consitancy and reusability.
   
**********************************************************************/

    double lwout;

    lwout = emiss*STEFAN_B*temp*temp*temp*temp;

    return(lwout);
}

double
calc_sensible_heat(double atmos_density, double t1, double t0, double Ra)
{
/**********************************************************************
   calc_sensible_heat      Joe Hamman        September 2, 2014

   This subroutine computes the sensible heat flux.
   The simple calculation is performed in a number of places within VIC and 
   has been moved to a function to improve consitancy and reusability.
   
**********************************************************************/

    double sensible;

    sensible = Cp * atmos_density * (t1 - t0) / Ra;

    return(sensible);
}
