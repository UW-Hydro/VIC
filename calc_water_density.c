/*
 * SUMMARY:      WaterDensity.c - Calculate the density of water
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-1996
 * LAST-MOD:     23-Aug-1996 at 16:58:14 by DHSVM Project Account
 * DESCRIPTION:  Calculate the density of water as a function of temperature
 * DESCRIP-END.
 * FUNCTIONS:    WaterDensity()
 * COMMENTS:     
 */

#ifndef lint
static char vcid[] = "nil";
#endif /* lint */

#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

/*
#include "settings.h"
#include "density.h"
*/

/*****************************************************************************
  Function name: WaterDensity()

  Purpose      : Calculate the density of water as a function of temperature 
                 using the Thiesen-Scheel-Diesselhorst equation

  Required     : 
    double T - Temperature in C

  Returns      : rho - water density in kg/m3

  Modifies     : NA

  Comments     : Reference: McCutcheon, S. C., J. L. Martin, and 
                   T. O. Barnwell, Jr, Water quality, 
                   In: Maidment, D. R. (ed.), Handbook of hydrology, 
                   1993, McGraw-Hill, New York, etc.. (p. 11.6, Fig. 11.1.1)
*****************************************************************************/
double WaterDensity(double T)
{
  const char *Routine = "WaterDensity";
  double rho;

  rho = 1000 * (1 - (T + 288.9414)/(508929.2 * (T + 68.12963)) * 
                pow((double) (T - 3.9863), (double) 2.0));

  return rho;
}
