/*
 * SUMMARY:      WaterDensity.c - Calculate the density of water
 * USAGE:        Part of DHSVM
 *
 * AUTHOR:       Bart Nijssen
 * ORG:          University of Washington, Department of Civil Engineering
 * E-MAIL:       nijssen@u.washington.edu
 * ORIG-DATE:    Apr-1996
 * LAST-MOD: Thu Oct  1 18:18:27 1998 by VIC Administrator <vicadmin@u.washington.edu>
 * DESCRIPTION:  Calculate the density of water as a function of temperature
 * DESCRIP-END.
 * FUNCTIONS:    WaterDensity()
 * COMMENTS:     
 */

#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

static char vcid[] = "$Id$";

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
  double rho;

  rho = 1000 * (1 - (T + 288.9414)/(508929.2 * (T + 68.12963)) * 
                pow((double) (T - 3.9863), (double) 2.0));

  return rho;
}
