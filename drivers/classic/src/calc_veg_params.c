#include <stdio.h>
#include <vicNl.h>
#include <math.h>

static char vcid[] = "$Id$";

double calc_veg_height(double displacement) {
/**********************************************************************
  calc_veg_height		Keith Cherkauer		March 3, 1997

  This subroutine backs the vegetation height out of the given
  displacement using the reverse procedure from cal_veg_displacement.
**********************************************************************/

  double value;

  value = displacement / 0.67;

  return (value);

}
