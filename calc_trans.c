#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

static char vcid[] = "$Id$";

double calc_trans(double deltat, double elevation)
/**********************************************************************
	calc_trans	Dag Lohmann		January 1996

  This routine computes the transmissity of the atmosphere for the
  current time step.

**********************************************************************/
{
  double trans, trans_clear;

  trans_clear = (A1_TRANS + A2_TRANS * elevation);
  trans = trans_clear * (1 - exp(B_TRANS * pow(deltat, C_TRANS)));

  return trans;
}
