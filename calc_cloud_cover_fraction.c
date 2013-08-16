#include <stdio.h>
#include <vicNl.h>

static char vcid[] = "$Id: calc_cloud_cover_fraction.c,v 4.1.2.1 2004/05/06 00:37:51 tbohn Exp $";

void calc_cloud_cover_fraction(atmos_data_struct *atmos,
			       dmy_struct        *dmy,
			       int                nrecs,
			       int                Ndays,
			       int                stepspday,
			       double            *tskc) {
/********************************************************************
  calc_cloud_cover_fraction     Keith Cherkauer    January 12, 2000

  This routine is designed to estimate cloud cover when observations
  of shortwave radiation are available.

*********************************************************************/

  nrerror("The function to estimate cloud cover from observed solar radiation does not yeat work.");

}
