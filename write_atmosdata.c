#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
 
static char vcid[] = "$Id$";

void write_atmosdata(atmos_data_struct *atmos, int nrecs)
/**********************************************************************
	write_atmosdata		Dag Lohmann	Januray 1996

  This routine writes atmospheric data to the screen.

**********************************************************************/
{
  int i;

  printf("Atmospheric Data:\n");
  for (i = 0; i < 10; i++) {
    printf("\ti = %d\n", i);
    printf("\tMelt            = %lf\n", atmos[i].melt); 
    printf("\tPrec            = %lf\n", atmos[i].prec); 
    printf("\tAir_temp        = %lf\n", atmos[i].air_temp); 
    printf("\trainonly        = %lf\n", atmos[i].rainonly); 
    printf("\tWind            = %lf\n", atmos[i].wind); 
    printf("\talbedo          = %lf\n", atmos[i].albedo); 
    printf("\ttmin            = %lf\n", atmos[i].tmin); 
    printf("\ttmax            = %lf\n\n", atmos[i].tmax); 
  }
}


