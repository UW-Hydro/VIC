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
  extern debug_struct debug;

  int i;

  for (i = 0; i < nrecs; i++) {
    fprintf(debug.fg_atmos,"%d",  i);
    fprintf(debug.fg_atmos,"\t%lf", atmos[i].melt); 
    fprintf(debug.fg_atmos,"\t%lf", atmos[i].prec); 
    fprintf(debug.fg_atmos,"\t%lf", atmos[i].air_temp); 
    fprintf(debug.fg_atmos,"\t%lf", atmos[i].rainonly); 
    fprintf(debug.fg_atmos,"\t%lf", atmos[i].wind); 
    fprintf(debug.fg_atmos,"\t%lf", atmos[i].rad); 
    fprintf(debug.fg_atmos,"\t%lf", atmos[i].vpd); 
    fprintf(debug.fg_atmos,"\t%lf", atmos[i].vp); 
    fprintf(debug.fg_atmos,"\t%lf", atmos[i].pressure); 
    fprintf(debug.fg_atmos,"\t%lf", atmos[i].density); 
    fprintf(debug.fg_atmos,"\t%lf", atmos[i].rel_humid); 
    fprintf(debug.fg_atmos,"\t%lf", atmos[i].spec_humid); 
    fprintf(debug.fg_atmos,"\t%lf", atmos[i].tmin); 
    fprintf(debug.fg_atmos,"\t%lf", atmos[i].tmax); 
    fprintf(debug.fg_atmos,"\t%lf", atmos[i].priest); 
    fprintf(debug.fg_atmos,"\t%lf", atmos[i].penman_temp); 
    fprintf(debug.fg_atmos,"\t%lf", atmos[i].tskc); 
    fprintf(debug.fg_atmos,"\t%lf", atmos[i].trans); 
    fprintf(debug.fg_atmos,"\t%lf", atmos[i].shortwave); 
    fprintf(debug.fg_atmos,"\t%lf", atmos[i].longwave); 
    fprintf(debug.fg_atmos,"\t%lf", atmos[i].albedo); 
    fprintf(debug.fg_atmos,"\t%i",  atmos[i].rise_hour); 
    fprintf(debug.fg_atmos,"\t%i",  atmos[i].set_hour); 
    fprintf(debug.fg_atmos,"\n");
  }

  fflush(debug.fg_atmos);

}


