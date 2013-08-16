#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
 
static char vcid[] = "$Id: write_atmosdata.c,v 4.1.2.1 2004/05/10 18:44:40 tbohn Exp $";

void write_atmosdata(atmos_data_struct *atmos, int nrecs)
/**********************************************************************
	write_atmosdata		Dag Lohmann	Januray 1996

  This routine writes atmospheric data to the screen.

  Modifications:
    28-Aug-99 Changed to reflect the new atmos_data_struct.	Bart Nijssen
    07-May-04 No longer close the debug file, since the next cell
	      must write to it.					TJB
**********************************************************************/
{
#if LINK_DEBUG
  extern debug_struct debug;

  int i;
  int j;

  /*  first write all the SNOW_STEP data  - only write if the modelstep !=
      SNOWSTEP */
  if (NR > 0) {
    for (i = 0; i < nrecs; i++) {
      for (j = 0; j < NF; j++) {
	fprintf(debug.fg_snowstep_atmos,"%d\t%d",  i, j);
	fprintf(debug.fg_snowstep_atmos,"\t%f", atmos[i].prec[j]); 
	fprintf(debug.fg_snowstep_atmos,"\t%f", atmos[i].air_temp[j]); 
	fprintf(debug.fg_snowstep_atmos,"\t%f", atmos[i].wind[j]); 
	fprintf(debug.fg_snowstep_atmos,"\t%f", atmos[i].vpd[j]); 
	fprintf(debug.fg_snowstep_atmos,"\t%f", atmos[i].vp[j]); 
	fprintf(debug.fg_snowstep_atmos,"\t%f", atmos[i].pressure[j]); 
	fprintf(debug.fg_snowstep_atmos,"\t%f", atmos[i].density[j]); 
	fprintf(debug.fg_snowstep_atmos,"\t%f", atmos[i].shortwave[j]); 
	fprintf(debug.fg_snowstep_atmos,"\t%f", atmos[i].longwave[j]); 
	fprintf(debug.fg_snowstep_atmos,"\n");
      }
    }
  /* don't close the debug output file, as we need to write to it for the next cell as well */
/*  fclose(debug.fg_snowstep_atmos);*/
  }
  
  /* then write all the dt data */
  for (i = 0; i < nrecs; i++) {
    fprintf(debug.fg_modelstep_atmos,"%d",  i);
    fprintf(debug.fg_modelstep_atmos,"\t%f", atmos[i].prec[NR]); 
    fprintf(debug.fg_modelstep_atmos,"\t%f", atmos[i].air_temp[NR]); 
    fprintf(debug.fg_modelstep_atmos,"\t%f", atmos[i].wind[NR]); 
    fprintf(debug.fg_modelstep_atmos,"\t%f", atmos[i].vpd[NR]); 
    fprintf(debug.fg_modelstep_atmos,"\t%f", atmos[i].vp[NR]); 
    fprintf(debug.fg_modelstep_atmos,"\t%f", atmos[i].pressure[NR]); 
    fprintf(debug.fg_modelstep_atmos,"\t%f", atmos[i].density[NR]); 
    fprintf(debug.fg_modelstep_atmos,"\t%f", atmos[i].shortwave[NR]); 
    fprintf(debug.fg_modelstep_atmos,"\t%f", atmos[i].longwave[NR]); 
    fprintf(debug.fg_modelstep_atmos,"\n");
  }
  fclose(debug.fg_modelstep_atmos);
#endif

}


