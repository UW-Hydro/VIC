#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
 
out_data_struct *make_out_data(int nrecs)
/**********************************************************************
	make_out_data	Dag Lohmann		January 1996

  This routine creates an array of out data structures, one for each
  time step.

**********************************************************************/
{
  out_data_struct *temp;

  temp = (out_data_struct*) calloc(nrecs, sizeof(out_data_struct));
  return temp;
}
