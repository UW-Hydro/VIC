#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

void write_vegvar(veg_var_struct veg, int n)
/**********************************************************************
	write_vegvar		Keith Cherkauer		May 29, 1996

  This routine writes vegetation variables to stdout.  Used primarily
  for debugging purposes.

  Modifications:
  5/21/96	Routine was modified to allow for variable
		number of layers				KAC

**********************************************************************/
{
  int i, j, l;

  printf("Vegetation Variables: vegtype %i\n",n);
  printf("\tcanopyevap  = %lf\n", veg.canopyevap);
  printf("\tWdew        = %lf\n", veg.Wdew);
  printf("\tthroughfall = %lf\n", veg.throughfall);
}

