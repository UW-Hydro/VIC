
#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

void write_vegparam(veg_con_struct *veg_con)
/**********************************************************************
	write_vegparam		Dag Lohmann	January 1996

  This routine writes vegitation parameters to the screen.

  Modifications:
  5/21/96	Routine was modified to allow for variable
		number of layers				KAC

**********************************************************************/
{
  extern veg_lib_struct *veg_lib;
  extern option_struct options;

  int i, j, l;

  printf("Vegitation Parameters:\n");
  for (i = 0; i < veg_con[0].vegetat_type_num; i++) {
    printf("\tvegetat_type_num = %d\n",  veg_con[i].vegetat_type_num);
    printf("\tveg_class = %d\n",  veg_lib[veg_con[i].veg_class].veg_class);
    printf("\tCv = %lf\n", veg_con[i].Cv);
    printf("\trarc = %lf\n", veg_lib[veg_con[i].veg_class].rarc);
    for(l=0;l<options.Nlayer;l++)
      printf("\troot_percent%d = %lf\n",l+1,
          veg_lib[veg_con[i].veg_class].root[l]);
    for (j = 0; j < 12; j++) 
      printf("\tLAI[%d] = %lf\n",j,veg_lib[veg_con[i].veg_class].LAI[j]);
    for (j = 0; j < 12; j++) 
      printf("\talbedo[%d] = %lf\n",j,
          veg_lib[veg_con[i].veg_class].albedo[j]);
    for (j = 0; j < 12; j++) 
      printf("\tdisplacement[%d] = %lf\n\n",j,
          veg_lib[veg_con[i].veg_class].displacement[j]);
    for (j = 0; j < 12; j++) 
      printf("\troughness[%d] = %lf\n\n",j,
          veg_lib[veg_con[i].veg_class].roughness[j]);
  }
}

