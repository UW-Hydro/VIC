#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

void write_vegparam(veg_con_struct *veg_con)
/**********************************************************************
	write_vegparam		Dag Lohmann	January 1996

  This routine writes vegetation parameters to stdout, used primarily 
  for debugging, and making sure the model is reading the proper 
  parameters..

  Modifications:
  5/21/96	Routine was modified to allow for variable
		number of layers				KAC
  4-12-98  Updated for new standard vegetation parameters       KAC

**********************************************************************/
{
  extern veg_lib_struct *veg_lib;
  extern option_struct options;

  int i, j, l;
  int vegclass;

  printf("Vegetation Parameters:\n");
  for (i = 0; i < veg_con[0].vegetat_type_num; i++) {
    vegclass = veg_con[i].veg_class;
    printf("\tvegetat_type_num = %d\n",  veg_con[i].vegetat_type_num);
    printf("\tveg_class        = %d\n",  veg_lib[vegclass].veg_class);
    printf("\tCv               = %lf\n", veg_con[i].Cv);
    if(veg_lib[vegclass].overstory)
      printf("\tOverstory        = TRUE\n");
    else 
      printf("\tOverstory        = FALSE\n");
    printf("\trarc             = %lf s/m\n", veg_lib[vegclass].rarc);
    printf("\trmin             = %lf s/m\n", veg_lib[vegclass].rmin);
    for(l=0;l<options.Nlayer;l++)
      printf("\troot_percent%d   = %lf\n",l+1,veg_con[i].root[l]);
    for (j = 0; j < 12; j++) 
      printf("\tLAI[%d]          = %lf\n",j,veg_lib[vegclass].LAI[j]);
    for (j = 0; j < 12; j++) 
      printf("\talbedo[%d]       = %lf\n",j,veg_lib[vegclass].albedo[j]);
    for (j = 0; j < 12; j++) 
      printf("\tdisplacement[%d] = %lf m\n",j,
	     veg_lib[vegclass].displacement[j]);
    for (j = 0; j < 12; j++) 
      printf("\troughness[%d]    = %lf m\n",j,veg_lib[vegclass].roughness[j]);
    printf("\twind_h           = %lf s/m\n", veg_lib[vegclass].wind_h);
  }
}

