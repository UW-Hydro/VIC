
#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

veg_con_struct *read_vegparam(FILE *vegparam,
                              int gridcel,
                              int Nveg_type)
/**********************************************************************
  This routine reads in vegetation parameters for the current grid cell.
  It also relates each vegetation class in the cell to the appropriate
  parameters in the vegetation library.
**********************************************************************/
{
  extern veg_lib_struct *veg_lib;
  extern option_struct options;
  extern debug_struct debug;

  veg_con_struct *temp;
  int    vegcel, i, j, vegetat_type_num;
  int    startlayer;
  char   str[500];
  rewind(vegparam);
      
  if(options.FULL_ENERGY || options.FROZEN_SOIL) startlayer=1;
  else startlayer=0;

  while ((fscanf(vegparam, "%d %d", &vegcel, &vegetat_type_num)) == 2 &&
          vegcel != gridcel) {
    for (i = 0; i <= vegetat_type_num; i++)
      fgets(str, 500, vegparam);
  }
  if (vegcel != gridcel) {
    fprintf(stderr, "Error in vegetation file.  Grid cell %d not found\n",
            gridcel);
    exit(99);
  }
  temp = (veg_con_struct*) calloc(vegetat_type_num, 
                                  sizeof(veg_con_struct));
  temp[0].Cv_sum = 0.0;
  for (i = 0; i < vegetat_type_num; i++) {
    temp[i].vegetat_type_num = vegetat_type_num;
    fscanf(vegparam, "%d",  &temp[i].veg_class);
    fscanf(vegparam, "%lf", &temp[i].Cv);
    
    for(j=0;j<Nveg_type;j++)
      if(temp[i].veg_class == veg_lib[j].veg_class)
        temp[i].veg_class = j;

    temp[0].Cv_sum += temp[i].Cv;

  }
  if(temp[0].Cv_sum>1.0){
    fprintf(stderr,"WARNING: Cv exceeds 1.0 at grid cell %d, fractions being adjusted to equal 1\n", gridcel);
    for(j=0;j<Nveg_type;j++)
      temp[i].Cv = temp[i].Cv / temp[0].Cv_sum;
    temp[0].Cv_sum = 1.;
  }
  return temp;
} 


