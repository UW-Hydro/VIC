#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

veg_con_struct *read_vegparam(FILE *vegparam,
                              int   gridcel,
                              int   Nveg_type)
/**********************************************************************
  read_vegparam.c    Keith Cherkauer and Dag Lohmann       1997

  This routine reads in vegetation parameters for the current grid cell.
  It also relates each vegetation class in the cell to the appropriate
  parameters in the vegetation library.

  Modifications:
  09-24-98  Modified to read root zone distribution information so
           that soil layer root fractions can be computed for new 
	   soil layer depths - see calc_root_fractions.c           KAC

**********************************************************************/
{
  extern veg_lib_struct *veg_lib;
  extern option_struct   options;
  extern debug_struct    debug;

  veg_con_struct *temp;
  int             vegcel, i, j, vegetat_type_num;
  float           depth_sum;
  float           sum;
  char            str[500];

  rewind(vegparam);
      
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

  /** Allocate memory for vegetation grid cell parameters **/
  if(vegetat_type_num>0)
    temp = (veg_con_struct*) calloc(vegetat_type_num, 
                                    sizeof(veg_con_struct));
  else
    temp = (veg_con_struct*) calloc(1, sizeof(veg_con_struct));
  temp[0].Cv_sum = 0.0;

  for (i = 0; i < vegetat_type_num; i++) {
    temp[i].zone_depth = calloc(options.ROOT_ZONES,sizeof(float));
    temp[i].zone_fract = calloc(options.ROOT_ZONES,sizeof(float));
    temp[i].vegetat_type_num = vegetat_type_num;
    fscanf(vegparam, "%d",  &temp[i].veg_class);
    fscanf(vegparam, "%lf", &temp[i].Cv);
    
    depth_sum = 0;
    sum = 0.;
    for(j=0;j<options.ROOT_ZONES;j++) {
      fscanf(vegparam,"%f %f",&temp[i].zone_depth[j], &temp[i].zone_fract[j]);
      depth_sum += temp[i].zone_depth[j];
      sum += temp[i].zone_fract[j];
    }
    if(depth_sum <= 0) {
      sprintf(str,"Root zone depths must sum to a value greater than 0.");
      nrerror(str);
    }
    if(sum != 1.) {
      for(j=0;j<options.ROOT_ZONES;j++) {
	temp[i].zone_fract[j] /= sum;
      }
    }

    for(j=0;j<Nveg_type;j++)
      if(temp[i].veg_class == veg_lib[j].veg_class)
        temp[i].veg_class = j;

    temp[0].Cv_sum += temp[i].Cv;

  }
  if(temp[0].Cv_sum>1.0){
    fprintf(stderr,"WARNING: Cv exceeds 1.0 at grid cell %d, fractions being adjusted to equal 1\n", gridcel);
    for(j=0;j<vegetat_type_num;j++)
      temp[j].Cv = temp[j].Cv / temp[0].Cv_sum;
    temp[0].Cv_sum = 1.;
  }
  if(temp[0].Cv_sum>0.99 && temp[0].Cv_sum<1.0){
    fprintf(stderr,"WARNING: Cv > 0.99 and Cv < 1.0 at grid cell %d, model assuming that bare soil is not to be run - fractions being adjusted to equal 1\n", gridcel);
    for(j=0;j<vegetat_type_num;j++)
      temp[j].Cv = temp[j].Cv / temp[0].Cv_sum;
    temp[0].Cv_sum = 1.;
  }
  return temp;
} 


