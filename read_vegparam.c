#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

static char vcid[] = "$Id$";

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
  07-15-99 Modified to read LAI values from a new line in the vegetation
           parameter file.  Added specifically to work with the new
	   global LAI files.
  11-18-02 Added code to read in blowing snow parameters.          LCB
  03-27-03 Modified code to update Wdmax based on LAI values read in
           for the current grid cell.  If LAI is not obtained from this
           function, then the values cacluated in read_veglib.c are
           left unchanged.                                   DP & KAC

**********************************************************************/
{
  extern veg_lib_struct *veg_lib;
  extern option_struct   options;
#if LINK_DEBUG
  extern debug_struct    debug;
#endif

  veg_con_struct *temp;
  int             vegcel, i, j, vegetat_type_num, skip, veg_class;
  int             MaxVeg;
  float           depth_sum;
  float           sum;
  char            str[500];
  char            ErrStr[MAXSTRING];

  if(options.GLOBAL_LAI) skip=2;
  else skip=1;

#if !NO_REWIND
  rewind(vegparam);
#endif  
    
  while ((fscanf(vegparam, "%d %d", &vegcel, &vegetat_type_num)) == 2 &&
          vegcel != gridcel) {
    for (i = 0; i <= vegetat_type_num * skip; i++)
      fgets(str, 500, vegparam);
  }
  if (vegcel != gridcel) {
    fprintf(stderr, "Error in vegetation file.  Grid cell %d not found\n",
            gridcel);
    exit(99);
  }
  if(vegetat_type_num >= MAX_VEG) {
    sprintf(ErrStr,"Vegetation parameter file wants more vegetation types in grid cell %i (%i) than are defined by MAX_VEG (%i) [NOTE: bare soil class is assumed].  Edit vicNl_def.h and recompile.",gridcel,vegetat_type_num+1,MAX_VEG);
    nrerror(ErrStr);
  }

  // Make sure to allocate extra memory for wetlands vegetation
  if ( options.LAKES )
    MaxVeg = vegetat_type_num+1;
  else
    MaxVeg = vegetat_type_num;

  /** Allocate memory for vegetation grid cell parameters **/
  if ( MaxVeg > 0 )
    temp = (veg_con_struct*) calloc( MaxVeg, sizeof(veg_con_struct));
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
      fprintf(stderr,"WARNING: Root zone fractions sum to more than 1 ( = %f), normalizing fractions.  If the sum is large, check that your vegetation parameter file is in the form - <zone 1 depth> <zone 1 fract> <zone 2 depth> <zone 2 fract> ...\n", sum);
      for(j=0;j<options.ROOT_ZONES;j++) {
	temp[i].zone_fract[j] /= sum;
      }
    }

    if(options.BLOWING) {
       fscanf(vegparam,"%f %f %f",&temp[i].sigma_slope, &temp[i].lag_one, &temp[i].fetch);
      if( temp[i].sigma_slope <= 0. || temp[i].lag_one <= 0.) {
      sprintf(str,"Deviation of terrain slope must be greater than 0.");
       nrerror(str);
      }
    }

    veg_class = MISSING;
    for(j=0;j<Nveg_type;j++)
      if(temp[i].veg_class == veg_lib[j].veg_class)
	veg_class = j;
    if(veg_class == MISSING) {
      sprintf(ErrStr,"Vegetation class %i from cell %i is not defined in the vegetation library file.", temp[i].veg_class, gridcel);
      nrerror(ErrStr);
    }
    else
      temp[i].veg_class = veg_class;

    temp[0].Cv_sum += temp[i].Cv;

    if ( options.GLOBAL_LAI )
      for ( j = 0; j < 12; j++ ) {
	fscanf(vegparam,"%lf",&veg_lib[temp[i].veg_class].LAI[j]);
	veg_lib[temp[i].veg_class].Wdmax[j] = 
	  LAI_WATER_FACTOR * veg_lib[temp[i].veg_class].LAI[j];
      }
    
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


