#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

static char vcid[] = "$Id$";

veg_lib_struct *read_veglib(FILE *veglib, int *Ntype)
/**********************************************************************
  read_veglib.c               Keith Cherkauer                 1997

  This routine reads in a library of vegetation parameters for all
  vegetation classes used in the model.  The veg class number is used
  to reference the information in this library.

  Modifications:
  09-24-98 Modified to remove root fractions from the library file.
           See read_vegparam.c and calc_root_fraction.c for new
           root fraction distribution information.               KAC

**********************************************************************/
{
  extern option_struct options;
#if LINK_DEBUG
  extern debug_struct debug;
#endif

  veg_lib_struct *temp;
  int    i, j;
  int    tmpflag;
  int    Nveg_type;
  char   str[MAXSTRING];
  char   ErrStr[MAXSTRING];
  double maxd;

  rewind(veglib);
  fgets(str,MAXSTRING,veglib);
  Nveg_type = 0;
  while(!feof(veglib)) {
    if(str[0]<=57 && str[0]>=48) Nveg_type++;
    fgets(str,MAXSTRING,veglib);
  }
  rewind(veglib);
      
  temp = (veg_lib_struct *)calloc(Nveg_type,sizeof(veg_lib_struct));

  fscanf(veglib, "%s", str);
  i=0;
  while (!feof(veglib)) {
    if(str[0]<=57 && str[0]>=48) {
      temp[i].veg_class = atoi(str);
      fscanf(veglib, "%d",  &tmpflag);
      if(tmpflag==0) temp[i].overstory = FALSE;
      else temp[i].overstory = TRUE;
      fscanf(veglib, "%lf", &temp[i].rarc);
      fscanf(veglib, "%lf", &temp[i].rmin);
      for (j = 0; j < 12; j++) {
        fscanf(veglib, "%lf", &temp[i].LAI[j]);
        temp[i].Wdmax[j] = LAI_WATER_FACTOR * temp[i].LAI[j];
      }
      for (j = 0; j < 12; j++) {
        fscanf(veglib, "%lf", &temp[i].albedo[j]);
      }
      for (j = 0; j < 12; j++) {
        fscanf(veglib, "%lf", &temp[i].roughness[j]);
      }
      temp[i].wind_h = 0.;
      maxd = 0;
      for (j = 0; j < 12; j++) {
        fscanf(veglib, "%lf", &temp[i].displacement[j]);
        if(temp[i].displacement[j] > maxd) maxd = temp[i].displacement[j];
        if(temp[i].LAI[j] > 0 && temp[i].displacement[j] <= 0) {
          sprintf(str,"Vegetation has leaves (LAI = %f), but no displacement (%f)",
	          temp[i].LAI[j], temp[i].displacement[j]);
          nrerror(str);
        }
        if(temp[i].albedo[j] < 0 || temp[i].albedo[j] > 1) {
          sprintf(str,"Albedo must be between 0 and 1 (%f)",
	          temp[i].albedo[j]);
          nrerror(str);
        }
      }
      fscanf(veglib, "%lf", &temp[i].wind_h);
      if(temp[i].wind_h < maxd && temp[i].overstory) {
        sprintf(str,"Vegetation reference height (%f) for vegetation class %d, must be greater than the maximum displacement height (%f) when OVERSTORY has been set TRUE.",
                temp[i].wind_h,temp[i].veg_class,maxd);
        nrerror(str);
      }
      fscanf(veglib, "%f",  &temp[i].RGL);         /* minimum value of incoming
						    solar radiation at which there
						   will still be transpiration */
      if(temp[i].RGL < 0) {
        sprintf(str,"Minimum value of incoming solar radiation at which there is transpiration (RGL) must be greater than 0 for vegetation class %d.  Check that the vegetation library has the correct number of columns.",
                temp[i].veg_class);
        nrerror(str);
      }
      fscanf(veglib, "%lf", &temp[i].rad_atten);   /* vegetation radiation 
						      attenuation factor */
      if(temp[i].rad_atten < 0 || temp[i].rad_atten > 1) {
        sprintf(str,"The vegetation radiation attenuation factor must be greater than 0, and less than 1 for vegetation class %d.  Check that the vegetation library has the correct number of columns.",
                temp[i].veg_class);
        nrerror(str);
      }
      fscanf(veglib, "%lf", &temp[i].wind_atten);  /* canopy wind speed
						      attenuation factor */
      fscanf(veglib, "%lf", &temp[i].trunk_ratio); /* ratio of tree height that
						      is trunk */
      fgets(str, MAXSTRING, veglib);	/* skip over end of line comments */
      i++;
    }
    else fgets(str, MAXSTRING, veglib);
    fscanf(veglib, "%s", str);
  }
  if(i!=Nveg_type) {
    sprintf(ErrStr,"ERROR: Problem reading vegetation library file - make sure the file has the right number of columns.\n");
    nrerror(ErrStr);
  }
  *Ntype = Nveg_type;

  return temp;
} 


