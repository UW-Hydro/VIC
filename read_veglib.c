#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

veg_lib_struct *read_veglib(FILE *veglib, int *Ntype)
/**********************************************************************
  This routine reads in a library of vegetation parameters for all
  vegetation classes used in the model.  The veg class number is used
  to reference the information in this library.
**********************************************************************/
{
  extern option_struct options;
  extern debug_struct debug;

  veg_lib_struct *temp;
  int    i, j;
  int    tmpflag;
  int    Nveg_type;
  char   str[MAXSTRING];
  double dum;

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
      fscanf(veglib, "%i",  &tmpflag);
      if(tmpflag==0) temp[i].overstory = FALSE;
      else temp[i].overstory = TRUE;
      fscanf(veglib, "%lf", &temp[i].rarc);
      fscanf(veglib, "%lf", &temp[i].rmin);
      dum=0.;
      for (j=0;j<options.Nlayer;j++) {
        fscanf(veglib, "%lf", &temp[i].root[j]);
        dum+=temp[i].root[j];
      }
      for (j=0;j<options.Nlayer;j++) {
        temp[i].root[j] /= dum;
      }
      if(dum == 0.0){
        sprintf(str,"Root fractions sum equals zero: %f , Vege Class: %d\n",
	        dum, temp[i].veg_class);
        nrerror(str);
      }
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
      for (j = 0; j < 12; j++) {
        fscanf(veglib, "%lf", &temp[i].displacement[j]);

        if(temp[i].LAI[j] > 0 && temp[i].displacement[j] <= 0) {
          sprintf(str,"Vegetation has leaves (LAI = %lf), but no displacement (%lf)",
	          temp[i].LAI[j], temp[i].displacement[j]);
          nrerror(str);
        }
        if(temp[i].albedo[j] < 0 || temp[i].albedo[j] > 1) {
          sprintf(str,"Albedo must be between 0 and 1 (%lf)",
	          temp[i].albedo[j]);
          nrerror(str);
        }
        temp[i].wind_h += calc_veg_height(temp[i].displacement[j]);
      }
      temp[i].wind_h *= 3./12.;
      fscanf(veglib, "%s", str);	/* skip over end of line comments */
      i++;
    }
    else fgets(str, MAXSTRING, veglib);
    fscanf(veglib, "%s", str);
  }
  if(i!=Nveg_type) nrerror("Error reading vegetation library file.\n");
  *Ntype = Nveg_type;

  return temp;
} 


