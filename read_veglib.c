#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

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
      }
      fscanf(veglib, "%lf", &temp[i].wind_h);
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


