#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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
           root fraction distribution information.               	KAC
  2009-Jun-09 Modified to use extension of veg_lib structure to contain
	      bare soil information, as well as other land cover types
	      used in potential evap calcualtions.			TJB
  2009-Oct-01 Added error message for case of LAI==0 and overstory==1.	TJB
  2010-Apr-28 Replaced GLOBAL_LAI with VEGPARAM_LAI and LAI_SRC.	TJB
  2012-Jan-16 Removed LINK_DEBUG code					BN
  2013-Jul-25 Added photosynthesis parameters.				TJB
  2014-Apr-25 Improved validation for LAI and albedo values.		TJB
  2014-Apr-25 Added partial vegcover fraction.				TJB
**********************************************************************/
{
  extern option_struct options;
  veg_lib_struct *temp;
  int    i, j;
  int    tmpflag;
  int    Nveg_type;
  char   str[MAXSTRING];
  char   ErrStr[MAXSTRING];
  double maxd;
  char   tmpstr[MAXSTRING];

  rewind(veglib);
  fgets(str,MAXSTRING,veglib);
  Nveg_type = 0;
  while(!feof(veglib)) {
    if(str[0]<=57 && str[0]>=48) Nveg_type++;
    fgets(str,MAXSTRING,veglib);
  }
  rewind(veglib);
      
  temp = (veg_lib_struct *)calloc(Nveg_type+N_PET_TYPES_NON_NAT,sizeof(veg_lib_struct));

  fscanf(veglib, "%s", str);
  i=0;
  while (!feof(veglib)) {
    if(str[0]<=57 && str[0]>=48) {
      temp[i].NVegLibTypes = Nveg_type;
      temp[i].veg_class = atoi(str);
      fscanf(veglib, "%d",  &tmpflag);
      if(tmpflag==0) temp[i].overstory = FALSE;
      else temp[i].overstory = TRUE;
      fscanf(veglib, "%lf", &temp[i].rarc);
      fscanf(veglib, "%lf", &temp[i].rmin);
      for (j = 0; j < 12; j++) {
        fscanf(veglib, "%lf", &temp[i].LAI[j]);
        if (options.LAI_SRC == LAI_FROM_VEGLIB && temp[i].overstory && temp[i].LAI[j] == 0) {
          sprintf(ErrStr,"ERROR: veg library: the specified veg class (%d) is listed as an overstory class, but the LAI given for this class for month %d is 0\n", temp[i].veg_class, j);
          nrerror(ErrStr);
        }
        temp[i].Wdmax[j] = LAI_WATER_FACTOR * temp[i].LAI[j];
      }
      if (options.VEGLIB_VEGCOVER) {
        for (j = 0; j < 12; j++) {
          fscanf(veglib, "%lf", &temp[i].vegcover[j]);
          if(temp[i].vegcover[j] < 0 || temp[i].vegcover[j] > 1) {
            sprintf(str,"Veg cover fraction must be between 0 and 1 (%f)", temp[i].vegcover[j]);
            nrerror(str);
          }
          if (temp[i].vegcover[j] < 0.01) temp[i].vegcover[j] = 0.01;
        }
      }
      else {
        for (j = 0; j < 12; j++) {
          temp[i].vegcover[j] = 1.00;
        }
      }
      for (j = 0; j < 12; j++) {
        fscanf(veglib, "%lf", &temp[i].albedo[j]);
        if(temp[i].albedo[j] < 0 || temp[i].albedo[j] > 1) {
          sprintf(str,"Albedo must be between 0 and 1 (%f)", temp[i].albedo[j]);
          nrerror(str);
        }
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
      /* Carbon-cycling parameters */
      if (options.VEGLIB_PHOTO) {
        fscanf(veglib, "%s", tmpstr); /* photosynthetic pathway */
        if (!strcmp(tmpstr,"C3")) temp[i].Ctype = PHOTO_C3;
        else if (!strcmp(tmpstr,"C4")) temp[i].Ctype = PHOTO_C4;
        fscanf(veglib, "%lf", &temp[i].MaxCarboxRate); /* Maximum carboxylation rate at 25 deg C */
        if (temp[i].Ctype == PHOTO_C3) {
          fscanf(veglib, "%lf", &temp[i].MaxETransport); /* Maximum electron transport rate at 25 deg C */
          temp[i].CO2Specificity = 0;
        }
        else if (temp[i].Ctype == PHOTO_C4) {
          fscanf(veglib, "%lf", &temp[i].CO2Specificity); /* CO2 Specificity */
          temp[i].MaxETransport = 0;
        }
        fscanf(veglib, "%lf", &temp[i].LightUseEff); /* Light-use efficiency */
        fscanf(veglib, "%s", tmpstr); /* Nitrogen-scaling flag */
        temp[i].NscaleFlag = atoi(tmpstr); /* Nitrogen-scaling flag */
        fscanf(veglib, "%lf", &temp[i].Wnpp_inhib); /* Moisture level in top soil layer above which photosynthesis begins experiencing inhibition due to saturation */
        fscanf(veglib, "%lf", &temp[i].NPPfactor_sat); /* photosynthesis multiplier when top soil layer is saturated */
      }
      else {
        temp[i].Wnpp_inhib = 1.0;
        temp[i].NPPfactor_sat = 1.0;
      }
      /* Irrigation parameters */
      if (options.VEGLIB_IRR) {
        for (j = 0; j < 12; j++) {
          fscanf(veglib, "%d", &temp[i].irr_active[j]);
          if(temp[i].irr_active[j] != 0 && temp[i].irr_active[j] != 1) {
            sprintf(str,"Irrigation active flag must be either 0 or 1 (%d)", temp[i].irr_active[j]);
            nrerror(str);
          }
        }
        fscanf(veglib, "%s", tmpstr); /* irrigation soil moisture threshold */
        if (!strcmp(tmpstr,"SAT")) temp[i].irr_sm_thresh = IRR_SAT;
        else if (!strcmp(tmpstr,"FC")) temp[i].irr_sm_thresh = IRR_FC;
        else if (!strcmp(tmpstr,"CR")) temp[i].irr_sm_thresh = IRR_CR;
              fscanf(veglib, "%s", tmpstr); /* irrigation soil moisture target */
        if (!strcmp(tmpstr,"SAT")) temp[i].irr_sm_target = IRR_SAT;
        else if (!strcmp(tmpstr,"FC")) temp[i].irr_sm_target = IRR_FC;
        else if (!strcmp(tmpstr,"CR")) temp[i].irr_sm_target = IRR_CR;
      }
      else {
        // Initialize but don't use
        for (j = 0; j < 12; j++) {
          temp[i].irr_active[j] = 0;
        }
        temp[i].irr_sm_thresh = IRR_SAT;
        temp[i].irr_sm_target = IRR_SAT;
      }

      // Default crop fractions
      for (j = 0; j < 12; j++) {
        temp[i].crop_frac[j] = 1.0;
      }


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
  for (i=0; i<N_PET_TYPES_NON_NAT; i++) {
    temp[Nveg_type+i].NVegLibTypes = Nveg_type;
    temp[Nveg_type+i].veg_class = Nveg_type+i+1;
    temp[Nveg_type+i].overstory = ref_veg_over[i];
    temp[Nveg_type+i].rarc = ref_veg_rarc[i];
    temp[Nveg_type+i].rmin = ref_veg_rmin[i];
    for (j=0; j<12; j++) {
      temp[Nveg_type+i].LAI[j] = ref_veg_lai[i];
      temp[Nveg_type+i].Wdmax[j] = LAI_WATER_FACTOR*ref_veg_lai[i];
      temp[Nveg_type+i].albedo[j] = ref_veg_albedo[i];
      temp[Nveg_type+i].roughness[j] = ref_veg_rough[i];
      temp[Nveg_type+i].displacement[j] = ref_veg_displ[i];
      temp[Nveg_type+i].vegcover[j] = ref_veg_vegcover[i];
      temp[Nveg_type+i].crop_frac[j] = ref_veg_crop_frac[i];
    }
    temp[Nveg_type+i].wind_h = ref_veg_wind_h[i];
    temp[Nveg_type+i].RGL = ref_veg_RGL[i];
    temp[Nveg_type+i].rad_atten = ref_veg_rad_atten[i];
    temp[Nveg_type+i].wind_atten = ref_veg_wind_atten[i];
    temp[Nveg_type+i].trunk_ratio = ref_veg_trunk_ratio[i];
  }

  return temp;
} 

/**********************************************************************
  free_veglib              Ted Bohn		April 2007 

  This routine frees the veglib structure.

  Modifications:

**********************************************************************/
void free_veglib(veg_lib_struct **veg_lib) {

  free((char*)(*veg_lib));
}
