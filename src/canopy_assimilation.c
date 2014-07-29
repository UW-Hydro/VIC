/********************************************************************************
  filename  : canopy_assimilation.c
  purpose   : Calculate GPP, Raut, and NPP for veg cover with multi-layer canopy
  interface : - input :
                - Ctype:          photosynthesis pathway (PHOTO_C3 or PHOTO_C4)
                - MaxCarboxRate:  maximum carboxlyation rate at 25 deg C       (mol(CO2)/m2 leaf area s)
                - MaxETransport:  maximum electron transport rate at 25 deg C  (mol(CO2)/m2 leaf area s) (C3 plants)
                - CO2Specificity: CO2 specificity at 25 deg C                  (mol(CO2)/m2 leaf area s) (C4 plants)
                - NscaleFactor[]: array of layer-specific nitrogen scaling
                                  factors at max carbox rate (Vm) and max
                                  electron transport rate (Jm)
                - Tfoliage:       vegetation temperature                       (deg C)
                - SWdown:         incoming shortwave radiation                 (W/m2)
                - aPAR[]:         array of layer-specific absorbed
                                  photosynthetically active radiation (mol
                                  photons/m2 leaf area s)
                - pz:          near-surface atmospheric pressure            (Pa)
                - Catm:           CO2 mixing ratio in the atmosphere           (mol(CO2)/mol(air))
                - CanopLayerBnd[]:array of canopy layer upper boundaries,
                                  expressed in terms of fraction of total LAI
                - LAItotal:       Total LAI of the entire canopy
                - mode:           'ci': assume a default Ci, and compute photosynthesis and rs
                                  'rs': take rs as an input, and compute photosynthesis and Ci

              - input or output, depending on value of mode:
                - rsLayer[]:      array of layer-specific stomatal resistance  (s/m) (mode = 'rs': input; mode = 'ci': output)
                                  NOTE: this is per layer leaf area
                - rc:             whole-canopy stomatal resistance             (s/m) (mode = 'rs': ignored; mode = 'ci': output)

              - output:
                - Ci:             whole-canopy leaf-internal CO2 mixing ratio      (mol(CO2)/mol(air))
                - GPP:            whole-canopy gross assimilation (photosynthesis) (mol(CO2)/m2 ground area s)
                - Rdark:          whole-canopy 'dark' of leaf respiration          (mol(CO2)/m2 ground area s)
                - Rphoto:         whole-canopy photorespiration                    (mol(CO2)/m2 ground area s)
                - Rmaint:         whole-forest maintenance respiration             (mol(CO2)/m2 ground area s)
                - Rgrowth:        whole-forest plant growth respiration            (mol(CO2)/m2 ground area s)
                - Raut:           whole-forest plant respiration                   (mol(CO2)/m2 ground area s)
                - NPP:            net primary productivity                         (mol(CO2)/m2 ground area s)

  programmer: Ted Bohn
  date      : October 20, 2006
  changes   :
  references: 
********************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vicNl.h>

static char vcid[] = "$Id: $";

void canopy_assimilation(char    Ctype,
                         double  MaxCarboxRate,
                         double  MaxETransport,
                         double  CO2Specificity,
                         double *NscaleFactor,
                         double  Tfoliage,
                         double  SWdown,
                         double *aPAR,
                         double  elevation,
                         double  Catm,
                         double *CanopLayerBnd,
                         double  LAItotal,
                         char   *mode,
                         double *rsLayer,
                         double *rc,
                         double *Ci,
                         double *GPP,
                         double *Rdark,
                         double *Rphoto,
                         double *Rmaint,
                         double *Rgrowth,
                         double *Raut,
                         double *NPP)
{
  extern option_struct options;
  double  h;
  double  pz;
  int     cidx;
  double  dLAI;
  double *CiLayer;
  double  AgrossLayer;
  double  RdarkLayer;
  double  RphotoLayer;
  double  gc;                  /* 1/rs */

  /* calculate scale height based on average temperature in the column */
  h  = 287/9.81 * ((Tfoliage + 273.15) + 0.5 * (double)elevation * T_LAPSE);

  /* use hypsometric equation to calculate p_z, assume that virtual
     temperature is equal air_temp */
  pz = PS_PM * exp(-(double)elevation/h);

  CiLayer = (double*)calloc(options.Ncanopy,sizeof(double));

  if (!strcasecmp(mode,"ci")) {

    /* Assume a default leaf-internal CO2; compute assimilation, respiration, and stomatal resistance */

    /* Default leaf-internal CO2 */
    for (cidx = 0; cidx < options.Ncanopy; cidx++) {
      if (Ctype == PHOTO_C3)
        CiLayer[cidx] = FCI1C3 * Catm;
      else if (Ctype == PHOTO_C4)
        CiLayer[cidx] = FCI1C4 * Catm;
    }
    if (Ctype == PHOTO_C3)
      *Ci = FCI1C3 * Catm;
    else if (Ctype == PHOTO_C4)
      *Ci = FCI1C4 * Catm;

    /* Sum over canopy layers */
    *GPP    = 0.0;
    *Rdark  = 0.0;
    *Rphoto = 0.0;
    gc      = 0.0;
    for (cidx = 0; cidx < options.Ncanopy; cidx++) {

      photosynth(Ctype,
                 MaxCarboxRate,
                 MaxETransport,
                 CO2Specificity,
                 NscaleFactor[cidx],
                 Tfoliage,
                 SWdown/Epar, /* note: divide by Epar to convert from W/m2 to mol(photons)/m2s */
                 aPAR[cidx],
                 pz,
                 Catm,
                 mode,
                 &(rsLayer[cidx]),
                 &(CiLayer[cidx]),
                 &RdarkLayer,
                 &RphotoLayer,
                 &AgrossLayer);

      if (cidx > 0)
        dLAI = LAItotal * (CanopLayerBnd[cidx] - CanopLayerBnd[cidx-1]);
      else
        dLAI = LAItotal * CanopLayerBnd[cidx];

      *GPP    += AgrossLayer * dLAI;
      *Rdark  += RdarkLayer * dLAI;
      *Rphoto += RphotoLayer * dLAI;
      gc      += (1/rsLayer[cidx]) * dLAI;

    }

    if (gc < SMALL) gc = SMALL;
    *rc = 1/gc;
    if (*rc > HUGE_RESIST) *rc = HUGE_RESIST;

  }
  else {

    /* Stomatal resistance given; compute assimilation, respiration, and leaf-internal CO2 */

    /* Sum over canopy layers */
    *GPP    = 0.0;
    *Rdark  = 0.0;
    *Rphoto = 0.0;
    *Ci     = 0.0;
    for (cidx = 0; cidx < options.Ncanopy; cidx++) {

      photosynth(Ctype,
                 MaxCarboxRate,
                 MaxETransport,
                 CO2Specificity,
                 NscaleFactor[cidx],
                 Tfoliage,
                 SWdown/Epar,
                 aPAR[cidx],
                 pz,
                 Catm,
                 mode,
                 &(rsLayer[cidx]),
                 &(CiLayer[cidx]),
                 &RdarkLayer,
                 &RphotoLayer,
                 &AgrossLayer);

      if (cidx > 0)
        dLAI = LAItotal * (CanopLayerBnd[cidx] - CanopLayerBnd[cidx-1]);
      else
        dLAI = LAItotal * CanopLayerBnd[cidx];

      *GPP    += AgrossLayer * dLAI;
      *Rdark  += RdarkLayer * dLAI;
      *Rphoto += RphotoLayer * dLAI;
      *Ci     += CiLayer[cidx] * dLAI;

    }

  }

  /* Compute whole-plant respiration terms and NPP */
  *Rmaint = *Rdark/FRLeaf;
  *Rgrowth = (FRGrowth/(1+FRGrowth))*((*GPP)-(*Rmaint));
//  if (*Rgrowth < 0) *Rgrowth = 0;
  *Raut = *Rmaint + *Rgrowth;
  *NPP = *GPP - *Raut;

  free((char*)CiLayer);

}

