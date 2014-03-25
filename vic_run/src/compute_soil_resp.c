/********************************************************************************
  filename  : compute_soil_resp.c
  purpose   : Calculate soil respiration (heterotrophic respiration, or Rh)
  interface : - input :
                - Nnodes:  number of points (vertically) at which to evaluate soil respiration
                - dZ[]:    array of layer thicknesses [m] pertaining to nodes
                - dZTot:   total depth [m] of all layers (total depth of soil containing carbon)
                - dt:      timestep length [h]
                - T[]:     array of node temperatures [K]
                - w[]:     array of node moisture fractions [fraction]
                - CLitter: carbon storage in litter pool [gC/m2]
                - CInter:  carbon storage in intermediate pool [gC/m2]
                - CSlow:   carbon storage in slow pool [gC/m2]
              Note: the litter pool is concentrated at the soil surface, while the intermediate and
                    slow pools are both assumed to be distributed uniformly throughout soil column.

              - constants:
                - E0_LT:     Lloyd-Taylor E0 parameter [K]
                - T0_LT:     Lloyd-Taylor E0 parameter [K]
                - wminFM:    minimum soil moisture (fraction) at which soil respiration can occur
                - wmaxFM:    maximum soil moisture (fraction) at which soil respiration can occur
                - woptFM:    soil moisture (fraction) at which maximum soil respiration occurs
                - Rhsat:     ratio of soil respiration rate under saturated conditions (w=wmaxFM) to
                             that under optimal conditions (w=woptFM)
                - Rfactor:   scaling factor to account for other (non-moisture) sources of inhibition
                             of respiration
                - tauLitter: Litter pool turnover time [y]
                - tauInter:  Intermediate pool turnover time [y]
                - tauSlow:   Slow pool turnover time [y]

              - output:
                - RhLitter:     Soil respiration from litter pool [gC/m2]
                - RhInterTot:   Soil respiration from intermediate pool [gC/m2]
                - RhSlowTot:    Soil respiration from slow pool [gC/m2]

  programmer: Ted Bohn
  date      : July 25, 2013
  changes   :
  references: 
********************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vicNl.h>

static char vcid[] = "$Id: $";

void compute_soil_resp(int Nnodes,
                       double *dZ,
                       double dZTot,
                       double dt,
                       double *T,
                       double *w,
                       double CLitter,
                       double CInter,
                       double CSlow,
                       double *RhLitter,
                       double *RhInterTot,
                       double *RhSlowTot)
{
  extern option_struct options;
  int i;
  double Tref;
  double *TK;
  double fTLitter;
  double *fTSoil;
  double fMLitter;
  double *fMSoil;
  double *CInterNode;
  double *CSlowNode;
  double *RhInter;
  double *RhSlow;

  /* Allocate temp arrays */
  TK = (double*)calloc(Nnodes,sizeof(double));
  fTSoil = (double*)calloc(Nnodes,sizeof(double));
  fMSoil = (double*)calloc(Nnodes,sizeof(double));
  CInterNode = (double*)calloc(Nnodes,sizeof(double));
  CSlowNode = (double*)calloc(Nnodes,sizeof(double));
  RhInter = (double*)calloc(Nnodes,sizeof(double));
  RhSlow = (double*)calloc(Nnodes,sizeof(double));

  /* Compute Lloyd-Taylor temperature dependence */
  Tref = 10+KELVIN; /* reference temperature of 10 C */
  for (i=0; i<Nnodes; i++) {
    TK[i] = T[i]+KELVIN;
    if (TK[i] < T0_LT) TK[i] = T0_LT;
  }
  fTLitter = exp( E0_LT * ( 1/(Tref-T0_LT) - 1/(TK[0]-T0_LT) ) );
  for (i=0; i<Nnodes; i++) {
    fTSoil[i] = exp( E0_LT * ( 1/(Tref-T0_LT) - 1/(TK[i]-T0_LT) ) );
  }

  /* Compute moisture dependence */
  for (i=0; i<Nnodes; i++) {
    if (w[i] < wminFM) w[i] = wminFM;
    if (w[i] > wmaxFM) w[i] = wmaxFM;
  }
  if (w[0] <= woptFM)
   fMLitter = (w[0]-wminFM)*(w[0]-wmaxFM)/((w[0]-wminFM)*(w[0]-wmaxFM)-(w[0]-woptFM)*(w[0]-woptFM));
  else
   fMLitter = Rhsat + (1-Rhsat)*(w[0]-wminFM)*(w[0]-wmaxFM)/((w[0]-wminFM)*(w[0]-wmaxFM)-(w[0]-woptFM)*(w[0]-woptFM));
  if (fMLitter > 1.0) fMLitter = 1.0;
  if (fMLitter < 0.0) fMLitter = 0.0;
  for (i=0; i<Nnodes; i++) {
    if (w[i] <= woptFM)
     fMSoil[i] = (w[i]-wminFM)*(w[i]-wmaxFM)/((w[i]-wminFM)*(w[i]-wmaxFM)-(w[i]-woptFM)*(w[i]-woptFM));
    else
     fMSoil[i] = Rhsat + (1-Rhsat)*(w[i]-wminFM)*(w[i]-wmaxFM)/((w[i]-wminFM)*(w[i]-wmaxFM)-(w[i]-woptFM)*(w[i]-woptFM));
    if (fMSoil[i] > 1.0) fMSoil[i] = 1.0;
    if (fMSoil[i] < 0.0) fMSoil[i] = 0.0;
  }

  /* Compute C per node */
  for (i=0; i<Nnodes; i++) {
    CInterNode[i] = CInter * dZ[i]/dZTot;
    CSlowNode[i] = CSlow * dZ[i]/dZTot;
  }

  /* Compute Rh for various pools, nodes; C fluxes in [gC/m2d] */
  *RhLitter = Rfactor*(fTLitter*fMLitter/(tauLitter*365.25*24/dt))*CLitter;
  *RhInterTot = 0;
  *RhSlowTot = 0;
  for (i=0; i<Nnodes; i++) {
    RhInter[i] = Rfactor*(fTSoil[i]*fMSoil[i]/(tauInter*365.25*24/dt))*CInterNode[i];
    RhSlow[i] = Rfactor*(fTSoil[i]*fMSoil[i]/(tauSlow*365.25*24/dt))*CSlowNode[i];
    *RhInterTot += RhInter[i];
    *RhSlowTot += RhSlow[i];
  }

  free((char *)TK);
  free((char *)fTSoil);
  free((char *)fMSoil);
  free((char *)CInterNode);
  free((char *)CSlowNode);
  free((char *)RhInter);
  free((char *)RhSlow);

}

