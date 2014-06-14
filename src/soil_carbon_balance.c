/********************************************************************************
  filename  : soil_carbon_balance.c
  purpose   : Calculate soil carbon balance
  interface : - input :
#                - Nnodes:  number of points (vertically) at which to evaluate soil respiration
#                - dZ[]:    array of layer thicknesses [m] pertaining to nodes
#                - dZTot:   total depth [m] of all layers (total depth of soil containing carbon)
#                - T[]:     array of node temperatures [K]
#                - w[]:     array of node moisture fractions [fraction]
#                - CLitter: carbon storage in litter pool [gC/m2]
#                - CInter:  carbon storage in intermediate pool [gC/m2]
#                - CSlow:   carbon storage in slow pool [gC/m2]
#              Note: the litter pool is concentrated at the soil surface, while the intermediate and
#                    slow pools are both assumed to be distributed uniformly throughout soil column.
#
#              - constants:
#                - fAir:      Fraction of respired carbon from litter pool that is lost to atmosphere
#                - fInter:    Fraction of [respired carbon from litter pool that goes to soil] that goes
#                             to intermediate pool
#
#              - output:
#                - RhLitter:     Soil respiration from litter pool [gC/m2]
#                - RhLitter2Atm: Soil respiration from litter pool that is lost to atmosphere [gC/m2]
#                - RhInter:      Soil respiration from intermediate pool [gC/m2]
#                - RhSlow:       Soil respiration from slow pool [gC/m2]
#                - RhTot:        Total soil respiration from all pools [gC/m2]
#                - Litterfall:   Flux of carbon from living biomass into soil [gC/m2]
#                - CLitter:      Updated carbon storage in litter pool [gC/m2]
#                - CInter:       Updated carbon storage in intermediate pool [gC/m2]
#                - CSlow:        Updated carbon storage in slow pool [gC/m2]
#              Note: We assume that respiration in the litter pool results in a mix of end products: CO2,
#                    which is lost to the atmosphere, and larger molecules, which become part of the soil
#                    carbon pools (intermediate and slow).  The fraction of respired carbon (by mass) that
#                    ends up as CO2 is given by fAir.  All respiration from the intermediate and slow pools
#                    is assumed to go to the atmosphere.
                    

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

void soil_carbon_balance(soil_con_struct *soil_con,
                         energy_bal_struct *energy,
                         cell_data_struct *cell,
                         veg_var_struct  *veg_var)
{
  extern option_struct options;
  extern global_param_struct global_param;
  int i;
  int lidx;
  int Nnodes;
  double *dZ;
  double *dZCum;
  double dZTot;
  double *T;
  double *w;
  double tmp_double;
  double b;
  double wtd;
  double w0;
  double w1;

  // Find subset of thermal nodes that span soil hydrologic layers
  dZTot = 0;
  for (i=0; i<options.Nlayer; i++) {
    dZTot += soil_con->depth[i];
  }
  i=0;
  while (i<options.Nnode-1 && soil_con->Zsum_node[i]<dZTot) {
    i++;
  }
  Nnodes = i;
  if (soil_con->Zsum_node[i] > dZTot) {
    Nnodes--;
  }
  dZ = (double*)calloc(Nnodes,sizeof(double));
  dZCum = (double*)calloc(Nnodes,sizeof(double));
  T = (double*)calloc(Nnodes,sizeof(double));
  w = (double*)calloc(Nnodes,sizeof(double));

  // Assign node thicknesses and temperatures for subset
  dZTot = 0;
  for (i=0; i<Nnodes; i++) {
    dZ[i] = soil_con->dz_node[i]*1000; // mm
    dZTot += dZ[i];
    dZCum[i] = dZTot;
    T[i] = energy->T[i];
  }

  // Compute node relative moistures based on lumped water table depth
  for (i=0; i<Nnodes; i++) {
    b = 0.5*(soil_con->expt_node[i]-3);
    wtd = -(cell->zwt_lumped)*10; // mm, positive downwards
    if (wtd > dZCum[i]) {
      if (i>0)
        w0 = pow((wtd+soil_con->bubble_node[i]-dZCum[i-1])/soil_con->bubble_node[i], -1/b);
      else
        w0 = pow((wtd+soil_con->bubble_node[i])/soil_con->bubble_node[i], -1/b);
      w1 = pow((wtd+soil_con->bubble_node[i]-dZCum[i])/soil_con->bubble_node[i], -1/b);
      w[i] = 0.5*(w0+w1);
    }
    else if ((i==0 && wtd>0) || (i>0 && wtd>dZCum[i-1]) ) {
      if (i>0) {
        w0 = pow((wtd+soil_con->bubble_node[i]-dZCum[i-1])/soil_con->bubble_node[i], -1/b);
        tmp_double = 0.5*(dZCum[i-1]+wtd);
        w1 = pow((wtd+soil_con->bubble_node[i]-tmp_double)/soil_con->bubble_node[i], -1/b);
        w[i] = ( 0.5*(w0+w1)*(tmp_double-dZCum[i-1]) + 0.5*(w1+1.0)*(wtd-tmp_double) + 1.0*(dZCum[i]-wtd) ) / (dZCum[i]-dZCum[i-1]);
      }
      else {
        w0 = pow((wtd+soil_con->bubble_node[i])/soil_con->bubble_node[i], -1/b);
        tmp_double = 0.5*(0+wtd);
        w1 = pow((wtd+soil_con->bubble_node[i]-tmp_double)/soil_con->bubble_node[i], -1/b);
        w[i] = ( 0.5*(w0+w1)*(tmp_double-0) + 0.5*(w1+1.0)*(wtd-tmp_double) + 1.0*(dZCum[i]-wtd) ) / (dZCum[i]-0);
      }
    }
    else {
      w[i] = 1.0;
    }
  }

  // Compute carbon fluxes out of soil (evaluate fluxes at subset of thermal nodes and sum them)
  compute_soil_resp(Nnodes,dZ,dZTot,global_param.dt,T,w,cell->CLitter,cell->CInter,cell->CSlow,&(cell->RhLitter),&(cell->RhInter),&(cell->RhSlow));
  cell->RhLitter2Atm = fAir*cell->RhLitter;
  cell->RhTot = cell->RhLitter2Atm + cell->RhInter + cell->RhSlow;

  // Compute balances of soil carbon pools
  veg_var->Litterfall = veg_var->AnnualNPPPrev/(365.25*24/global_param.dt); // Assume previous year's NPP enters soil evenly throughout current year
  cell->CLitter += veg_var->Litterfall - cell->RhLitter;
  cell->CInter += (1-fAir)*cell->RhLitter*fInter - cell->RhInter;
  cell->CSlow += (1-fAir)*cell->RhLitter*(1-fInter) - cell->RhSlow;

  // Free temporary dynamic arrays
  free((char *)dZ);
  free((char *)dZCum);
  free((char *)T);
  free((char *)w);

}

