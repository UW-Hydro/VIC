#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

static char vcid[] = "$Id$";

double frozen_soil_conductivity(double moist, 
                                double Wu, 
                                double soil_density, 
                                double bulk_density,
                                double quartz) {
/**********************************************************************
  Soil thermal conductivity calculated using Johansen's method.

  Reference: Farouki, O.T., "Thermal Properties of Soils" 1986
	Chapter 7: Methods for Calculating the Thermal Conductivity 
		of Soils

  H.B.H. - refers to the handbook of hydrology.

  porosity = n = porosity
  ratio = Sr = fractionaldegree of saturation
  All K values are conductivity in W/mK
  Wu is the fractional volume of unfrozen water

  The follow was extracted from Xu Liang's code:
    Ke = 0.7*log10(Sr)+1.0; * for coarse unfrozen soil *
    Ke = log10(Sr)+1.0;     * for fine unfrozen soil *
    Ksolid = 2.289; * for silty clay with 10% sand, in unit of W/mK *

  UNITS: input in m, kg, s

  Returns K in W/m/K

**********************************************************************/
  double Ke;
  double Ki = 2.2;	/* thermal conductivity of ice (W/mK) */
  double Kw = 0.57;	/* thermal conductivity of water (W/mK) */
  double Ksat;
  double Ks;		/* thermal conductivity of solid (W/mK)
			   function of quartz content */
  double Kdry;
  double Sr;		/* fractional degree of saturation */
  double K;
  double porosity;

  if(moist>0.) {
    porosity = 1.0 - bulk_density / soil_density;

    Sr = moist/porosity;

    Ks = pow(7.7,quartz) * pow(2.2,1.0-quartz);

    Ksat = pow(Ks,1.0-porosity) * pow(Ki,porosity-Wu) * pow(Kw,Wu);

    Ke = Sr;
    Kdry = (0.135*bulk_density+64.7)/(soil_density-0.947*bulk_density);
    K = (Ksat-Kdry)*Ke+Kdry;
    if(K<Kdry) K=Kdry;
  }
  else K=0.;

  return (K); 
}


double soil_conductivity(double moist, 
                         double soil_density, 
                         double bulk_density,
                         double quartz) {
/**********************************************************************
	soil_conductivity	

  This routine calculates the conductivity of unfrozen soils

  UNITS: returns W/m/K

**********************************************************************/
  double Ke;
  double Kw = 0.57;	/* thermal conductivity of water (W/mK) */
  double Ksat;
  double Ks;		/* thermal conductivity of solid (W/mK)
			   function of quartz content */
  double Kdry;
  double Sr;		/* fractional degree of saturation */
  double K;
  double porosity;

  if(moist>0.) {
    porosity = 1.0 - bulk_density / soil_density;

    Sr = moist/porosity;

    Ks = pow(7.7,quartz) * pow(2.2,1.0-quartz);

    Ksat = pow(Ks,1.0-porosity) * pow(Kw,porosity);

    Ke = 0.7 * log10(Sr) + 1.0;
    Kdry = (0.135*bulk_density+64.7)/(soil_density-0.947*bulk_density);
    K = (Ksat-Kdry)*Ke+Kdry;
    if(K<Kdry) K=Kdry;
  }
  else K=0.;

  return (K);

}

#define organic_fract 0.00

double volumetric_heat_capacity(double soil_fract,
                                double water_fract,
                                double ice_fract) {
/**********************************************************************
  This subroutine calculates the soil volumetric heat capacity based 
  on the fractional volume of its component parts.

  Constant values are volumetric heat capacities in J/m^3/K
	Soil value is for clay or quartz - assumed for all other types
**********************************************************************/

  double Cs;

  Cs = 2.0e6 * (soil_fract - organic_fract);
  Cs += 4.2e6 * water_fract;
  Cs += 1.9e6 * ice_fract;
  Cs += 2.7e6 * organic_fract;
  Cs += 1.3e3 * (1. - (soil_fract + water_fract + ice_fract + organic_fract));

  return (Cs);

}

#undef organic_fract

void distribute_soil_property(double  *dz,
                              double   fdepth,
                              double   tdepth,
                              double **l_param,
                              int      Nlayer,
                              int      Nnodes,
                              double  *depth,
                              double  *param) {
/**********************************************************************
  This subroutine distributes soil parameters calculated for freezing,
  thawing, upper, and lower soil layers, for all layers used by the
  frozen soil model.

  Note: this subroutine assumes that the values of dz are layer 
  thicknesses, centered at the temperature solution depth, therefore
  the upper and lower most values of dz, are twice their actual
  thickness since they are at the upper and lower bounds of the
  soil column.
**********************************************************************/

  int    lindex, l2index, zindex;
  double Zsum, Lsum;
  double Ltmp, Ftmp, Ttmp, Ztmp;

  dz[0] /= 2.;
  while(dz[Nnodes-1]<=0) Nnodes--;
  dz[Nnodes-1] /= 2.;

  Zsum=0.;
  Lsum = depth[0];
  lindex=l2index=0;

  for(zindex=0;zindex<Nnodes;zindex++) {
    Ltmp = (double)((int)(((Lsum - Zsum)+0.0005)*1000.))/1000.;
    Ftmp = (double)((int)(((fdepth - Zsum)+0.0005)*1000.))/1000.;
    Ttmp = (double)((int)(((tdepth - Zsum)+0.0005)*1000.))/1000.;
    Ztmp = (double)((int)(((dz[zindex])+0.0005)*1000.))/1000.;
    if(lindex == Nlayer-1 && Ztmp>Ltmp) Ztmp = Ltmp;

    if(Ttmp>=Ztmp && Ltmp>=Ztmp)
      param[zindex] = l_param[lindex][0];
    else if(Ttmp<0 && Ftmp>=Ztmp && Ltmp>=Ztmp)
      param[zindex] = l_param[lindex][1];
    else if(Ttmp<0 && Ftmp<0 && Ltmp>=Ztmp)
      param[zindex] = l_param[lindex][2];
    else if(Ttmp>=0 && Ttmp<Ztmp && Ftmp>=Ztmp && Ltmp>=Ztmp)
      param[zindex] = (l_param[lindex][0]*Ttmp 
          + l_param[lindex][1]*(dz[zindex]-Ttmp)) / dz[zindex];
    else if(Ttmp<0 && Ftmp>=0 && Ftmp<Ztmp && Ltmp>=Ztmp)
      param[zindex] = (l_param[lindex][1]*Ftmp 
          + l_param[lindex][2]*(dz[zindex]-Ftmp)) / dz[zindex];
    else if(Ttmp<0 && Ftmp<0 && Ltmp>=0 && Ltmp<Ztmp)
      param[zindex] = (l_param[lindex][2]*Ltmp 
		       + l_param[lindex+1][2]*(dz[zindex]-Ltmp)) / dz[zindex];
    else if(Ttmp<0 && Ftmp>=Ztmp && Ltmp>=0 && Ltmp<Ztmp)
      param[zindex] = (l_param[lindex][1]*Ltmp 
		       + l_param[lindex+1][1]*(dz[zindex]-Ltmp)) / dz[zindex];
    else if(Ttmp>=Ztmp && Ftmp>=Ztmp && Ltmp>=0 && Ltmp<Ztmp)
      param[zindex] = (l_param[lindex][0]*Ltmp 
		       + l_param[lindex+1][0]*(dz[zindex]-Ltmp)) / dz[zindex];
    else if(Ttmp>=0 && Ftmp>=Ztmp && Ltmp>=Ttmp && Ltmp<Ztmp)
      param[zindex] = (l_param[lindex][0]*Ttmp 
		       + l_param[lindex][1]*(Ltmp-Ttmp) 
		       + l_param[lindex+1][1]*(dz[zindex]-Ltmp)) / dz[zindex];
    else if(Ttmp>=Ltmp && Ttmp<Ztmp && Ftmp>=Ztmp && Ltmp>=0)
      param[zindex] = (l_param[lindex][0]*Ltmp 
		       + l_param[lindex+1][0]*(Ttmp-Ltmp) 
		       + l_param[lindex+1][1]*(dz[zindex]-Ttmp)) / dz[zindex];
    else if(Ttmp<0 && Ftmp>=0 && Ltmp>=Ftmp && Ltmp<Ztmp)
      param[zindex] = (l_param[lindex][1]*Ftmp 
		       + l_param[lindex][2]*(Ltmp-Ftmp) 
		       + l_param[lindex+1][2]*(dz[zindex]-Ltmp)) / dz[zindex];
    else if(Ttmp<0 && Ftmp>=Ltmp && Ftmp<Ztmp && Ltmp>=0)
      param[zindex] = (l_param[lindex][1]*Ltmp 
		       + l_param[lindex+1][1]*(Ftmp-Ltmp) 
		       + l_param[lindex+1][2]*(dz[zindex]-Ftmp)) / dz[zindex];
    else if(Ttmp>=0 && Ftmp>=Ttmp && Ftmp<Ztmp && Ltmp>=Ztmp)
      param[zindex] = (l_param[lindex][0]*Ttmp 
		       + l_param[lindex][1]*(Ftmp-Ttmp) 
		       + l_param[lindex][2]*(dz[zindex]-Ftmp)) / dz[zindex];
    else if(Ttmp>=0 && Ftmp>=Ttmp && Ltmp>=Ftmp && Ltmp<Ztmp)
      param[zindex] = (l_param[lindex][0]*Ttmp 
		       + l_param[lindex][1]*(Ftmp-Ttmp) 
		       + l_param[lindex][2]*(Ltmp-Ftmp) 
		       + l_param[lindex+1][2]*(dz[zindex]-Ltmp)) / dz[zindex];
    else if(Ttmp>=0 && Ftmp>=Ltmp && Ftmp<Ztmp && Ltmp>=Ttmp)
      param[zindex] = (l_param[lindex][0]*Ttmp 
		       + l_param[lindex][1]*(Ltmp-Ttmp) 
		       + l_param[lindex+1][1]*(Ftmp-Ltmp) 
		       + l_param[lindex+1][2]*(dz[zindex]-Ftmp)) / dz[zindex];
    else if(Ttmp>=Ltmp && Ftmp>=Ttmp && Ftmp<Ztmp && Ltmp>=0)
      param[zindex] = (l_param[lindex][0]*Ltmp 
		       + l_param[lindex+1][0]*(Ttmp-Ltmp) 
		       + l_param[lindex+1][1]*(Ftmp-Ttmp) 
		       + l_param[lindex+1][2]*(dz[zindex]-Ftmp)) / dz[zindex];

    if(Ltmp>=0 && Ltmp<Ztmp) {
      lindex++;
      Lsum+=depth[lindex];
    }
    Zsum += dz[zindex];
  }

  dz[0] *= 2.;
  dz[Nnodes-1] *= 2.;

}

void soil_thermal_calc(soil_con_struct    soil_con,
                       layer_data_struct *layer,
                       energy_bal_struct  energy,
                       double            *kappa,
                       double            *Cs, 
                       double            *moist,
                       int                Nlayer,
                       int                Nnodes) {
/**********************************************************************
  This subroutine will calculate thermal conductivity and volumetric
  heat capacity for frozen, thawed, and unfrozen sublayers, of each
  soil layer.

  kappa and Cs should be allocated before calling this routine.
**********************************************************************/

  extern option_struct options;
  extern debug_struct debug;

  int      zindex, index;
  double   Lsum, Zsum;
  double **K_layer;
  double **Cs_layer;
  double **moist_layer;
  double  *depth_mm;
  double   unfrozen;

  depth_mm = (double *)calloc(options.Nlayer,sizeof(double));
  for(index=0;index<options.Nlayer;index++) 
    depth_mm[index] = soil_con.depth[index] * 1000.;

  if(kappa != NULL && Cs != NULL && options.FROZEN_SOIL) {

    moist_layer = (double **)calloc(options.Nlayer,sizeof(double*));
    for(index=0;index<options.Nlayer;index++) {
      moist_layer[index]     = (double *)calloc(3,sizeof(double));
      moist_layer[index][0]  = layer[index].moist_thaw / depth_mm[index];
      moist_layer[index][1]  = layer[index].moist_froz / depth_mm[index];
      moist_layer[index][1] += layer[index].ice / depth_mm[index];
      moist_layer[index][2]  = layer[index].moist / depth_mm[index];
    }
    distribute_soil_property(energy.dz,energy.fdepth[0],energy.fdepth[1],
        moist_layer,options.Nlayer,Nnodes,soil_con.depth,moist);
    for(index=0;index<options.Nlayer;index++) free((char*)moist_layer[index]);
    free((char*)moist_layer);

    Lsum = 0.;
    Zsum = 0.;
    index = 0;
    for(zindex=0;zindex<Nnodes;zindex++) {

      if(energy.T[zindex]<0.) {
        unfrozen = maximum_unfrozen_water(energy.T[zindex],
					  soil_con.max_moist_node[zindex], 
					  soil_con.bubble, 
					  soil_con.expt_node[zindex]);
        if(unfrozen>soil_con.max_moist_node[zindex] || unfrozen<0.)
          unfrozen = soil_con.max_moist_node[zindex];
        if(unfrozen > moist[zindex]) unfrozen = moist[zindex];
        kappa[zindex] = frozen_soil_conductivity(moist[zindex],unfrozen,
            soil_con.soil_density,soil_con.bulk_density[index],
            soil_con.quartz);
        Cs[zindex] = volumetric_heat_capacity(soil_con.bulk_density[index]
            / soil_con.soil_density, unfrozen, moist[zindex]-unfrozen);
      }
      else {
        kappa[zindex] = soil_conductivity(moist[zindex],
					  soil_con.soil_density,
					  soil_con.bulk_density[index],
					  soil_con.quartz);
        Cs[zindex] = volumetric_heat_capacity(soil_con.bulk_density[index]
            / soil_con.soil_density, moist[zindex], 0.);
      }
  
      if(zindex<Nnodes-1) {
        Zsum += (energy.dz[zindex] + energy.dz[zindex+1]) / 2.;
        if(((int)(((Zsum - Lsum)+0.0005)*1000.) > (int)(depth_mm[index]))
           && index<options.Nlayer-1) {
          Lsum += soil_con.depth[index];
          index++;
        }
      }

    }

    if(kappa[0]==0) {
      kappa[0]=kappa[1];
      Cs[0] = Cs[1];
    }

  }

  K_layer = (double **)calloc(Nlayer,sizeof(double *));
  Cs_layer = (double **)calloc(Nlayer,sizeof(double *));

  for(index=0;index<Nlayer;index++) {
  
    K_layer[index] = (double *)calloc(Nlayer,sizeof(double));
    Cs_layer[index] = (double *)calloc(Nlayer,sizeof(double));

    /** Compute Thermal Conductivity of VIC Layers **/
    K_layer[index][2] = soil_conductivity(layer[index].moist/depth_mm[index],
        soil_con.soil_density,soil_con.bulk_density[index],soil_con.quartz);
    if(layer[index].tdepth > 0.0)
      K_layer[index][0] = soil_conductivity(layer[index].moist_thaw
          /depth_mm[index],soil_con.soil_density,
          soil_con.bulk_density[index],soil_con.quartz);
    else K_layer[index][0] = K_layer[index][2];
    if(layer[index].fdepth > 0.0)
      K_layer[index][1] = frozen_soil_conductivity((layer[index].moist_froz
          + layer[index].ice) / depth_mm[index],
          layer[index].moist_froz/depth_mm[index],
          soil_con.soil_density,soil_con.bulk_density[index],soil_con.quartz);
    else K_layer[index][1] = K_layer[index][2];

    /** Compute Volumetric Heat Capacities of VIC Layers **/
    Cs_layer[index][2] = volumetric_heat_capacity(soil_con.bulk_density[index]
        / soil_con.soil_density, layer[index].moist/depth_mm[index],0.0);
    if(layer[index].tdepth > 0.0)
      Cs_layer[index][0] = volumetric_heat_capacity(soil_con.bulk_density[index]
          / soil_con.soil_density,layer[index].moist_thaw/depth_mm[index],0.0);
    else Cs_layer[index][0] = Cs_layer[index][2];
    if(layer[index].fdepth > 0.0)
      Cs_layer[index][1] 
	= volumetric_heat_capacity(soil_con.bulk_density[index]
				   / soil_con.soil_density,
				   layer[index].moist_froz/depth_mm[index],
				   layer[index].ice/depth_mm[index]);
    else Cs_layer[index][1] = Cs_layer[index][2];
  
  } 

  /** Compute Average Thermal Conductivities and Volumetric
      Heat Capacities for Model Layers **/

  for(index=0;index<Nlayer;index++) {
    layer[index].Cs = Cs_layer[index][0] * layer[index].tdepth;
    layer[index].Cs += Cs_layer[index][1] * (layer[index].fdepth 
        - layer[index].tdepth);
    layer[index].Cs += Cs_layer[index][2] * (soil_con.depth[index]
        - layer[index].fdepth);
    layer[index].Cs /= soil_con.depth[index];
    layer[index].kappa = K_layer[index][0] * layer[index].tdepth;
    layer[index].kappa += K_layer[index][1] * (layer[index].fdepth 
        - layer[index].tdepth);
    layer[index].kappa += K_layer[index][2] * (soil_con.depth[index] 
        - layer[index].fdepth);
    layer[index].kappa /= soil_con.depth[index];
    free((char*)K_layer[index]);
    free((char*)Cs_layer[index]);
  }
  free((char*)depth_mm);
  free((char*)K_layer);
  free((char*)Cs_layer);

}

double maximum_unfrozen_water(double T,
                              double max_moist,
                              double bubble,
                              double expt) {
/**********************************************************************
  This subroutine computes the maximum amount of unfrozen water that
  can exist at the current temperature.
**********************************************************************/

  double unfrozen;

  unfrozen = max_moist * pow((-Lf * T) / (T
      + 273.16) / (9.18 * bubble / 100.), -(2.0 / (expt - 3.0)));

  return (unfrozen);
}

void find_0_degree_fronts(energy_bal_struct *energy,
                          layer_data_struct *layer,
                          double dp,
                          double *depth,
                          double *T,
                          int Nlayer,
                          int Nnodes) {
/**********************************************************************
  This subroutine finds the freezing and thawing fronts, and locates
  them within the fixed soil layers.
**********************************************************************/

  int    index;
  double Zsum, Lsum;
  double tdepth, fdepth;
  double MINLAYER = 0.001;

  /** Calculate New Layer Depths **/
  Lsum = 0.;
  for(index=0;index<Nlayer;index++) Lsum += depth[index];
  if(dp > Lsum) {
    Nnodes--;
    Zsum = Lsum;
  }
  else Zsum=dp;
  index=Nnodes-1;
  fdepth=tdepth=0;
  if(T[index]>0.0) {
    /***** Find Freezing Front Depth *****/
    while(index>0 && T[index]>0.0) {
      index--;
      Zsum -= (energy->dz[index+1]+energy->dz[index])/2.0;
    }
    if(index==0) Zsum=0.0;

    if(index>=0 && T[index]<0.0) {
      fdepth = linear_interp(0.0,T[index],T[index+1],
          Zsum,Zsum+(energy->dz[index+1]+energy->dz[index])/2.0);
      /***** Find Thawing Front Depth *****/
      while(index>0 && T[index]<=0.0) {
        index--;
        Zsum -= (energy->dz[index+1]+energy->dz[index])/2.0;
      }
      if(index==0) Zsum=0.0;

      if(index>=0 && T[index]>0.0) tdepth = linear_interp(0.0,
          T[index],T[index+1],Zsum,Zsum+(energy->dz[index+1]
          +energy->dz[index])/2.0);
      else tdepth=0.0;
    }
    else {
      fdepth=0.0;
      tdepth=0.0;
    }
  }
  else {
    /***** Freezing Front Penetrates Lower Layer *****/
    fprintf(stderr,"WARNING: Frozen layer extends through entire lower layer.\n");
    fdepth=Zsum;
  }

  if(fdepth-tdepth < MINLAYER) fdepth = tdepth = 0.;
  energy->fdepth[0] = fdepth;
  energy->fdepth[1] = tdepth;

  Lsum = 0.;
  for(index=0;index<Nlayer;index++) {
    if(tdepth > Lsum && tdepth < Lsum + depth[index])
      layer[index].tdepth = tdepth - Lsum;
    else if(tdepth<=Lsum) layer[index].tdepth = 0.;
    else layer[index].tdepth = depth[index];
    if(fdepth > Lsum && fdepth < Lsum + depth[index])
      layer[index].fdepth = fdepth - Lsum;
    else if(fdepth<=Lsum) layer[index].fdepth = 0.;
    else layer[index].fdepth = depth[index];
    Lsum += depth[index];
  }
 
}

void redistribute_moisture(layer_data_struct *layer,
                           double *fdepth,
                           double *max_moist,
                           double *old_fdepth,
                           double *old_tdepth,
                           double *depth,
                           int Nlayer) {
/**********************************************************************
  This subroutine redistributes soil moisture between from the upper,
  frozen and thawed layers of the previous time step, to those
  computed for the current time step.
 
  Ice is melted and redistributed as liquid in this subroutine, it
  must be refrozen in a later step.
**********************************************************************/
 
  int index;
  double temp_moist;
  double Zsum;
  double last_td, last_fd;
  double next_td, next_fd;
 
  Zsum=0.;

  for(index=0;index<Nlayer;index++) {

    layer[index].moist_froz += layer[index].ice;
    layer[index].ice = 0.;
    layer[index].moist_thaw *= (old_tdepth[index]) / depth[index];
    layer[index].moist_froz *= (old_fdepth[index] - old_tdepth[index]) 
        / depth[index];
    layer[index].moist *= (depth[index] - old_fdepth[index]) / depth[index];

    last_td = old_tdepth[index];
    if(last_td < 0) last_td=0.;
    else if(last_td > depth[index]) last_td=depth[index];

    next_td = fdepth[1] - Zsum;
    if(next_td < 0) next_td=0.;
    else if(next_td > depth[index]) next_td=depth[index];

    last_fd = old_fdepth[index];
    if(last_fd < 0) last_fd=0.;
    else if(last_fd > depth[index]) last_fd=depth[index];

    next_fd = fdepth[0] - Zsum;
    if(next_fd < 0) next_fd=0.;
    else if(next_fd > depth[index]) next_fd=depth[index];

    if(next_td > last_td && (last_fd - last_td)> 0.) {
      /** Thawing layer advances **/
      layer[index].moist_thaw += (next_td - last_td)
          / (last_fd - last_td) * layer[index].moist_froz;
      layer[index].moist_froz -= (next_td - last_td)
          / (last_fd - last_td) * layer[index].moist_froz;
    }
    else if(next_td > last_td && (last_fd - last_td) == 0.) {
      /** Thawing and Freezing Layer Reset by New Storm **/
      layer[index].moist_thaw += (next_td - last_td)
          / depth[index] * layer[index].moist;
      layer[index].moist_froz -= (next_td - last_td)
          / depth[index] * layer[index].moist;
    }
    else if(next_td < last_td) {
      /** Thawing layer retreats **/
      layer[index].moist_froz += (last_td - next_td)
          / (last_td) * layer[index].moist_thaw;
      layer[index].moist_thaw -= (last_td - next_td)
          / (last_td) * layer[index].moist_thaw;
    }
 
    if(next_fd > last_fd) {
      /** Freezing layer advances **/
      layer[index].moist_froz += (next_fd - last_fd)
          / (depth[index] - last_fd) * layer[index].moist;
      layer[index].moist -= (next_fd - last_fd)
          / (depth[index] - last_fd) * layer[index].moist;
    }
    else if(next_fd < last_fd) {
      /** Thawing layer retreats **/
      layer[index].moist += (last_fd - next_fd)
          / (last_fd - next_td) * layer[index].moist_froz;
      layer[index].moist_froz -= (last_fd - next_fd)
          / (last_fd - next_td) * layer[index].moist_froz;
    }
   
    temp_moist = 0.0;
    if(next_td > 0.0) {
      layer[index].moist_thaw /= (next_td) / depth[index];
      if(layer[index].moist_thaw > max_moist[index]) {
        temp_moist = max_moist[index]-layer[index].moist_thaw;
        layer[index].moist_thaw = max_moist[index];
      }
    }
    else layer[index].moist_thaw = 0.;
    if(next_fd > 0.0 && next_td < depth[index]) {
      layer[index].moist_froz /= (next_fd
          - next_td) / depth[index];
      layer[index].moist_froz += temp_moist;
      if(layer[index].moist_froz > max_moist[index]) {
        temp_moist = max_moist[index]-layer[index].moist_froz;
        layer[index].moist_froz = max_moist[index];
      }
      else temp_moist = 0.0;
    }
    else layer[index].moist_froz = 0.;
    if(next_fd < depth[index]) {
      layer[index].moist /= (depth[index] - next_fd) / depth[index];
      layer[index].moist += temp_moist;
    }
    else layer[index].moist = 0.;

    if(layer[index].moist_thaw < 0.) layer[index].moist_thaw=0.;
    if(layer[index].moist_froz < 0.) layer[index].moist_froz=0.;
    if(layer[index].moist < 0.) layer[index].moist=0.;

    Zsum += depth[index];
   
  }
}

void find_sublayer_temperatures(layer_data_struct *layer,
                                double            *T,
                                double            *dz,
                                double            *depth,
                                double             fdepth,
                                double             tdepth,
                                int                Nlayer,
                                int                Nnodes) {
/**********************************************************************
  This subroutine computes the soil temperature for each layer and 
  it's subdivisions.
**********************************************************************/

  int i, lindex, zindex;
  double Zsum;   /** Total depth of soil thermal layers from surface **/
  double Lsum;   /** Total depth of soil moisture layers from surface **/
  double Ztmp;   /** Depth of current soil thermal layer **/
  double Ltmp;   /** Thickness of soil moisture layer below last 
		     thermal layer **/
  double Ftmp;   /** Depth from current soil thermal layer to frost depth **/ 
  double Ttmp;   /** Depth from current soil thermal layer to thaw depth **/
  double Tstr;   /** Stored temperature for current soil thermal layer **/
  double Tavg;   /** Averaged temperature for current soil moisture 
		     sublayer **/
  double tmpsum; /** Thickness of stored temperatures (Tstr) for current 
		     layer **/
  double T1;     /** Temperature of last thermal node **/
  double T2;     /** Temperature of current soil thermal node **/

  Zsum=0.;
  Lsum = depth[0];
  lindex=0;

  for(i=0;i<Nlayer;i++)
    layer[i].T = layer[i].T_froz = layer[i].T_thaw = -999.;
  Tstr = 0.;
  Tavg = 0.;
  tmpsum = 0.;

  for(zindex=0;zindex<Nnodes-1;zindex++) {
    if(lindex<Nlayer) {
      Ltmp = Lsum - Zsum;
      Ftmp = fdepth - Zsum;
      Ttmp = tdepth - Zsum;
      Ztmp = (dz[zindex]+dz[zindex+1])/2.;
      
      T1 = T[zindex];
      T2 = T[zindex+1];
      
      if(Ttmp>Ztmp && Ltmp>Ztmp) {
	Tstr = 0.5*(T1+T2)*Ztmp;
	tmpsum += Ztmp;
      }
      else if(Ttmp<=0 && Ftmp>Ztmp && Ltmp>Ztmp) {
	Tstr = 0.5*(T1+T2)*Ztmp;
	tmpsum += Ztmp;
      }
      else if(Ttmp<=0 && Ftmp<=0 && Ltmp>Ztmp) {
	Tstr = 0.5*(T1+T2)*Ztmp;
	tmpsum += Ztmp;
      }
      else if(Ttmp>0 && Ttmp<=Ztmp && Ftmp>Ztmp && Ltmp>Ztmp) {
	Tavg += 0.5*(T1+0.)*Ttmp;
	tmpsum += Ttmp;
	if(tmpsum>0.) Tavg /= tmpsum;
	else Tavg= -999.;
	layer[lindex].T_thaw = Tavg;
	Tavg = 0.;
	tmpsum = Ztmp - Ttmp;
	Tstr = 0.5*(0.+T2)*tmpsum;
      }
      else if(Ttmp<=0 && Ftmp>0 && Ftmp<=Ztmp && Ltmp>Ztmp) {
	Tavg += 0.5*(T1+0.)*Ftmp;
	tmpsum += Ftmp;
	if(tmpsum>0.) Tavg /= tmpsum;
	else Tavg= -999.;
	layer[lindex].T_froz = Tavg;
	Tavg = 0.;
	tmpsum = Ztmp - Ftmp;
	Tstr = 0.5*(0.+T2)*tmpsum;
      }
      else if(Ttmp<=0 && Ftmp<=0 && Ltmp>0 && Ltmp<=Ztmp) {
	Tavg += 0.5*(T1+linear_interp(Ltmp,0.,Ztmp,T1,T2))*Ltmp;
	tmpsum += Ltmp;
	if(tmpsum>0.) Tavg /= tmpsum;
	else Tavg= -999.;
	layer[lindex].T = Tavg;
	Tavg = 0.;
	tmpsum = Ztmp - Ltmp;
	Tstr = 0.5*(linear_interp(Ltmp,0.,Ztmp,T1,T2)+T2)*tmpsum;
      }
      else if(Ttmp<=0 && Ftmp>Ztmp && Ltmp>0 && Ltmp<=Ztmp) {
	Tavg += 0.5*(T1+linear_interp(Ltmp,0.,Ztmp,T1,T2))*Ltmp;
	tmpsum += Ltmp;
	if(tmpsum>0.) Tavg /= tmpsum;
	else Tavg= -999.;
	layer[lindex].T_froz = Tavg;
	Tavg = 0.;
	tmpsum = Ztmp - Ltmp;
	Tstr = 0.5*(linear_interp(Ltmp,0.,Ztmp,T1,T2)+T2)*tmpsum;
      }
      else if(Ttmp>Ztmp && Ftmp>Ztmp && Ltmp>0 && Ltmp<=Ztmp) {
	Tavg += 0.5*(T1+linear_interp(Ltmp,0.,Ztmp,T1,T2))*Ltmp;
	tmpsum += Ltmp;
	if(tmpsum>0.) Tavg /= tmpsum;
	else Tavg= -999.;
	layer[lindex].T_thaw = Tavg;
	Tavg = 0.;
	tmpsum = Ztmp - Ltmp;
	Tstr = 0.5*(linear_interp(Ltmp,0.,Ztmp,T1,T2)+T2)*tmpsum;
      }
      else if(Ttmp>0 && Ftmp>Ztmp && Ltmp>Ttmp && Ltmp<=Ztmp) {
	Tavg += 0.5*(T1+linear_interp(Ttmp,0.,Ztmp,T1,T2))*Ttmp;
	tmpsum += Ttmp;
	if(tmpsum>0.) Tavg /= tmpsum;
	else Tavg= -999.;
	layer[lindex].T_thaw = Tavg;
	layer[lindex].T_froz = 0.5*(linear_interp(Ttmp,0.,Ztmp,T1,T2)
				    +linear_interp(Ltmp,0.,Ztmp,T1,T2));
	Tavg = 0.;
	tmpsum = Ztmp - Ltmp;
	Tstr = 0.5*(linear_interp(Ltmp,0.,Ztmp,T1,T2)+T2)*tmpsum;
      }
      else if(Ttmp>Ltmp && Ttmp<=Ztmp && Ftmp>Ztmp && Ltmp>0) {
	Tavg += 0.5*(T1+linear_interp(Ltmp,0.,Ztmp,T1,T2))*Ltmp;
	tmpsum += Ltmp;
	if(tmpsum>0.) Tavg /= tmpsum;
	else Tavg= -999.;
	layer[lindex].T_thaw = Tavg;
	layer[lindex+1].T_thaw = 0.5*(linear_interp(Ltmp,0.,Ztmp,T1,T2)
				      +linear_interp(Ttmp,0.,Ztmp,T1,T2));
	Tavg = 0.;
	tmpsum = Ztmp - Ttmp;
	Tstr = 0.5*(linear_interp(Ttmp,0.,Ztmp,T1,T2)+T2)*tmpsum;
      }
      else if(Ttmp<=0 && Ftmp>0 && Ltmp>Ftmp && Ltmp<=Ztmp) {
	Tavg += 0.5*(T1+linear_interp(Ftmp,0.,Ztmp,T1,T2))*Ftmp;
	tmpsum += Ftmp;
	if(tmpsum>0.) Tavg /= tmpsum;
	else Tavg= -999.;
	layer[lindex].T_froz = Tavg;
	layer[lindex].T = 0.5*(linear_interp(Ftmp,0.,Ztmp,T1,T2)
			       +linear_interp(Ltmp,0.,Ztmp,T1,T2));
	Tavg = 0.;
	tmpsum = Ztmp - Ltmp;
	Tstr = 0.5*(linear_interp(Ltmp,0.,Ztmp,T1,T2)+T2)*tmpsum;
      }
      else if(Ttmp<=0 && Ftmp>Ltmp && Ftmp<=Ztmp && Ltmp>0) {
	Tavg += 0.5*(T1+linear_interp(Ltmp,0.,Ztmp,T1,T2))*Ltmp;
	tmpsum += Ltmp;
	if(tmpsum>0.) Tavg /= tmpsum;
	else Tavg= -999.;
	layer[lindex].T_froz = Tavg;
	layer[lindex+1].T_froz = 0.5*(linear_interp(Ltmp,0.,Ztmp,T1,T2)
				      +linear_interp(Ftmp,0.,Ztmp,T1,T2));
	Tavg = 0.;
	tmpsum = Ztmp - Ftmp;
	Tstr = 0.5*(linear_interp(Ftmp,0.,Ztmp,T1,T2)+T2)*tmpsum;
      }
      else if(Ttmp>0 && Ftmp>Ttmp && Ftmp<=Ztmp && Ltmp>Ztmp) {
	Tavg += 0.5*(T1+linear_interp(Ttmp,0.,Ztmp,T1,T2))*Ttmp;
	tmpsum += Ttmp;
	if(tmpsum>0.) Tavg /= tmpsum;
	else Tavg= -999.;
	layer[lindex].T_thaw = Tavg;
	layer[lindex].T_froz = 0.5*(linear_interp(Ttmp,0.,Ztmp,T1,T2)
				    +linear_interp(Ftmp,0.,Ztmp,T1,T2));
	Tavg = 0.;
	tmpsum = Ztmp - Ftmp;
	Tstr = 0.5*(linear_interp(Ftmp,0.,Ztmp,T1,T2)+T2)*tmpsum;
      }
      else if(Ttmp>0 && Ftmp>Ttmp && Ltmp>Ftmp && Ltmp<=Ztmp) {
	Tavg += 0.5*(T1+linear_interp(Ttmp,0.,Ztmp,T1,T2))*Ttmp;
	tmpsum += Ttmp;
	if(tmpsum>0.) Tavg /= tmpsum;
	else Tavg= -999.;
	layer[lindex].T_thaw = Tavg;
	layer[lindex].T_froz = 0.5*(linear_interp(Ttmp,0.,Ztmp,T1,T2)
				    +linear_interp(Ftmp,0.,Ztmp,T1,T2));
	layer[lindex].T = 0.5*(linear_interp(Ftmp,0.,Ztmp,T1,T2)
			       +linear_interp(Ltmp,0.,Ztmp,T1,T2));
	Tavg = 0.;
	tmpsum = Ztmp - Ltmp;
	Tstr = 0.5*(linear_interp(Ltmp,0.,Ztmp,T1,T2)+T2)*tmpsum;
      }
      else if(Ttmp>0 && Ftmp>Ltmp && Ftmp<=Ztmp && Ltmp>Ttmp) {
	Tavg += 0.5*(T1+linear_interp(Ttmp,0.,Ztmp,T1,T2))*Ttmp;
	tmpsum += Ttmp;
	if(tmpsum>0.) Tavg /= tmpsum;
	else Tavg= -999.;
	layer[lindex].T_thaw = Tavg;
	layer[lindex].T_froz = 0.5*(linear_interp(Ttmp,0.,Ztmp,T1,T2)
				    +linear_interp(Ltmp,0.,Ztmp,T1,T2));
	layer[lindex+1].T_froz = 0.5*(linear_interp(Ltmp,0.,Ztmp,T1,T2)
				      +linear_interp(Ftmp,0.,Ztmp,T1,T2));
	Tavg = 0.;
	tmpsum = Ztmp - Ftmp;
	Tstr = 0.5*(linear_interp(Ftmp,0.,Ztmp,T1,T2)+T2)*tmpsum;
      }
      else if(Ttmp>Ltmp && Ftmp>Ttmp && Ftmp<=Ztmp && Ltmp>0) {
	Tavg += 0.5*(T1+linear_interp(Ltmp,0.,Ztmp,T1,T2))*Ltmp;
	tmpsum += Ltmp;
	if(tmpsum>0.) Tavg /= tmpsum;
	else Tavg= -999.;
	layer[lindex].T_thaw = Tavg;
	layer[lindex].T_thaw = 0.5*(linear_interp(Ltmp,0.,Ztmp,T1,T2)
				    +linear_interp(Ttmp,0.,Ztmp,T1,T2));
	layer[lindex+1].T_froz = 0.5*(linear_interp(Ttmp,0.,Ztmp,T1,T2)
				      +linear_interp(Ftmp,0.,Ztmp,T1,T2));
	Tavg = 0.;
	tmpsum = Ztmp - Ftmp;
	Tstr = 0.5*(linear_interp(Ftmp,0.,Ztmp,T1,T2)+T2)*tmpsum;
      }

      Tavg += Tstr;

      if(Ltmp>0 && Ltmp<=Ztmp) {
	lindex++;
	if(lindex<Nlayer) Lsum+=depth[lindex];
      }
      Zsum += Ztmp;
    }
  }
  if(Lsum>Zsum && lindex<Nlayer)
    layer[lindex].T = ((Lsum-Zsum)*T[Nnodes-1] + Tavg) / depth[lindex];
  if(lindex<Nlayer-1) {
    for(i=lindex;i<Nlayer;i++)
      layer[lindex].T = T[Nnodes-1];
  }

  for(i=0;i<Nlayer;i++) {
    if(layer[i].tdepth == 0.) layer[i].T_thaw = -999.;
    if(layer[i].fdepth == 0. || layer[i].tdepth == depth[i])
      layer[i].T_froz = -999.;
    else if(layer[i].T_froz > 0.) { 
      layer[i].T_froz = 0.;
      fprintf(stderr,"WARNING: Frozen layer temperature for layer %i was set to 0. thickness = %.4lf\n",i,layer[i].fdepth-layer[i].tdepth);
    }
    if(layer[i].fdepth == depth[i]) layer[i].T = -999.;
  }

}

layer_data_struct find_average_layer(layer_data_struct wet,
				     layer_data_struct dry,
				     double            depth,
				     double            mu) {
/*************************************************************
  This subroutine computes the average soil layer moistures
  between the wet and dry fraction for use in computing 
  energy balance parameters.  Other layer variables are copied 
  from the wet fraction structure since they are they same for 
  wet and dry fractions.
**************************************************************/

  layer_data_struct layer;

  layer = wet;

  layer.moist_thaw = ((wet.moist_thaw * mu) 
		      + (dry.moist_thaw * (1. - mu)));
  layer.moist_froz = ((wet.moist_froz * mu)
		      + (dry.moist_froz * (1. - mu)));
  layer.ice = ((wet.ice * mu) + (dry.ice * (1. - mu)));
  layer.moist = ((wet.moist * mu)
	      + (dry.moist * (1. - mu)));

  return(layer);

}

double find_total_layer_moisture(layer_data_struct layer,
				 double            depth) {
/*****************************************************************
  This subroutine returns the total soil moisture in a layer 
  (including ice).
*****************************************************************/

  double moist;

  moist  = layer.moist_thaw * layer.tdepth / depth;
  moist += layer.moist_froz
         * (layer.fdepth - layer.tdepth) / depth;
  moist += layer.ice
         * (layer.fdepth - layer.tdepth) / depth;
  moist += layer.moist * (depth - layer.fdepth) 
         / depth;

  return (moist);
}

double find_total_layer_ice(layer_data_struct layer,
			    double            depth) {
/*****************************************************************
  This subroutine computes the total moisture in a layer
*****************************************************************/

  double ice;

  ice = layer.ice * (layer.fdepth - layer.tdepth) / depth;

  return (ice);

}
