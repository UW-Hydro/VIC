#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
#include <string.h>

void initialize_energy_bal (energy_bal_struct *energy, 
                            cell_data_struct *cell,
                            soil_con_struct soil_con,
                            double surf_temp,
		            int veg_num,
                            int Ulayer,
                            int Llayer,
                            FILE *fsoil)
/**********************************************************************
	initialize_energy_bal	Keith Cherkauer		July 30, 1996

  This routine initializes the energy balance data structure for use
  with the frozen soils model.

  Soil temperatures are initialized using a linear interpolation
  between the surface temperature (assumed = air temperature) and
  the average annual air temperature (assumed that bottom of 2nd soil
  layer is at const T = average annual air temp).

  UNITS: (m, s, kg, C, moisture in mm) unless otherwise specified

  variable modified:
	energy[veg].dz
	energy[veg].T
	energy[veg].fdepth
	cell[veg].layer[index].T
	cell[veg].layer[index].T_thaw
	cell[veg].layer[index].T_froz
	cell[veg].layer[index].moist
	cell[veg].layer[index].moist_thaw
	cell[veg].layer[index].moist_froz
	cell[veg].layer[index].ice
	cell[veg].layer[index].tdepth
	cell[veg].layer[index].fdepth
	cell[veg].layer[index].kappa
	cell[veg].layer[index].Cs

  Modifications:
  10-15-96  modified to reflect fixes in the frozen soils code	KAC
  07-10-97  modified to read water content in (m/m) instead of (mm) KAC
  09-18-97  modified to initialize with no files, and to initialize
		FULL_ENERGY case without FROZEN_SOIL		KAC

**********************************************************************/
{
  extern option_struct options;
  extern debug_struct debug;

  char tmpstr[MAXSTRING];
  char errorstr[MAXSTRING];
  int i, j, veg, index, Tlayer;
  int tmpint;
  double sum, Lsum, Zsum, dp, Ltotal, tmpdepth;
  double *kappa, *Cs, *M;
  double moist[MAXlayer], ice[MAXlayer];
  double unfrozen, frozen;
  double *thermdepths;
  double **layer_ice;

  Tlayer = Ulayer+Llayer+2;
  dp = soil_con.dp;
  Ltotal = 0;
  for(index=0;index<options.Nlayer;index++) Ltotal += soil_con.depth[index];

  /** SOIL INITIALIZATION **/
  if(options.INIT_SOIL && options.FROZEN_SOIL) {

    thermdepths = (double *)calloc(Tlayer,sizeof(double));

    for ( veg = 0 ; veg <= veg_num ; veg++) {
      rewind(fsoil);
      fgets(tmpstr,MAXSTRING,fsoil);

      fscanf(fsoil,"%*s %i\n",&tmpint);
      if(tmpint!=Ulayer) nrerror("Ulayer in soil initialization file, does not match that of the current model version");
      fscanf(fsoil,"%*s %i\n",&tmpint);
      if(tmpint!=Llayer) nrerror("Llayer in soil initialization file, does not match that of the current model version");

      fscanf(fsoil,"%*s %lf %lf\n",&energy[veg].fdepth[0],
          &energy[veg].fdepth[1]);
      fscanf(fsoil,"%*s");
      for(i=0;i<options.Nlayer;i++) fscanf(fsoil,"%lf",&moist[i]);
      fscanf(fsoil,"%*s");
      for(i=0;i<options.Nlayer;i++) fscanf(fsoil,"%lf",&ice[i]);

      fgets(tmpstr,MAXSTRING,fsoil);
      fgets(tmpstr,MAXSTRING,fsoil);
      for(j=0;j<Tlayer;j++) {
        fscanf(fsoil,"%lf",&thermdepths[j]);
        fscanf(fsoil,"%lf",&energy[veg].T[j]);
      }
      if(thermdepths[Ulayer+1] != dp) {
        sprintf(errorstr,"Thermal solution depth %i (Ulayer+1) must equal thermal damping depth %lf, but is equal to %lf",Ulayer+1,dp,thermdepths[Ulayer+1]);
        nrerror(errorstr);
      }
      if(thermdepths[Tlayer-1] != Ltotal) {
        sprintf(errorstr,"Last defined thermal depth (%lf) must equal bottom of soil layers (%lf)",thermdepths[Tlayer-1],Ltotal);
        nrerror(errorstr);
      }
      for(j=Tlayer-1;j>0;j--) {
        thermdepths[j] -= thermdepths[j-1];
        thermdepths[j] = (double)((int)(thermdepths[j] * 10000. + 0.5))
                       / 10000.;
      }
      for(j=1;j<Tlayer-1;j++) {
        if((int)(thermdepths[j]*1000.+0.5) == (int)(thermdepths[j+1]*1000.+0.5))
          energy[veg].dz[j] = (thermdepths[j] + thermdepths[j+1]) / 2.;
        else {
          energy[veg].dz[j] = thermdepths[j];
          j++;
          energy[veg].dz[j] = thermdepths[j+1];
          if((int)((energy[veg].dz[j-1]+energy[veg].dz[j])/2.*1000.+0.5)
              != (int)(thermdepths[j]*1000.+0.5)) {
            sprintf(errorstr,"Check spacing between thermal layers %i and %i\n",
                    j-1,j);
            nrerror(errorstr);
          }
        }
      }
      energy[veg].dz[0] = thermdepths[1];
      energy[veg].dz[Tlayer-1] = thermdepths[Tlayer-1];

      sum = energy[veg].dz[0]/2. + energy[veg].dz[Tlayer-1]/2.;
      for(index=1;index<Tlayer-1;index++) sum += energy[veg].dz[index];
      if((int)(Ltotal*1000.+0.5) != (int)(sum*1000.+0.5)) {
        sprintf(errorstr,"Sum of initial dz values (%lf) is not equal to the total depth (%lf).",sum*1000.+0.5,Ltotal*1000.+0.5);
        nrerror(errorstr);
      }

      Lsum=0.;
      for(index=0;index<options.Nlayer;index++) {
        Lsum += soil_con.depth[index];
        if(energy[veg].fdepth[1] > (Lsum - soil_con.depth[index])) {
          if(energy[veg].fdepth[1] < Lsum) {
            cell[veg].layer[index].tdepth = energy[veg].fdepth[1] 
                - (Lsum - soil_con.depth[index]);
            cell[veg].layer[index].moist_thaw = moist[index] 
                * soil_con.depth[index] * 1000.;
            if(energy[veg].fdepth[0] < Lsum) {
              cell[veg].layer[index].fdepth = energy[veg].fdepth[0] 
                  - (Lsum - soil_con.depth[index]);
              cell[veg].layer[index].moist_froz = soil_con.depth[index]
                  * moist[index] * 1000.;
              cell[veg].layer[index].ice = soil_con.depth[index]
                  * ice[index] * 1000.;
              cell[veg].layer[index].moist = soil_con.depth[index] 
                  * moist[index] * 1000.;
            }
            else {
              cell[veg].layer[index].fdepth = soil_con.depth[index];
              cell[veg].layer[index].moist_froz = soil_con.depth[index]
                  * moist[index] * 1000.;
              cell[veg].layer[index].ice = soil_con.depth[index]
                  * ice[index] * 1000.;
              cell[veg].layer[index].moist = 0.;
            }
          }
          else {
            cell[veg].layer[index].tdepth = soil_con.depth[index];
            cell[veg].layer[index].fdepth = soil_con.depth[index];
            cell[veg].layer[index].moist_thaw = soil_con.depth[index]
                * moist[index] * 1000.;
            cell[veg].layer[index].moist_froz = 0.;
            cell[veg].layer[index].ice = 0.;
            cell[veg].layer[index].moist = 0.;
          }
        }
        else if(energy[veg].fdepth[0] > (Lsum - soil_con.depth[index])) {
          if(energy[veg].fdepth[0] < Lsum) {
            cell[veg].layer[index].tdepth = 0.;
            cell[veg].layer[index].fdepth = energy[veg].fdepth[0]
                - (Lsum - soil_con.depth[index]);
            cell[veg].layer[index].moist_thaw = 0.;
            cell[veg].layer[index].moist_froz = soil_con.depth[index]
                * moist[index] * 1000.;
            cell[veg].layer[index].ice = soil_con.depth[index]
                * ice[index] * 1000.;
            cell[veg].layer[index].moist = soil_con.depth[index] 
                * moist[index] * 1000.;
          }
          else {
            cell[veg].layer[index].tdepth = 0.;
            cell[veg].layer[index].fdepth = soil_con.depth[index];
            cell[veg].layer[index].moist_thaw = 0.;
            cell[veg].layer[index].moist_froz = soil_con.depth[index]
                * moist[index] * 1000.;
            cell[veg].layer[index].ice = soil_con.depth[index]
                * ice[index] * 1000.;
            cell[veg].layer[index].moist = 0.;
          }
        }
        else {
          cell[veg].layer[index].tdepth = 0.;
          cell[veg].layer[index].fdepth = 0.;
          cell[veg].layer[index].moist_thaw = 0.;
          cell[veg].layer[index].moist_froz = 0.;
          cell[veg].layer[index].ice = 0.;
          cell[veg].layer[index].moist = soil_con.depth[index] 
              * moist[index] * 1000.;
        }
      }
      layer_ice = (double **)calloc(options.Nlayer,sizeof(double *));
      for(i=0;i<options.Nlayer;i++) {
        layer_ice[i] = (double *)calloc(3,sizeof(double));
        layer_ice[i][0] = 0.;
        layer_ice[i][1] = cell[veg].layer[i].ice / (soil_con.depth[i] * 1000.);
        layer_ice[i][2] = 0.;
      }
      distribute_soil_property(energy[veg].dz,energy[veg].fdepth[0],
          energy[veg].fdepth[1],
          layer_ice,options.Nlayer,Tlayer,soil_con.depth,energy[veg].ice);
      for(i=0;i<options.Nlayer;i++) free((char *)layer_ice[i]);
      free((char *)layer_ice);
  
    }

    free((char *)thermdepths);

  }
  else if(options.INIT_SOIL && options.FULL_ENERGY) {

    thermdepths = (double *)calloc(Tlayer,sizeof(double));

    for ( veg = 0 ; veg <= veg_num ; veg++) {
      rewind(fsoil);
      fgets(tmpstr,MAXSTRING,fsoil);

      fscanf(fsoil,"%*s %i\n",&tmpint);
      if(tmpint!=Ulayer) nrerror("Ulayer in soil initialization file, does not match that of the current model version");
      fscanf(fsoil,"%*s %i\n",&tmpint);
      if(tmpint!=Llayer) nrerror("Llayer in soil initialization file, does not match that of the current model version");

      energy[veg].fdepth[0]=energy[veg].fdepth[1]=0.;
      for(i=0;i<options.Nlayer;i++) {
        cell[veg].layer[i].fdepth = 0.;
        cell[veg].layer[i].tdepth = 0.;
      }
      fscanf(fsoil,"%*s");
      for(i=0;i<options.Nlayer;i++) {
        fscanf(fsoil,"%lf",&cell[veg].layer[i].moist);
        cell[veg].layer[i].moist *= soil_con.depth[i]*1000.;
        cell[veg].layer[i].moist_froz = 0.;
        cell[veg].layer[i].moist_thaw = 0.;
      }
      fgets(tmpstr,MAXSTRING,fsoil);
      fgets(tmpstr,MAXSTRING,fsoil);
      for(j=0;j<Tlayer;j++) {
        fscanf(fsoil,"%lf",&thermdepths[j]);
        fscanf(fsoil,"%lf",&energy[veg].T[j]);
      }
      if(thermdepths[Ulayer+1] != dp) {
        sprintf(errorstr,"Thermal solution depth %i (Ulayer+1) must equal thermal damping depth %lf, but is equal to %lf",Ulayer+1,dp,thermdepths[Ulayer+1]);
        nrerror(errorstr);
      }
      if(thermdepths[Tlayer-1] != Ltotal) {
        sprintf(errorstr,"Last defined thermal depth (%lf) must equal bottom of soil layers (%lf)",thermdepths[Tlayer-1],Ltotal);
        nrerror(errorstr);
      }
      for(j=Tlayer-1;j>0;j--) {
        thermdepths[j] -= thermdepths[j-1];
        thermdepths[j] = (double)((int)(thermdepths[j] * 10000. + 0.5))
                       / 10000.;
      }

      energy[veg].dz[0] = soil_con.depth[0];
      energy[veg].dz[1] = soil_con.depth[0];
      energy[veg].dz[2] = 2. * (dp - 1.5 * soil_con.depth[0]);
      energy[veg].dz[2] = 2. * (Ltotal - 1.5 * soil_con.depth[0]
                        - energy[veg].dz[2]);
    }

    free((char *)thermdepths);

  }

  else if(options.FULL_ENERGY && !options.FROZEN_SOIL) {

    for ( veg = 0 ; veg <= veg_num ; veg++) {
      energy[veg].fdepth[0] = 0.;
      energy[veg].fdepth[1] = 0.;
      for(index=0;index<options.Nlayer;index++) {
        cell[veg].layer[index].fdepth = 0.;
        cell[veg].layer[index].tdepth = 0.;
        cell[veg].layer[index].T = 5.;	/* just needs to be positive */
        cell[veg].layer[index].T_thaw = -999.;
        cell[veg].layer[index].T_froz = -999.;
        cell[veg].layer[index].moist_thaw = 0.;
        cell[veg].layer[index].moist_froz = 0.;
        cell[veg].layer[index].ice = 0.;
      }
      energy[veg].dz[0] = soil_con.depth[0];
      energy[veg].dz[1] = soil_con.depth[0];
      energy[veg].dz[2] = 2. * (dp - 1.5 * soil_con.depth[0]);
      energy[veg].dz[3] = 2. * (Ltotal - energy[veg].dz[2] 
                        - 1.5 * soil_con.depth[0]);
      energy[veg].T[0] = surf_temp;
      energy[veg].T[1] = surf_temp;
      energy[veg].T[2] = soil_con.avg_temp;
      energy[veg].T[3] = soil_con.avg_temp;
    }
  }

  /*****************************************************************
    Initialize Energy Balance Variables if Frozen Soils Activated,
    and no Initial Condition File Given 
  *****************************************************************/
  else if(options.FROZEN_SOIL) {
    for ( veg = 0 ; veg <= veg_num ; veg++) {

      energy[veg].T[0] = surf_temp;
      energy[veg].dz[0] = energy[veg].dz[1] = soil_con.depth[0];
      energy[veg].T[Tlayer-1] = soil_con.avg_temp;
      energy[veg].T[1] = linear_interp(soil_con.depth[0],0.,
                         dp,surf_temp,soil_con.avg_temp);

      Zsum = 0.;
      if(dp<Ltotal) tmpdepth = dp;
      else tmpdepth = Ltotal;
      dp -= soil_con.depth[0] * 1.5;
      tmpdepth -= soil_con.depth[0] * 1.5;
      for(index=2;index<=Ulayer+1;index++) {
        energy[veg].dz[index] = tmpdepth/(((double)Ulayer-0.5));
        Zsum += (energy[veg].dz[index]+energy[veg].dz[index-1])/2.;
        energy[veg].T[index] = soil_con.avg_temp;
      }
      dp += soil_con.depth[0] * 1.5;
      if(dp < Ltotal) {
	for(index=Ulayer+2;index<Tlayer;index++) {
	  energy[veg].dz[index] = (Ltotal - dp - energy[veg].dz[Ulayer+1]/2.)
	    / ((double)Llayer-0.5);
	  energy[veg].T[index] = soil_con.avg_temp;
	}
      }
      else {
	for(index=Ulayer+2;index<Tlayer;index++) {
	  energy[veg].dz[index] = 0;
	  energy[veg].T[index] = soil_con.avg_temp;
	}
      }

      find_0_degree_fronts(&energy[veg],cell[veg].layer,Ltotal,
                           soil_con.depth,energy[veg].T,options.Nlayer,
                           Tlayer);

      find_sublayer_temperatures(cell[veg].layer,energy[veg].T,energy[veg].dz,
                                 soil_con.depth,energy[veg].fdepth[0],
                                 energy[veg].fdepth[1],options.Nlayer,Tlayer);
 
      for(index=0;index<options.Nlayer;index++) {
        if(cell[veg].layer[index].tdepth > 0.)
          cell[veg].layer[index].moist_thaw = cell[veg].layer[index].moist;
        if(cell[veg].layer[index].tdepth == soil_con.depth[index])
          cell[veg].layer[index].moist = 0.;
        else {
          if(cell[veg].layer[index].fdepth > 0.)
            cell[veg].layer[index].moist_froz = cell[veg].layer[index].moist;
          if(cell[veg].layer[index].fdepth == soil_con.depth[index])
            cell[veg].layer[index].moist = 0.;
        }
        if(cell[veg].layer[index].T_froz<0.
            && cell[veg].layer[index].T_froz != -999.) {
          unfrozen = maximum_unfrozen_water(cell[veg].layer[index].T_froz,
              soil_con.max_moist[index], soil_con.bubble,
              soil_con.expt[index]);
          if(unfrozen>soil_con.max_moist[index] || unfrozen<0.)
            unfrozen = soil_con.max_moist[index];
          cell[veg].layer[index].unfrozen = unfrozen;
 
          frozen = cell[veg].layer[index].moist_froz - unfrozen;
          if(frozen < 0.0) {
            frozen = 0.0;
            unfrozen = cell[veg].layer[index].moist_froz;
          }
          cell[veg].layer[index].ice = frozen;
          cell[veg].layer[index].moist_froz = unfrozen;
        }
        else if(cell[veg].layer[index].T_froz == 0.) {
          cell[veg].layer[index].unfrozen = soil_con.max_moist[index];
          cell[veg].layer[index].ice = 0.;
        }
        else if(cell[veg].layer[index].T_froz != -999.)
          vicerror("ERROR: Frozen Layer Temperature > 0C");

      }

      layer_ice = (double **)calloc(options.Nlayer,sizeof(double *));
      for(i=0;i<options.Nlayer;i++) {
        layer_ice[i] = (double *)calloc(3,sizeof(double));
        layer_ice[i][0] = 0.;
        layer_ice[i][1] = cell[veg].layer[i].ice / (soil_con.depth[i] * 1000.);
        layer_ice[i][2] = 0.;
      }
      distribute_soil_property(energy[veg].dz,energy[veg].fdepth[0],
          energy[veg].fdepth[1],
          layer_ice,options.Nlayer,Tlayer,soil_con.depth,energy[veg].ice);
      for(i=0;i<options.Nlayer;i++) free((char *)layer_ice[i]);
      free((char *)layer_ice);

    }

  }
  else {
    for ( veg = 0 ; veg <= veg_num ; veg++) {
      for(index=0;index<options.Nlayer;index++) {
        energy[veg].dz[index] = 1.;
      }
    }
  }

  /** set initial soil properties for energy balance **/

  if(options.FULL_ENERGY || options.FROZEN_SOIL) {
    for ( veg = 0 ; veg <= veg_num ; veg++) {
      index = 0;
      Zsum = energy[veg].dz[index]/2.;
      while(index<=Ulayer && (int)((Zsum+energy[veg].dz[index+1]/2.)*1000.+0.5) 
                             < (int)((soil_con.depth[0])*1000.+0.5)) {
        index++;
        Zsum += energy[veg].dz[index];
      }
      if(index>Ulayer || (int)((Zsum+energy[veg].dz[index+1]/2.)*1000.+0.5) 
                             > (int)((soil_con.depth[0])*1000.+0.5)) {
        nrerror("One soil temperature depth must equal the depth of the top soil layer");
      }
      energy[veg].T1_index = index+1;

      find_sublayer_temperatures(cell[veg].layer,energy[veg].T,
          energy[veg].dz,soil_con.depth,energy[veg].fdepth[0],
          energy[veg].fdepth[1],options.Nlayer,Tlayer);

      kappa = NULL;
      Cs = NULL;
      M = NULL;
      soil_thermal_calc(soil_con, cell[veg].layer, energy[veg], kappa,
          Cs, M, M, M, options.Nlayer, Tlayer);

      /** Save Thermal Conductivities for Energy Balance **/
      energy[veg].kappa[0] = cell[veg].layer[0].kappa;
      energy[veg].Cs[0] = cell[veg].layer[0].Cs;
      energy[veg].kappa[1] = cell[veg].layer[1].kappa;
      energy[veg].Cs[1] = cell[veg].layer[1].Cs;

      if(energy[veg].fdepth[0]>0.) energy[veg].frozen=TRUE;
    }

    if(debug.DEBUG || debug.PRT_TEMP) {
      fprintf(debug.fg_temp,"%i\n",Tlayer);
      fprintf(debug.fg_temp,"Date - hour(REC)        \tAir T\tFdpth\tTdpth");
      for(i=0;i<options.Nlayer;i++) fprintf(debug.fg_temp,"\t%i Th T\t%i Fr T\t%i T",i,i,i);
      sum=0.0;
      fprintf(debug.fg_temp,"\t2.5cm");
      for(i=0;i<5;i++) {
        fprintf(debug.fg_temp,"\t%.0lf",sum*100.0);
        sum+=(energy[veg_num].dz[i]+energy[veg_num].dz[i+1])/2.0;
      }
      fprintf(debug.fg_temp,"\n");
    }
  }

}
