#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vicNl.h>

void write_debug(atmos_data_struct     atmos,
                 soil_con_struct       soil_con,
                 cell_data_struct     *cell,
                 energy_bal_struct    *energy,
                 snow_data_struct     *snow,
                 veg_var_struct       *veg_var,
                 dmy_struct            dmy,
                 global_param_struct   gp,
                 double                out_short,
                 double                mu,
                 int                   Nveg,
                 int                   veg,
                 int                   rec,
                 int                   gridcell,
                 int                   prcpdist,
                 char                  NEWCELL) {
/**********************************************************************
  write_debug		Keith Cherkauer		October 8, 1997

  This subroutine controls the output of selected debug data files.
  Which debugging files are output is controlled by flags set in the
  model control file.

  This subroutine is not essential to the operation of the model, and 
  can be excluded from the compiled code as long as references to it 
  from full_energy.c are removed.  Due to the number of static variables
  in this routine, it may contribut significantly to the size of the
  compiled code.

  Modifications:
  07-15-98 modified to work with elevation bands                 KAC

**********************************************************************/

  extern option_struct options;
  extern debug_struct debug;

  static short int   FIRST;
  static double    **MOIST_ERROR;
  static double     *INIT_MOIST;
  static double     *ENERGY_ERROR;
  static double     *ENERGY_ERROR_CALC;
  static double     *INFLOW;
  static double     *RUNOFF;
  static double     *BASEFLOW;
  static double     *EVAP;
  static double     *INSHORT;
  static double     *OUTSHORT;
  static double     *INLONG;
  static double     *OUTLONG;
  static double     *SENSIBLE;
  static double     *LATENT;
  static double     *GRND_FLUX;
  static double     *ADVECTION;
  static double     *DELTA_CC;
  static double     *SNOW_FLUX;
  static double     *REFREEZEENERGY;
  static double     *DELTA_H;

  int     i;
  int     Ntemp;
  int     band;
  int     Nbands;
  double  sum;
  double *Evap;
  double *curr_moist;
  double  curr_ice;
  double  Zsum;
  double  longwave;
  double  grnd_flux;
  double  advection;
  double  deltaH;
  double  deltaCC;
  double  Qnet;

  Nbands = options.SNOW_BAND;

  if(debug.PRT_FLUX && options.FULL_ENERGY) {

    /***** Record Hourly Energy Balance Terms *****/

    if(NEWCELL && prcpdist==0) {
      if(gridcell>0) {
        free((char *)ENERGY_ERROR);
        free((char *)ENERGY_ERROR_CALC);
        free((char *)INSHORT);
        free((char *)OUTSHORT);
        free((char *)INLONG);
        free((char *)OUTLONG);
        free((char *)SENSIBLE);
        free((char *)LATENT);
        free((char *)GRND_FLUX);
        free((char *)ADVECTION);
        free((char *)DELTA_H);
        free((char *)DELTA_CC);
        free((char *)SNOW_FLUX);
        free((char *)REFREEZEENERGY);
      }
      ENERGY_ERROR      = (double *)calloc(Nveg+1,sizeof(double));
      ENERGY_ERROR_CALC = (double *)calloc(Nveg+1,sizeof(double));
      INSHORT           = (double *)calloc(Nveg+1,sizeof(double));
      OUTSHORT          = (double *)calloc(Nveg+1,sizeof(double));
      INLONG            = (double *)calloc(Nveg+1,sizeof(double));
      OUTLONG           = (double *)calloc(Nveg+1,sizeof(double));
      SENSIBLE          = (double *)calloc(Nveg+1,sizeof(double));
      LATENT            = (double *)calloc(Nveg+1,sizeof(double));
      GRND_FLUX         = (double *)calloc(Nveg+1,sizeof(double));
      ADVECTION         = (double *)calloc(Nveg+1,sizeof(double));
      DELTA_H           = (double *)calloc(Nveg+1,sizeof(double));
      DELTA_CC          = (double *)calloc(Nveg+1,sizeof(double));
      SNOW_FLUX         = (double *)calloc(Nveg+1,sizeof(double));
      REFREEZEENERGY    = (double *)calloc(Nveg+1,sizeof(double));
    }
    if(rec==0 && prcpdist==0) {
      ENERGY_ERROR[veg]      = 0.;
      ENERGY_ERROR_CALC[veg] = 0.;
    }
    if(prcpdist==0) {
      INSHORT[veg]        = 0.;
      OUTSHORT[veg]       = 0.;
      INLONG[veg]         = 0.;
      OUTLONG[veg]        = 0.;
      SENSIBLE[veg]       = 0.;
      LATENT[veg]         = 0.;
      GRND_FLUX[veg]      = 0.;
      ADVECTION[veg]      = 0.;
      DELTA_H[veg]        = 0.;
      DELTA_CC[veg]       = 0.;
      SNOW_FLUX[veg]      = 0.;
      REFREEZEENERGY[veg] = 0.;
    } 

    for(band=0;band<Nbands;band++) {
      if(options.SNOW_MODEL) {
	deltaCC = energy[band].deltaCC;
	for(band=0;band<Nbands;band++)
	  Qnet = snow[band].Qnet / soil_con.AreaFract[band]; 
      }
      else {
	deltaCC = 0.;
	Qnet = 0.;
      }
      advection = energy[band].advection;
      deltaH = energy[band].deltaH;
      
      ENERGY_ERROR[veg] += energy[band].error * mu;
      ENERGY_ERROR_CALC[veg] 
	+= ((1.-energy[band].albedo)*energy[band].shortwave
	    + energy[band].longwave + energy[band].grnd_flux
	    + energy[band].latent + energy[band].sensible 
	    + energy[band].deltaH - energy[band].deltaCC 
	    - energy[band].snow_flux + energy[band].refreeze_energy 
	    + energy[band].advection) * mu / soil_con.AreaFract[band];
      
      INSHORT[veg]        += (1.-energy[band].albedo)*atmos.shortwave * mu 
	/ soil_con.AreaFract[band];
      INLONG[veg]         += atmos.longwave * mu / soil_con.AreaFract[band];
      SENSIBLE[veg]       += energy[band].sensible * mu 
	/ soil_con.AreaFract[band];
      LATENT[veg]         += energy[band].latent * mu 
	/ soil_con.AreaFract[band];
      GRND_FLUX[veg]      += energy[band].grnd_flux * mu 
	/ soil_con.AreaFract[band];
      ADVECTION[veg]      += energy[band].advection * mu 
	/ soil_con.AreaFract[band];
      DELTA_H[veg]        += energy[band].deltaH * mu 
	/ soil_con.AreaFract[band];
      DELTA_CC[veg]       += energy[band].deltaCC * mu 
	/ soil_con.AreaFract[band];
      SNOW_FLUX[veg]      += energy[band].snow_flux * mu 
	/ soil_con.AreaFract[band];
      REFREEZEENERGY[veg] += energy[band].refreeze_energy * mu 
	/ soil_con.AreaFract[band];
    }

    if(rec==0 && veg==0 && prcpdist==0) {
      fprintf(debug.fg_energy,"DATE\tNET SHT\tNET LNG\t");
      fprintf(debug.fg_energy,"GRND F\tLATENT\tSENSBL\tADVEC\tdel H\t");
      fprintf(debug.fg_energy,"del CC\tSNWFLX\tMELT\t");
      fprintf(debug.fg_energy,"ERROR\tERR CAL\tGRND T\tT_1\tWIND\n");
    }
    if((options.DIST_PRCP && prcpdist==1) || !options.DIST_PRCP) {
      fprintf(debug.fg_energy,"%7.4f\t%7.4lf\t%7.4lf",
	      (float)rec/24.0*(float)gp.dt, INSHORT[veg],
	      INLONG[veg]);
      fprintf(debug.fg_energy,"\t%7.4lf\t%7.4lf\t%7.4lf",
	      -GRND_FLUX[veg], LATENT[veg], SENSIBLE[veg]);
      fprintf(debug.fg_energy,"\t%7.4lf\t%7.4lf\t%7.4lf\t%7.4lf\t%7.4lf",
	      ADVECTION[veg], DELTA_H[veg], DELTA_CC[veg], SNOW_FLUX[veg], 
	      REFREEZEENERGY[veg]);
      fprintf(debug.fg_energy,"\t%7.4lf\t%7.4lf",
	      ENERGY_ERROR[veg],ENERGY_ERROR_CALC[veg]);
      fprintf(debug.fg_energy,"\t%7.4lf\t%7.4lf\t%7.4lf\n",
	      energy[0].T[0], energy[0].T[1],
	      atmos.wind);
    }
  }
 
  if(debug.PRT_SNOW) {

    /***** Record Hourly Snow Terms *****/
    
    for(band=0;band<Nbands;band++) {
      
      if(snow[band].snow) grnd_flux = energy->grnd_flux;
      else grnd_flux = 0.;
 
      if(rec==0 && veg==0 && prcpdist==0 && band==0) {
	/** Print File Header **/
	fprintf(debug.fg_snow,"Date\tBand\tSWE TOT\tSWE SRF\tSWE PCK\t");
	fprintf(debug.fg_snow,"GRND T\tlyr1 T\tSURF T\tPACK T\tMELT\t");
	fprintf(debug.fg_snow,"VPR FLX\tAIR T\tSNOW\tRAIN\tGRNDFLX\t");
	fprintf(debug.fg_snow,"DEPTH\tKAPPA\tCANOPY\tCNPYFLUX\n");
      }
      fprintf(debug.fg_snow,"%7.4f\t%7.4lf\t%7.4lf\t%7.4lf\t%7.3lf",
	      (float)rec/24.0*(float)gp.dt,snow[band].swq*1000.,
	      snow[band].surf_water*1000.,
	      snow[band].pack_water*1000.,energy->T[0]);
      fprintf(debug.fg_snow,"\t%7.3lf\t%7.3lf\t%7.3lf\t%7.3lf\t%7.3lf",
        energy->T[1],snow[band].surf_temp,snow[band].pack_temp,atmos.melt,
        snow[band].vapor_flux*1000.);
      fprintf(debug.fg_snow,"\t%7.3lf\t%7.3lf\t%7.4lf\t%7.4f",
	      atmos.air_temp,atmos.prec-atmos.rainonly,atmos.rainonly,
	      grnd_flux);
      fprintf(debug.fg_snow,"\t%7.4f\t%7.4f\t%7.4f\t%7.4f\n",
	      snow[band].depth,snow[band].density,
	      snow[band].snow_canopy*1000.,
	      snow[band].canopy_vapor_flux*1000.);
    }
  }
 
  if(debug.PRT_BALANCE) {
 
    /***** Compute Water Balance Error *****/

    if(NEWCELL && prcpdist==0) {
      if(gridcell>0) {
        for(i=0;i<=Nveg;i++) free((char *)MOIST_ERROR[i]);
        free((char *)INIT_MOIST);
        free((char *)MOIST_ERROR);
        free((char *)INFLOW);
        free((char *)RUNOFF);
        free((char *)BASEFLOW);
        free((char *)EVAP);
      }
      INIT_MOIST  = (double *)calloc(Nveg+1,sizeof(double));
      MOIST_ERROR = (double **)calloc(Nveg+1,sizeof(double*));
      INFLOW      = (double *)calloc(Nveg+1,sizeof(double));
      RUNOFF      = (double *)calloc(Nveg+1,sizeof(double));
      BASEFLOW    = (double *)calloc(Nveg+1,sizeof(double));
      EVAP        = (double *)calloc(Nveg+1,sizeof(double));
      for(i=0;i<=Nveg;i++)
        MOIST_ERROR[i] = (double *)calloc(options.Nlayer+3,sizeof(double));
    }
    if(rec==0 && prcpdist==0) {
      for(band=0;band<Nbands;band++) {
	INIT_MOIST[veg] = debug.store_moist[WET][band][options.Nlayer+2] * mu;
	INIT_MOIST[veg] += debug.store_moist[DRY][band][options.Nlayer+2] 
	  * (1. - mu);
      }
      for(i=0;i<options.Nlayer+3;i++) MOIST_ERROR[veg][i] = 0.;
    }
    if(prcpdist==0) {
      INFLOW[veg] = 0.;
      RUNOFF[veg] = 0.;
      BASEFLOW[veg] = 0.;
      EVAP[veg] = 0.;
    }
    Evap       = (double *)calloc(options.Nlayer+3,sizeof(double));
    curr_moist = (double *)calloc(options.Nlayer+3,sizeof(double));
    for(band=0;band<Nbands;band++) {
      Evap[options.Nlayer+2] = 0.;
      curr_moist[options.Nlayer+2] = 0.;
      if(veg<Nveg) {
	/** Vegetation Present **/
	Evap[0]       = veg_var[band].canopyevap;
	curr_moist[0] = veg_var[band].Wdew;
	if(options.FULL_ENERGY || options.SNOW_MODEL) { 
	  Evap[0]       += snow[band].canopy_vapor_flux * 1000.;
	  curr_moist[0] += (snow[band].snow_canopy) * 1000.;
	}
	Evap[options.Nlayer+2]       += Evap[0];
	curr_moist[options.Nlayer+2] += curr_moist[0];
      }
      else {
	/** No vegetation **/
	Evap[0]       = 0.;
	curr_moist[0] = 0.;
      }

      if(options.FULL_ENERGY || options.SNOW_MODEL) {
	/** Snow present **/
	Evap[1]                       = snow[band].vapor_flux * 1000.;
	Evap[options.Nlayer+2]       += Evap[1];
	curr_moist[1]                += snow[band].swq * 1000.;
	curr_moist[options.Nlayer+2] += curr_moist[1];
      }
      else {
	/** No snow **/
	Evap[1]       = 0.;
	curr_moist[1] = 0.;
      }

      for(i=0;i<options.Nlayer;i++) {
	/** All Soil Layers **/
	Evap[i+2]               = cell[band].layer[i].evap;
	Evap[options.Nlayer+2] += Evap[i+2];
	curr_moist[i+2] 
	  = find_total_layer_moisture(cell[band].layer[i],soil_con.depth[i]);
	curr_moist[options.Nlayer+2] += curr_moist[i+2];
      }

      /** Compute Moisture Balance Error **/
      for(i=0;i<options.Nlayer+3;i++) {
	if(prcpdist==0 && band==0) MOIST_ERROR[veg][i] = 0.;
	MOIST_ERROR[veg][i] 
	  += (debug.inflow[prcpdist][band][i] 
	      - (debug.outflow[prcpdist][band][i] + Evap[i]) 
	      - (curr_moist[i] - debug.store_moist[prcpdist][band][i])) 
	  * mu / soil_con.AreaFract[band];
	if(fabs(MOIST_ERROR[veg][i]) > 1.e-4) {
	  fprintf(stderr,"WARNING: Debug Layer %i has a Moisture Balance Error of %lf in rec %i, veg %i\n",i,MOIST_ERROR[veg][i],rec,veg);
	}
      }

      /** Store Variables **/
      INFLOW[veg]   += atmos.prec * soil_con.Pfactor[band] 
	/ soil_con.AreaFract[band];
      BASEFLOW[veg] += cell[band].baseflow * mu / soil_con.AreaFract[band];
      RUNOFF[veg]   += cell[band].runoff * mu / soil_con.AreaFract[band];
      EVAP[veg]     += Evap[options.Nlayer+2] * mu / soil_con.AreaFract[band];
    }

    if(rec==0 && veg==0 && prcpdist==0) {
      fprintf(debug.fg_balance,"Date\tVeg Num\tPrecip\tRunoff");
      fprintf(debug.fg_balance,"\tBFlow\tEvap\tdStor\tError\n");
    }
    if((options.DIST_PRCP && prcpdist==1) || !options.DIST_PRCP) {
      fprintf(debug.fg_balance,"%lf\t%i\t%lf\t%lf\t%lf",
	      (double)rec/24.0*(double)gp.dt,veg,INFLOW[veg],RUNOFF[veg],
	      BASEFLOW[veg]);
      fprintf(debug.fg_balance,"\t%lf\t%lf\t%lf\n",
          EVAP[veg],curr_moist[options.Nlayer+2]-INIT_MOIST[veg],
          MOIST_ERROR[veg][options.Nlayer+2]);
    }
    free((char*)Evap);
    free((char*)curr_moist);
  }
 
  if(debug.PRT_TEMP) {
  
    /***** Temperature Profile Debugging Output *****/
 
    if(rec==0 && veg==0 && prcpdist==0) {
      fprintf(debug.fg_temp,"%i\n",gp.Nnodes);
      fprintf(debug.fg_temp,"Date - hour(REC)\tveg\tband\tAir T\tFdpth\tTdpth");
      for(i=0;i<options.Nlayer;i++) 
	fprintf(debug.fg_temp,"\t%i Th T\t%i Fr T\t%i T",i,i,i);
      sum=0.0;
      fprintf(debug.fg_temp,"\t%.2lfcm",soil_con.depth[1]/2.*100.);
      for(i=0;i<gp.Nnodes;i++) {
	fprintf(debug.fg_temp,"\t%.0lf",sum*100.0);
	sum+=(energy->dz[i]+energy->dz[i+1])/2.0;
      }
      fprintf(debug.fg_temp,"\n");
    }
    for(band=0;band<Nbands;band++) {
      fprintf(debug.fg_temp,"%02i/%02i/%04i - %02i\t%7f\t%i\t%i",
	      dmy.day,dmy.month,dmy.year,dmy.hour,
	      (float)rec/24.0*(float)gp.dt,veg,band);
      fprintf(debug.fg_temp,"\t%6.2lf\t%6.2lf\t%6.2lf\t",
	      atmos.air_temp,energy->fdepth[0]*100.,energy->fdepth[1]*100.);
      
      if(options.FROZEN_SOIL) {
	Ntemp = 5;
	for(i=0;i<options.Nlayer;i++)
	  fprintf(debug.fg_temp,"%6.2lf\t%6.2lf\t%6.2lf\t",
		  cell[band].layer[i].T_thaw,cell[band].layer[i].T_froz,
		  cell[band].layer[i].T);
      }
      else Ntemp = 2;
      
      fprintf(debug.fg_temp,"%6.2lf",(energy->T[0]+energy->T[1])/2.);
      
      for(i=0;i<gp.Nnodes;i++) 
	fprintf(debug.fg_temp,"\t%6.2lf", energy->T[i]);
      fprintf(debug.fg_temp,"\n");
    }
 
  }
 
  if(debug.PRT_MOIST) {

    /***** Moisture Profile Debugging Output *****/

    if(FIRST != -999) {
      fprintf(debug.fg_moist,"Date - hour(REC)        \tVeg Num\tDist Num");
      fprintf(debug.fg_moist,"\tT Air");
      fprintf(debug.fg_moist,"\tInflow\tRunoff");
      for(i=0;i<options.Nlayer;i++)
        fprintf(debug.fg_moist,"\t%i Moist\t%i Ice",i,i);
      fprintf(debug.fg_moist,"\n");
    }

    fprintf(debug.fg_moist,"%02i/%02i/%04i - %02i\t%7f\t%i\t%i",
        dmy.day,dmy.month,dmy.year,dmy.hour,
        (float)rec/24.0*(float)gp.dt,veg,prcpdist);
    fprintf(debug.fg_moist,"\t%6.2lf\t%6.4lf\t%6.4lf",
        atmos.air_temp,cell[band].inflow,cell[band].runoff);
 
    curr_moist = (double *)calloc(1,sizeof(double));
    for(i=0;i<options.Nlayer;i++) {
      curr_moist[0] = cell[band].layer[i].moist_thaw
                 * cell[band].layer[i].tdepth / soil_con.depth[i];
      curr_moist[0] += cell[band].layer[i].moist_froz
                  * (cell[band].layer[i].fdepth - cell[band].layer[i].tdepth)
                  / soil_con.depth[i];
      curr_ice = cell[band].layer[i].ice * (cell[band].layer[i].fdepth
                - cell[band].layer[i].tdepth) / soil_con.depth[i];
      curr_moist[0] += cell[band].layer[i].moist * (soil_con.depth[i]
                  - cell[band].layer[i].fdepth) / soil_con.depth[i];
      fprintf(debug.fg_moist,"\t%6.4lf\t%6.4lf",
          curr_moist[0] / soil_con.depth[i] / 1000.,
          curr_ice / soil_con.depth[i] / 1000.);
    }
    fprintf(debug.fg_moist,"\n");
    free((char*)curr_moist);
  }
 
  if(debug.PRT_KAPPA) {
 
    /***** Soil Thermal Properties Profile Debugging Output *****/

    fprintf(debug.fg_kappa,"%02i/%02i/%04i - %02i\t%7f\t%i",
        dmy.day,dmy.month,dmy.year,dmy.hour,
        (float)rec/24.0*(float)gp.dt,veg);
 
    for(i=0;i<options.Nlayer;i++) {
      fprintf(debug.fg_kappa,"\t%6.2lf\t%6.0lf",
          cell[band].layer[i].kappa,
          cell[band].layer[i].Cs);
    }
    fprintf(debug.fg_kappa,"\n");
  }
 
  if(debug.PRT_GRID) {

    /***** Soil Themperature GMT Grided Profile Output *****/
 
    Zsum=0.;
    for(i=0;i<gp.Nnodes;i++) {
      fprintf(debug.fg_grid,"%7f\t%lf\t%lf\n",
          (float)rec/24.0*(float)gp.dt,
          Zsum,energy->T[i]);
      if(i<gp.Nnodes)
        Zsum+=(energy->dz[i]+energy->dz[i+1])/2.;
    }
 
  }

  FIRST = -999;
 
}
