#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vicNl.h>

void write_debug(atmos_data_struct atmos,
                 soil_con_struct soil_con,
                 cell_data_struct cell,
                 energy_bal_struct *energy,
                 snow_data_struct *snow,
                 veg_var_struct *veg_var,
                 dmy_struct dmy,
                 global_param_struct gp,
                 double out_short,
                 double mu,
                 int Nveg,
                 int veg,
                 int rec,
                 int gridcell,
                 int prcpdist,
                 char NEWCELL) {
/**********************************************************************
  write_debug		Keith Cherkauer		October 8, 1997

  This subroutine controls the output of selected debug data files.

**********************************************************************/

  extern option_struct options;
  extern debug_struct debug;

  static short int FIRST;
  static double **MOIST_ERROR;
  static double *INIT_MOIST;
  static double *ENERGY_ERROR;
  static double *ENERGY_ERROR_CALC;
  static double *INFLOW;
  static double *RUNOFF;
  static double *BASEFLOW;
  static double *EVAP;
  static double *INSHORT;
  static double *OUTSHORT;
  static double *INLONG;
  static double *OUTLONG;
  static double *SENSIBLE;
  static double *LATENT;
  static double *GRND_FLUX;
  static double *ADVECTION;
  static double *DELTA_CC;
  static double *QNET;
  static double *MELTENERGY;
  static double *DELTA_H;
  static double *SWQ;
  static double *SURF_WATER;
  static double *PACK_WATER;
  static double *SURF_TEMP;
  static double *TEMP_1;
  static double *SNOW_SURF_TEMP;
  static double *SNOW_PACK_TEMP;
  static double *SNOW_MELT;
  static double *VAPOR_FLUX;
  static double *CANOPY_VAPOR_FLUX;
  static double *RAINFALL;
  static double *SNOWFALL;
  static double *SNOW_DEPTH;
  static double *SNOW_DENSITY;
  static double *SNOW_CANOPY;

  int i;
  int Ntemp;
  double *Evap;
  double *curr_moist;
  double curr_ice;
  double Zsum;
  double longwave;
  double grnd_flux;
  double advection;
  double deltaH;
  double deltaCC;
  double meltenergy;
  double Qnet;

  if(debug.PRT_FLUX) {

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
        free((char *)QNET);
        free((char *)MELTENERGY);
      }
      ENERGY_ERROR = (double *)calloc(Nveg,sizeof(double));
      ENERGY_ERROR_CALC = (double *)calloc(Nveg,sizeof(double));
      INSHORT = (double *)calloc(Nveg,sizeof(double));
      OUTSHORT = (double *)calloc(Nveg,sizeof(double));
      INLONG = (double *)calloc(Nveg,sizeof(double));
      OUTLONG = (double *)calloc(Nveg,sizeof(double));
      SENSIBLE = (double *)calloc(Nveg,sizeof(double));
      LATENT = (double *)calloc(Nveg,sizeof(double));
      GRND_FLUX = (double *)calloc(Nveg,sizeof(double));
      ADVECTION = (double *)calloc(Nveg,sizeof(double));
      DELTA_H = (double *)calloc(Nveg,sizeof(double));
      DELTA_CC = (double *)calloc(Nveg,sizeof(double));
      QNET = (double *)calloc(Nveg,sizeof(double));
      MELTENERGY = (double *)calloc(Nveg,sizeof(double));
    }
    if(rec==0 && prcpdist==0) {
      ENERGY_ERROR[veg] = 0.;
      ENERGY_ERROR_CALC[veg] = 0.;
    }
    if(prcpdist==0) {
      INSHORT[veg] = 0.;
      OUTSHORT[veg] = 0.;
      INLONG[veg] = 0.;
      OUTLONG[veg] = 0.;
      SENSIBLE[veg] = 0.;
      LATENT[veg] = 0.;
      GRND_FLUX[veg] = 0.;
      ADVECTION[veg] = 0.;
      DELTA_H[veg] = 0.;
      DELTA_CC[veg] = 0.;
      QNET[veg] = 0.;
      MELTENERGY[veg] = 0.;
    } 

    if(snow->snow) {
      longwave = STEFAN_B*pow(snow->surf_temp + KELVIN,4.0);
      deltaCC = energy->deltaCC;
      meltenergy = snow->melt_energy;
      Qnet = snow->Qnet;
    }
    else {
      longwave = STEFAN_B*pow(energy->T[0] + KELVIN,4.0);
      deltaCC = 0.;
      meltenergy = 0.;
      Qnet = 0.;
    }
    advection = energy->advection;
    deltaH = energy->deltaH;

    ENERGY_ERROR[veg] += energy->error;
    ENERGY_ERROR_CALC[veg] += (1.-energy->albedo)*energy->shortwave
                            + energy->longwave + energy->grnd_flux
                            + energy->latent + energy->sensible 
                            + energy->deltaH - deltaCC + meltenergy;

    INSHORT[veg] += (1.-energy->albedo)*atmos.shortwave * mu;
    INLONG[veg] += atmos.longwave * mu;
    SENSIBLE[veg] += energy->sensible * mu;
    LATENT[veg] += energy->latent * mu;
    GRND_FLUX[veg] += energy->grnd_flux * mu;
    ADVECTION[veg] += advection * mu;
    DELTA_H[veg] += deltaH * mu;
    DELTA_CC[veg] += deltaCC * mu;
    QNET[veg] += Qnet * mu;
    MELTENERGY[veg] += meltenergy * mu;

    if(rec==0 && veg==0 && prcpdist==0) {
      fprintf(debug.fg_energy,"DATE\tNET SHT\tNET LNG\t");
      fprintf(debug.fg_energy,"GRND F\tLATENT\tSENSBL\tADVEC\tdel H\tdel CC\tMELT\t");
      fprintf(debug.fg_energy,"ERROR\tERR CAL\tGRND T\tT_1\tWIND\n");
    }
    if((options.DIST_PRCP && prcpdist==1) || !options.DIST_PRCP) {
      fprintf(debug.fg_energy,"%7.4f\t%7.2lf\t%7.2lf",
          (float)rec/24.0*(float)gp.dt, INSHORT[veg],
          INLONG[veg]);
      fprintf(debug.fg_energy,"\t%7.2lf\t%7.2lf\t%7.2lf",
          -GRND_FLUX[veg], LATENT[veg], SENSIBLE[veg]);
      fprintf(debug.fg_energy,"\t%7.2lf\t%7.2lf\t%7.2lf\t%7.2lf\t%7.2lf",
          ADVECTION[veg], DELTA_H[veg], DELTA_CC[veg], MELTENERGY[veg], QNET[veg]);
      fprintf(debug.fg_energy,"\t%7.2lf\t%7.2lf",
          ENERGY_ERROR[veg],ENERGY_ERROR_CALC[veg]);
      fprintf(debug.fg_energy,"\t%7.2lf\t%7.2lf\t%7.2lf\n",energy->T[0], energy->T[1],
          atmos.wind);
    }
  }
 
  if(debug.PRT_SNOW) {

    /***** Record Hourly Snow Terms *****/

    if(snow->snow) grnd_flux = energy->grnd_flux;
    else grnd_flux = 0.;
 
    fprintf(debug.fg_snow,"%7.4f\t%7.4lf\t%7.4lf\t%7.4lf\t%7.3lf",
        (float)rec/24.0*(float)gp.dt,snow->swq*1000.,snow->surf_water*1000.,
        snow->pack_water*1000.,energy->T[0]);
    fprintf(debug.fg_snow,"\t%7.3lf\t%7.3lf\t%7.3lf\t%7.3lf\t%7.3lf",
        energy->T[1],snow->surf_temp,snow->pack_temp,atmos.melt,
        snow->vapor_flux*1000.);
    fprintf(debug.fg_snow,"\t%7.3lf\t%7.3lf\t%7.4lf\t%7.4f",
        atmos.air_temp,atmos.prec-atmos.rainonly,atmos.rainonly,
        grnd_flux);
    fprintf(debug.fg_snow,"\t%7.4f\t%7.4f\t%7.4f\t%7.4f\n",
        snow->depth,snow->density,snow->snow_canopy*1000.,
        snow->canopy_vapor_flux*1000.);
  }
 
  if(debug.PRT_BALANCE) {
 
    /***** Compute Water Balance Error *****/

    if(NEWCELL && prcpdist==0) {
      if(gridcell>0) {
        for(i=0;i<Nveg;i++) free((char *)MOIST_ERROR[i]);
        free((char *)INIT_MOIST);
        free((char *)MOIST_ERROR);
        free((char *)INFLOW);
        free((char *)RUNOFF);
        free((char *)BASEFLOW);
        free((char *)EVAP);
      }
      INIT_MOIST = (double *)calloc(Nveg,sizeof(double));
      MOIST_ERROR = (double **)calloc(Nveg,sizeof(double*));
      INFLOW = (double *)calloc(Nveg,sizeof(double));
      RUNOFF = (double *)calloc(Nveg,sizeof(double));
      BASEFLOW = (double *)calloc(Nveg,sizeof(double));
      EVAP = (double *)calloc(Nveg,sizeof(double));
      for(i=0;i<Nveg;i++)
        MOIST_ERROR[i] = (double *)calloc(options.Nlayer+3,sizeof(double));
    }
    if(rec==0 && prcpdist==0) {
      INIT_MOIST[veg] = debug.store_moist[options.Nlayer+2];
      for(i=0;i<options.Nlayer+3;i++) MOIST_ERROR[veg][i] = 0.;
    }
    else if(rec==0) INIT_MOIST[veg] += debug.store_moist[options.Nlayer+2];
    if(prcpdist==0) {
      INFLOW[veg] = 0.;
      RUNOFF[veg] = 0.;
      BASEFLOW[veg] = 0.;
      EVAP[veg] = 0.;
    }
    Evap = (double *)calloc(options.Nlayer+3,sizeof(double));
    curr_moist = (double *)calloc(options.Nlayer+3,sizeof(double));
    Evap[options.Nlayer+2] = 0.;
    curr_moist[options.Nlayer+2] = 0.;
    if(veg<Nveg) {
      Evap[0] = veg_var->canopyevap;
      if(options.FULL_ENERGY || options.SNOW_MODEL)
        Evap[0] += snow->canopy_vapor_flux * 1000.;
      Evap[options.Nlayer+2] += Evap[0];
      curr_moist[0] = veg_var->Wdew;
      if(options.FULL_ENERGY || options.SNOW_MODEL)
        curr_moist[0] += snow->snow_canopy*1000.;
      curr_moist[options.Nlayer+2] += curr_moist[0];
    }
    else {
      Evap[0] = 0.;
      curr_moist[0] = 0.;
    }
    if(options.FULL_ENERGY || options.SNOW_MODEL) {
      Evap[1] = snow->vapor_flux * 1000.;
      Evap[options.Nlayer+2] += Evap[1];
      curr_moist[1] = snow->swq*1000.;
      curr_moist[options.Nlayer+2] += curr_moist[1];
    }
    else {
      Evap[1] = 0.;
      curr_moist[1] = 0.;
    }
    for(i=0;i<options.Nlayer;i++) {
      Evap[i+2] = cell.layer[i].evap;
      Evap[options.Nlayer+2] += Evap[i+2];
      curr_moist[i+2] = cell.layer[i].moist_thaw
                  * cell.layer[i].tdepth / soil_con.depth[i];
      curr_moist[i+2] += cell.layer[i].moist_froz
                  * (cell.layer[i].fdepth - cell.layer[i].tdepth)
                  / soil_con.depth[i];
      curr_moist[i+2] += cell.layer[i].ice * (cell.layer[i].fdepth
                  - cell.layer[i].tdepth) / soil_con.depth[i];
      curr_moist[i+2] += cell.layer[i].moist * (soil_con.depth[i]
                  - cell.layer[i].fdepth) / soil_con.depth[i];
      curr_moist[options.Nlayer+2] += curr_moist[i+2];
    }
    for(i=0;i<options.Nlayer+3;i++) {
      MOIST_ERROR[veg][i] += (debug.inflow[i] - (debug.outflow[i]
                          + Evap[i]) - (curr_moist[i]
                          - debug.store_moist[i])) * mu;
    }
    INFLOW[veg] += atmos.prec;
    BASEFLOW[veg] += cell.baseflow * mu;
    RUNOFF[veg] += cell.runoff * mu;
    EVAP[veg] += Evap[options.Nlayer+2] * mu;

    if(rec==0 && veg==0 && prcpdist==0) {
      fprintf(debug.fg_balance,"Date\tVeg Num\tPrecip\tRunoff");
      fprintf(debug.fg_balance,"\tBFlow\tEvap\tdStor\tError\n");
    }
    if((options.DIST_PRCP && prcpdist==1) || !options.DIST_PRCP) {
      fprintf(debug.fg_balance,"%lf\t%i\t%lf\t%lf\t%lf",
          (double)rec/24.0*(double)gp.dt,veg,INFLOW[veg],RUNOFF[veg],BASEFLOW);
      fprintf(debug.fg_balance,"\t%lf\t%lf\t%lf\n",
          EVAP[veg],curr_moist[options.Nlayer+2]-INIT_MOIST[veg],
          MOIST_ERROR[veg][options.Nlayer+2]);
    }
    free((char*)Evap);
    free((char*)curr_moist);
  }
 
  if(debug.PRT_TEMP) {
  
    /***** Temperature Profile Debugging Output *****/
 
    fprintf(debug.fg_temp,"%02i/%02i/%04i - %02i\t%7f\t%i",
        dmy.day,dmy.month,dmy.year,dmy.hour,
        (float)rec/24.0*(float)gp.dt,veg);
    fprintf(debug.fg_temp,"\t%6.2lf\t%6.2lf\t%6.2lf\t",
        atmos.air_temp,energy->fdepth[0]*100.,energy->fdepth[1]*100.);
 
    if(options.FROZEN_SOIL) {
      Ntemp = 5;
      for(i=0;i<options.Nlayer;i++)
        fprintf(debug.fg_temp,"%6.2lf\t%6.2lf\t%6.2lf\t",
            cell.layer[i].T_thaw,cell.layer[i].T_froz,
            cell.layer[i].T);
    }
    else Ntemp = 2;
 
    fprintf(debug.fg_temp,"%6.2lf",(energy->T[0]+energy->T[1])/2.);
 
    for(i=0;i<Ntemp;i++) fprintf(debug.fg_temp,"\t%6.2lf", energy->T[i]);
    fprintf(debug.fg_temp,"\n");
 
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
        atmos.air_temp,cell.inflow,cell.runoff);
 
    curr_moist = (double *)calloc(1,sizeof(double));
    for(i=0;i<options.Nlayer;i++) {
      curr_moist[0] = cell.layer[i].moist_thaw
                 * cell.layer[i].tdepth / soil_con.depth[i];
      curr_moist[0] += cell.layer[i].moist_froz
                  * (cell.layer[i].fdepth - cell.layer[i].tdepth)
                  / soil_con.depth[i];
      curr_ice = cell.layer[i].ice * (cell.layer[i].fdepth
                - cell.layer[i].tdepth) / soil_con.depth[i];
      curr_moist[0] += cell.layer[i].moist * (soil_con.depth[i]
                  - cell.layer[i].fdepth) / soil_con.depth[i];
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
          cell.layer[i].kappa,
          cell.layer[i].Cs);
    }
    fprintf(debug.fg_kappa,"\n");
  }
 
  if(debug.PRT_GRID) {

    /***** Soil Themperature GMT Grided Profile Output *****/
 
    Zsum=0.;
    for(i=0;i<gp.Ulayer+gp.Llayer+2;i++) {
      fprintf(debug.fg_grid,"%7f\t%lf\t%lf\n",
          (float)rec/24.0*(float)gp.dt,
          Zsum,energy->T[i]);
      if(i<gp.Ulayer+gp.Llayer+1)
        Zsum+=(energy->dz[i]+energy->dz[i+1])/2.;
    }
 
  }

  FIRST = -999;
 
}
