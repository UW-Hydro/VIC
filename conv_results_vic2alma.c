#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

static char vcid[] = "$Id$";

void conv_results_vic2alma(out_data_struct *out_data, int dt, double *depth, out_data_alma_struct *out_data_alma, int rec)
/**********************************************************************
	conv_results_vic2alma	Ted Bohn	March 24, 2005

  This routine converts standard VIC output variables into ALMA-
  compliant output variables and stores them in the structure out_data_alma.
  
  modifications:

**********************************************************************/
{

  extern option_struct options;

  static double prev_SoilMoist[MAX_LAYERS];
  static double prev_SWE;
  static double prev_SurfStor;
  static double prev_Canopint;

  int index;
  double cumulative_depth;

  /* Energy Balance Variables */
  out_data_alma->SWnet = out_data->net_short;
  out_data_alma->LWnet = out_data->net_long;
  out_data_alma->Qle = out_data->latent;
  out_data_alma->Qh = out_data->sensible + out_data->advected_sensible[0];
  out_data_alma->Qg = out_data->grnd_flux;
  out_data_alma->Qf = out_data->fusion - out_data->refreeze_energy[0];
  out_data_alma->Qv = out_data->latent_sub[0];
  out_data_alma->Qa = out_data->advection[0];
  out_data_alma->DelSurfHeat = out_data->deltaH * (dt*3600);
  out_data_alma->DelColdCont = out_data->deltaCC[0] * (dt*3600);

  /* Water Balance Variables */
  out_data_alma->Rainf = out_data->rain;
  out_data_alma->Snowf = out_data->snow;
  out_data_alma->Evap = out_data->evap / (dt*3600);
  out_data_alma->Qs = out_data->runoff / (dt*3600);
  out_data_alma->Qsb = out_data->baseflow / (dt*3600);
  out_data_alma->Qsm = out_data->melt[0] / (dt*3600);

  /* Surface State Variables */
  out_data_alma->SnowT = out_data->snow_surf_temp[0] + KELVIN;
  out_data_alma->VegT = out_data->tfoliage;
  out_data_alma->BaresoilT = out_data->baresoilt;
  out_data_alma->AvgSurfT = out_data->AvgSurfT;
  out_data_alma->RadT = out_data->rad_temp;
  out_data_alma->Albedo = out_data->albedo;
  out_data_alma->SWE = out_data->swq[0];
  out_data_alma->SWEVeg = out_data->snow_canopy[0];
  out_data_alma->SurfStor = out_data->surfstor;

  /* Subsurface State Variables */
  cumulative_depth = 0.0;
  for (index=0; index<options.Nlayer; index++) {
    // depth[index] appears to really be the layer thickness
    cumulative_depth += depth[index];
    out_data_alma->HLayerDepth[index] = cumulative_depth;
    out_data_alma->SoilMoist[index] = out_data->moist[index] + out_data->ice[index];
    out_data_alma->LSoilMoist[index] = out_data->moist[index];
    out_data_alma->SoilTemp[index] = out_data->layer_temp[index] + KELVIN;
    out_data_alma->SMLiqFrac[index] = out_data_alma->LSoilMoist[index] / out_data_alma->SoilMoist[index];
    out_data_alma->SMFrozFrac[index] = 1 - out_data_alma->SMLiqFrac[index];
  }
  out_data_alma->SoilWet = out_data->soil_wetness;

  /* Evaporation Variables */
  out_data_alma->ECanop = out_data->evap_canop / (dt*3600);
  out_data_alma->TVeg = out_data->evap_veg / (dt*3600);
  out_data_alma->ESoil = out_data->evap_bare / (dt*3600);
  out_data_alma->EWater = 0.0;
#if LAKE_MODEL
  if (options.LAKES) {
    out_data_alma->EWater = out_data->evap_lake / (dt*3600);
  }
#endif
  out_data_alma->RootMoist = out_data->rootmoist;
  out_data_alma->Canopint = out_data->Wdew;
  out_data_alma->SubSnow = out_data->sub_total;
  out_data_alma->SubSurf = 0.0;
  if (out_data->aero_resist > SMALL) {
    out_data_alma->ACond = 1 / out_data->aero_resist;
  }
  else {
    out_data_alma->ACond = HUGE_RESIST;
  }

  /* Cold-Season Processes */
  out_data_alma->SnowFrac = out_data->coverage[0];
  out_data_alma->IceFrac = 0.0;
  out_data_alma->IceT = 0.0;
#if LAKE_MODEL
  if (options.LAKES) {
    out_data_alma->IceFrac = out_data->lake_ice_fract;
    out_data_alma->IceT = out_data->lake_ice_height / 100.;
  }
#endif
  for (index = 0; index<MAX_FRONTS; index++) {
    out_data_alma->Fdepth[index] = out_data->fdepth[index] / 100.;
    out_data_alma->Tdepth[index] = out_data->tdepth[index] / 100.;
  }
  out_data_alma->SAlbedo = out_data->snow_albedo;
  out_data_alma->SnowTProf = out_data->snow_pack_temp[0] + KELVIN;
  out_data_alma->SnowDepth = out_data->snow_depth[0] / 100.;

  /* More Water Balance Variables */
  /* these track change from last tstep to current */
  if (rec == 0) {
    out_data_alma->DelSoilMoist = 0;
    out_data_alma->DelSWE = 0;
    out_data_alma->DelSurfStor = 0;
    out_data_alma->DelIntercept = 0;
  }
  else {
    out_data_alma->DelSoilMoist = 0;
    for (index=0; index<options.Nlayer; index++) {
      out_data_alma->DelSoilMoist += ( out_data_alma->SoilMoist[index] - prev_SoilMoist[index] )
                                     * ( depth[index] / cumulative_depth );
    }
    out_data_alma->DelSWE = out_data_alma->SWE + out_data_alma->SWEVeg - prev_SWE;
    out_data_alma->DelSurfStor = out_data_alma->SurfStor - prev_SurfStor;
    out_data_alma->DelIntercept = out_data_alma->Canopint - prev_Canopint;
  }

  /* Save values for next time step */
  for (index=0; index<options.Nlayer; index++) {
    prev_SoilMoist[index] = out_data_alma->SoilMoist[index];
  }
  prev_SWE = out_data_alma->SWE + out_data_alma->SWEVeg;
  prev_SurfStor = out_data_alma->SurfStor;
  prev_Canopint = out_data_alma->Canopint;

}

