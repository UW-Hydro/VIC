#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
#include <math.h>

void initialize_new_storm(dist_prcp_struct *prcp,
                          int Nveg,
                          int rec,
                          double old_mu,
                          double new_mu,
                          soil_con_struct soil_con,
                          int Tlayer) {
/**********************************************************************
  initialize_new_storm.c	Keith Cherkauer		January 13, 1998

  This subroutine redistributes soil moisture, and the snow pack (if
  present) before a new strom.

**********************************************************************/
 
  extern option_struct options;

  unsigned char error;
  char          ErrorString[MAXSTRING];
  int           veg;
  int           layer;
  int           i;
  int           dist;
  double        temp;
  double        temp_wet;
  double        temp_dry;
  double        frozen;
  double        unfrozen;
  double        Ltotal;
  double       *old_fdepth;
  double       *old_tdepth;
  double       *M;
  double       *Cs;
  double       *kappa;
  double      **layer_ice;

  for(veg=0;veg<=Nveg;veg++) {

    /** Redistribute Soil Moisture **/
    Ltotal = 0.;
    for(layer=0;layer<options.Nlayer;layer++) {
      Ltotal += soil_con.depth[layer];

      temp_wet = (prcp[0].dist[0].cell[veg].layer[layer].moist_thaw
          * prcp[0].dist[0].cell[veg].layer[layer].tdepth
          / soil_con.depth[layer]);
      temp_wet += (prcp[0].dist[0].cell[veg].layer[layer].moist_froz
          * (prcp[0].dist[0].cell[veg].layer[layer].fdepth
          - prcp[0].dist[0].cell[veg].layer[layer].tdepth)
          / soil_con.depth[layer]);
      temp_wet += (prcp[0].dist[0].cell[veg].layer[layer].moist
          * (soil_con.depth[layer]
          - prcp[0].dist[0].cell[veg].layer[layer].fdepth)
          / soil_con.depth[layer]);
      temp_wet += (prcp[0].dist[0].cell[veg].layer[layer].ice
          * (prcp[0].dist[0].cell[veg].layer[layer].fdepth
          - prcp[0].dist[0].cell[veg].layer[layer].tdepth)
          / soil_con.depth[layer]);

      temp_dry = (prcp[0].dist[1].cell[veg].layer[layer].moist_thaw
          * prcp[0].dist[1].cell[veg].layer[layer].tdepth
          / soil_con.depth[layer]);
      temp_dry += (prcp[0].dist[1].cell[veg].layer[layer].moist_froz
          * (prcp[0].dist[1].cell[veg].layer[layer].fdepth
          - prcp[0].dist[1].cell[veg].layer[layer].tdepth)
          / soil_con.depth[layer]);
      temp_dry += (prcp[0].dist[1].cell[veg].layer[layer].moist
          * (soil_con.depth[layer]
          - prcp[0].dist[1].cell[veg].layer[layer].fdepth)
          / soil_con.depth[layer]);
      temp_dry += (prcp[0].dist[1].cell[veg].layer[layer].ice
          * (prcp[0].dist[1].cell[veg].layer[layer].fdepth
          - prcp[0].dist[1].cell[veg].layer[layer].tdepth)
          / soil_con.depth[layer]);

 
      error = redistribute_new_storm(&temp_wet, &temp_dry, old_mu, new_mu);
      if(error) {
        sprintf(ErrorString,"Water balance does not balance after new storm: %lf -> %lf record %i\n",
            prcp[0].dist[0].cell[veg].layer[layer].moist*new_mu
            + prcp[0].dist[1].cell[veg].layer[layer].moist*(1.-new_mu),
            temp_wet+temp_dry,rec);
        vicerror(ErrorString);
      }

      prcp[0].dist[0].cell[veg].layer[layer].moist = temp_wet;
      prcp[0].dist[1].cell[veg].layer[layer].moist = temp_dry;
      prcp[0].dist[0].cell[veg].layer[layer].moist_thaw = 0.;
      prcp[0].dist[1].cell[veg].layer[layer].moist_thaw = 0.;
      prcp[0].dist[0].cell[veg].layer[layer].moist_froz = 0.;
      prcp[0].dist[1].cell[veg].layer[layer].moist_froz = 0.;
      prcp[0].dist[0].cell[veg].layer[layer].ice = 0.;
      prcp[0].dist[1].cell[veg].layer[layer].ice = 0.;

    }

    if(options.FULL_ENERGY || options.FROZEN_SOIL) {

      for(i=0;i<Tlayer;i++) {
        temp_wet = prcp[0].dist[0].energy[veg].T[i];
        temp_dry = prcp[0].dist[1].energy[veg].T[i];
        error = redistribute_new_storm(&temp_wet, &temp_dry, old_mu, new_mu);
        if(error) {
          sprintf(ErrorString,"Layer temperature does not balance after new storm: %lf -> %lf record %i, layer %i\n",
              prcp[0].dist[0].energy[veg].T[i]*new_mu
              + prcp[0].dist[1].energy[veg].T[i]*(1.-new_mu),
              temp_wet+temp_dry,rec,i);
          vicerror(ErrorString);
        }
        prcp[0].dist[0].energy[veg].T[i] = temp_wet;
        prcp[0].dist[1].energy[veg].T[i] = temp_dry;
      }

      /***** Find New Layer Temperatures *****/
      for(dist=0;dist<2;dist++) {

        if(options.FROZEN_SOIL && (prcp[0].dist[0].energy[veg].frozen
            || prcp[0].dist[1].energy[veg].frozen)) {

          /** Calculate New Layer Depths **/
          old_fdepth = (double *)calloc(options.Nlayer,sizeof(double));
          old_tdepth = (double *)calloc(options.Nlayer,sizeof(double));
          for(layer=0;layer<options.Nlayer;layer++) {
            old_fdepth[layer] = 0.;
            old_tdepth[layer] = 0.;
          }
          find_0_degree_fronts(prcp[0].dist[dist].energy,
              prcp[0].dist[dist].cell[veg].layer, Ltotal, soil_con.depth,
              prcp[0].dist[dist].energy[veg].T, options.Nlayer, Tlayer);

          /** Calculate New Layer Temperatures **/
          find_sublayer_temperatures(prcp[0].dist[dist].cell[veg].layer,
              prcp[0].dist[dist].energy[veg].T,
              prcp[0].dist[dist].energy[veg].dz,soil_con.depth,
              prcp[0].dist[dist].energy[veg].fdepth[0],
              prcp[0].dist[dist].energy[veg].fdepth[1],options.Nlayer,Tlayer);
  
          if(prcp[0].dist[dist].energy[veg].fdepth[0]>0)
            prcp[0].dist[dist].energy[veg].frozen = TRUE;
          else prcp[0].dist[dist].energy[veg].frozen = FALSE;

          /** Redistribute Soil Properties for New Frozen Soil Layer Size **/

          if(prcp[0].dist[dist].energy[veg].frozen)
            redistribute_moisture(prcp[0].dist[dist].cell[veg].layer,
                prcp[0].dist[dist].energy[veg].fdepth,
                soil_con.max_moist, old_fdepth, old_tdepth, soil_con.depth,
                options.Nlayer);
          free((char*)old_fdepth);
          free((char*)old_tdepth);

          /** Compute Amount of Unfrozen Moisture in Frozen Layer **/
          for(layer=0;layer<options.Nlayer;layer++) {
            if(prcp[0].dist[dist].cell[veg].layer[layer].fdepth > 0.0) {
              if(prcp[0].dist[dist].cell[veg].layer[layer].T_froz<0.
                  && prcp[0].dist[dist].cell[veg].layer[layer].T_froz != -999.) {
                unfrozen = maximum_unfrozen_water(
                           prcp[0].dist[dist].cell[veg].layer[layer].T_froz,
                           soil_con.max_moist[layer], soil_con.bubble,
                           soil_con.expt[layer]);
                if(unfrozen>soil_con.max_moist[layer] || unfrozen<0.)
                  unfrozen = soil_con.max_moist[layer];
                prcp[0].dist[dist].cell[veg].layer[layer].unfrozen = unfrozen;
  
                frozen = prcp[0].dist[dist].cell[veg].layer[layer].moist_froz
                       - unfrozen;
                if(frozen < 0.0) {
                  frozen = 0.0;
                  unfrozen = prcp[0].dist[dist].cell[veg].layer[layer].moist_froz;
                }
                prcp[0].dist[dist].cell[veg].layer[layer].ice = frozen;
                prcp[0].dist[dist].cell[veg].layer[layer].moist_froz = unfrozen;
              }
              else if(prcp[0].dist[dist].cell[veg].layer[layer].T_froz == 0.) {
                prcp[0].dist[dist].cell[veg].layer[layer].unfrozen
                    = soil_con.max_moist[layer];
                prcp[0].dist[dist].cell[veg].layer[layer].ice = 0.;
              }
              else if(prcp[0].dist[dist].cell[veg].layer[layer].T_froz != -999.) {
                sprintf(ErrorString,"ERROR: Frozen Layer Temperature > 0C (%lf)",
                    prcp[0].dist[dist].cell[veg].layer[layer].T_froz);
                vicerror(ErrorString);
              }
            }
            else prcp[0].dist[dist].cell[veg].layer[layer].ice=0.0;
          }

          layer_ice = (double **)calloc(options.Nlayer,sizeof(double *));
          for(i=0;i<options.Nlayer;i++) {
            layer_ice[i] = (double *)calloc(3,sizeof(double));
            layer_ice[i][0] = 0.;
            layer_ice[i][1] = prcp[0].dist[dist].cell[veg].layer[i].ice
                / (soil_con.depth[i] * 1000.);
            layer_ice[i][2] = 0.;
          }
          distribute_soil_property(prcp[0].dist[dist].energy[veg].dz,
              prcp[0].dist[dist].energy[veg].fdepth[0],
              prcp[0].dist[dist].energy[veg].fdepth[1],
              layer_ice,options.Nlayer,Tlayer,soil_con.depth,
              prcp[0].dist[dist].energy[veg].ice);
          for(i=0;i<options.Nlayer;i++) free((char *)layer_ice[i]);
          free((char *)layer_ice);
        }

        kappa = NULL;
        Cs = NULL;
        M = NULL;
        soil_thermal_calc(soil_con, prcp[0].dist[dist].cell[veg].layer,
            prcp[0].dist[dist].energy[veg], kappa,
            Cs, M, M, M, options.Nlayer, Tlayer);

        /** Save Thermal Conductivities for Energy Balance **/
        prcp[0].dist[dist].energy[veg].kappa[0]
            = prcp[0].dist[dist].cell[veg].layer[0].kappa;
        prcp[0].dist[dist].energy[veg].Cs[0]
            = prcp[0].dist[dist].cell[veg].layer[0].Cs;
        prcp[0].dist[dist].energy[veg].kappa[1]
            = prcp[0].dist[dist].cell[veg].layer[1].kappa;
        prcp[0].dist[dist].energy[veg].Cs[1]
            = prcp[0].dist[dist].cell[veg].layer[1].Cs;

      }
    }

    if(veg<Nveg) {
      temp_wet = prcp[0].dist[0].veg_var[veg].Wdew;
      temp_dry = prcp[0].dist[1].veg_var[veg].Wdew;
      error = redistribute_new_storm(&temp_wet, &temp_dry, old_mu, new_mu);
      if(error) {
        sprintf(ErrorString,"Wdew does not balance after new storm: %lf -> %lf record %i\n",
            prcp[0].dist[0].veg_var[veg].Wdew*new_mu
            + prcp[0].dist[1].veg_var[veg].Wdew*(1.-new_mu),
            temp_wet+temp_dry,rec);
        vicerror(ErrorString);
      }
      prcp[0].dist[0].veg_var[veg].Wdew = temp_wet;
      prcp[0].dist[1].veg_var[veg].Wdew = temp_dry;
    }

    if(options.FULL_ENERGY || options.SNOW_MODEL) {
  
      if(prcp[0].dist[0].snow[veg].snow || prcp[0].dist[1].snow[veg].snow) {
        temp_wet = (double)prcp[0].dist[0].snow[veg].last_snow;
        temp_dry = (double)prcp[0].dist[1].snow[veg].last_snow;
        error = redistribute_new_storm(&temp_wet, &temp_dry, old_mu, new_mu);
        prcp[0].dist[0].snow[veg].last_snow = (int)temp_wet;
        prcp[0].dist[1].snow[veg].last_snow = (int)temp_dry;

        temp_wet = prcp[0].dist[0].snow[veg].swq;
        temp_dry = prcp[0].dist[1].snow[veg].swq;
        error = redistribute_new_storm(&temp_wet, &temp_dry, old_mu, new_mu);
        if(error) {
          sprintf(ErrorString,"swq does not balance after new storm: %lf -> %lf record %i\n",
              prcp[0].dist[0].snow[veg].swq*new_mu
              + prcp[0].dist[1].snow[veg].swq*(1.-new_mu),
              temp_wet+temp_dry,rec);
          vicerror(ErrorString);
        }
        prcp[0].dist[0].snow[veg].swq = temp_wet;
        prcp[0].dist[1].snow[veg].swq = temp_dry;

        temp_wet = prcp[0].dist[0].snow[veg].surf_water;
        temp_dry = prcp[0].dist[1].snow[veg].surf_water;
        error = redistribute_new_storm(&temp_wet, &temp_dry, old_mu, new_mu);
        if(error) {
          sprintf(ErrorString,"surf_water does not balance after new storm: %lf -> %lf record %i\n",
              prcp[0].dist[0].snow[veg].surf_water*new_mu
              + prcp[0].dist[1].snow[veg].surf_water*(1.-new_mu),
              temp_wet+temp_dry,rec);
          vicerror(ErrorString);
        }
        prcp[0].dist[0].snow[veg].surf_water = temp_wet;
        prcp[0].dist[1].snow[veg].surf_water = temp_dry;

        temp_wet = prcp[0].dist[0].snow[veg].pack_water;
        temp_dry = prcp[0].dist[1].snow[veg].pack_water;
        error = redistribute_new_storm(&temp_wet, &temp_dry, old_mu, new_mu);
        if(error) {
          sprintf(ErrorString,"pack_water does not balance after new storm: %lf -> %lf record %i\n",
              prcp[0].dist[0].snow[veg].pack_water*new_mu
              + prcp[0].dist[1].snow[veg].pack_water*(1.-new_mu),
              temp_wet+temp_dry,rec);
          vicerror(ErrorString);
        }
        prcp[0].dist[0].snow[veg].pack_water = temp_wet;
        prcp[0].dist[1].snow[veg].pack_water = temp_dry;

        temp_wet = prcp[0].dist[0].snow[veg].density;
        temp_dry = prcp[0].dist[1].snow[veg].density;
        error = redistribute_new_storm(&temp_wet, &temp_dry, old_mu, new_mu);
        if(error) {
          sprintf(ErrorString,"density does not balance after new storm: %lf -> %lf record %i\n",
              prcp[0].dist[0].snow[veg].density*new_mu
              + prcp[0].dist[1].snow[veg].density*(1.-new_mu),
              temp_wet+temp_dry,rec);
          vicerror(ErrorString);
        }
        prcp[0].dist[0].snow[veg].density = temp_wet;
        prcp[0].dist[1].snow[veg].density = temp_dry;

        temp_wet = prcp[0].dist[0].snow[veg].depth;
        temp_dry = prcp[0].dist[1].snow[veg].depth;
        error = redistribute_new_storm(&temp_wet, &temp_dry, old_mu, new_mu);
        if(error) {
          sprintf(ErrorString,"depth does not balance after new storm: %lf -> %lf record %i\n",
              prcp[0].dist[0].snow[veg].depth*new_mu
              + prcp[0].dist[1].snow[veg].depth*(1.-new_mu),
              temp_wet+temp_dry,rec);
          vicerror(ErrorString);
        }
        prcp[0].dist[0].snow[veg].depth = temp_wet;
        prcp[0].dist[1].snow[veg].depth = temp_dry;

        temp_wet = prcp[0].dist[0].snow[veg].snow_canopy;
        temp_dry = prcp[0].dist[1].snow[veg].snow_canopy;
        error = redistribute_new_storm(&temp_wet, &temp_dry, old_mu, new_mu);
        if(error) {
          sprintf(ErrorString,"snow_canopy does not balance after new storm: %lf -> %lf record %i\n",
              prcp[0].dist[0].snow[veg].snow_canopy*new_mu
              + prcp[0].dist[1].snow[veg].snow_canopy*(1.-new_mu),
              temp_wet+temp_dry,rec);
          vicerror(ErrorString);
        }
        prcp[0].dist[0].snow[veg].snow_canopy = temp_wet;
        prcp[0].dist[1].snow[veg].snow_canopy = temp_dry;

        temp_wet = prcp[0].dist[0].snow[veg].tmp_int_storage;
        temp_dry = prcp[0].dist[1].snow[veg].tmp_int_storage;
        error = redistribute_new_storm(&temp_wet, &temp_dry, old_mu, new_mu);
        if(error) {
          sprintf(ErrorString,"tmp_int_storage does not balance after new storm: %lf -> %lf record %i\n",
              prcp[0].dist[0].snow[veg].tmp_int_storage*new_mu
              + prcp[0].dist[1].snow[veg].tmp_int_storage*(1.-new_mu),
              temp_wet+temp_dry,rec);
          vicerror(ErrorString);
        }
        prcp[0].dist[0].snow[veg].tmp_int_storage = temp_wet;
        prcp[0].dist[1].snow[veg].tmp_int_storage = temp_dry;

        if(prcp[0].dist[0].snow[veg].swq > 0) prcp[0].dist[0].snow[veg].snow = TRUE;
        else prcp[0].dist[0].snow[veg].snow = FALSE;
        if(prcp[0].dist[1].snow[veg].swq > 0) prcp[0].dist[1].snow[veg].snow = TRUE;
        else prcp[0].dist[1].snow[veg].snow = FALSE;

      }
    }
  }
}

unsigned char redistribute_new_storm(double *wet_value,
                                     double *dry_value,
                                     double old_mu,
                                     double new_mu) {
/**********************************************************************
  This subroutine redistributes the given parameter between wet and
  dry cell fractions when the precipitation changes in intensity.
**********************************************************************/

  unsigned char error;
  double temp_wet;
  double temp_dry;
  double diff;

  temp_wet = *wet_value * old_mu;
  temp_dry = *dry_value * (1. - old_mu);
  if(old_mu>new_mu && (new_mu!=1. && new_mu!=0.)) {
    *wet_value
        = (((1. - (old_mu-new_mu)/old_mu)*temp_wet) / new_mu);
    *dry_value
        = (((old_mu-new_mu)/old_mu*temp_wet + temp_dry) / (1.-new_mu));
  }
  else if(new_mu!=1. && new_mu!=0.) {
    *wet_value
        = (((new_mu-old_mu)/(1.-old_mu)*temp_dry + temp_wet) / new_mu);
    *dry_value
        = (((1. - (new_mu-old_mu)/(1.-old_mu))*temp_dry) / (1.-new_mu));
  }
  else {
    *wet_value
        = (temp_dry + temp_wet);
    *dry_value
        = (temp_dry + temp_wet);
  }
  diff = (temp_wet+temp_dry) - (*wet_value*new_mu + *dry_value*(1.-new_mu));
  if(fabs(diff) > 1.e-10) error = TRUE;
  else error = FALSE;

  return (error);
}
