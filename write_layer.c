
#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

static char vcid[] = "$Id$";

void write_layer(layer_data_struct *layer,
                 int veg,
                 int Nlayer,
                 double *depth)
/**********************************************************************
	write_soilvar		Keith Cherkauer		July 17, 1997

  This routine writes soil variables to stdout.  It creates a table 
  of soil moisture values which shows how much liquid water and ice 
  are contained in the thawed, frozen and unfrozen sublayers of each 
  soil layer.  It also gives the total soil moisture for each layer.

**********************************************************************/
{
  extern option_struct options;

  int index;
  double layer_moist;
  double sum_moist;

  printf("Layer Data for Vegetation Type #%i\n",veg);
  printf("Layer:\t");
  for(index=0;index<Nlayer;index++) printf("\t\t%i",index+1);
  printf("\nEvaporation:\t");
  for(index=0;index<Nlayer;index++) printf("\t%lf",layer[index].evap);
  printf("\n      Kappa:\t");
  for(index=0;index<Nlayer;index++) printf("\t%lf",layer[index].kappa);
  printf("\n         Cs:\t");
  for(index=0;index<Nlayer;index++) printf("\t%lf",layer[index].Cs);
  printf("\n\nMoisture Table\n---------------------------------------------------------------------------\n Thaw Moist:\t");
  for(index=0;index<Nlayer;index++) printf("\t%lf",layer[index].moist_thaw);
  printf("\n Froz Moist:\t");
  for(index=0;index<Nlayer;index++) printf("\t%lf",layer[index].moist_froz);
  printf("\n      Moist:\t");
  for(index=0;index<Nlayer;index++) printf("\t%lf",layer[index].moist);
  printf("\n        Ice:\t");
  for(index=0;index<Nlayer;index++) printf("\t%lf",layer[index].ice);
  printf("\n Thaw Depth:\t");
  for(index=0;index<Nlayer;index++) printf("\t%lf",layer[index].tdepth);
  printf("\n Froz Depth:\t");
  for(index=0;index<Nlayer;index++) printf("\t%lf",layer[index].fdepth);
  printf("\n---------------------------------------------------------------------------\nLayer Moist:\t");
  sum_moist = 0.;
  for(index=0;index<Nlayer;index++) {
    layer_moist = layer[index].moist_thaw * layer[index].tdepth / depth[index];
    layer_moist += layer[index].moist_froz * (layer[index].fdepth 
        - layer[index].tdepth) / depth[index];
    layer_moist += layer[index].ice * (layer[index].fdepth 
        - layer[index].tdepth) / depth[index];
    layer_moist += layer[index].moist * (depth[index] - layer[index].fdepth)
        / depth[index];
    sum_moist += layer_moist;
    printf("\t%lf",layer_moist);
  }
  printf("\n\n-----> Total Moisture = %lf\n\n",sum_moist);
}




