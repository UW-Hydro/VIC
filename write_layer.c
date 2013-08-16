
#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

static char vcid[] = "$Id: write_layer.c,v 4.1 2000/05/16 21:07:16 vicadmin Exp $";

void write_layer(layer_data_struct *layer,
                 int                veg,
                 int                Nlayer,
                 double            *depth)
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
  for(index=0;index<Nlayer;index++) printf("\t%f",layer[index].evap);
  printf("\n      Kappa:\t");
  for(index=0;index<Nlayer;index++) printf("\t%f",layer[index].kappa);
  printf("\n         Cs:\t");
  for(index=0;index<Nlayer;index++) printf("\t%f",layer[index].Cs);
  printf("\n\nMoisture Table\n---------------------------------------------------------------------------\n Moist:\t");
  for(index=0;index<Nlayer;index++) printf("\t%f",layer[index].moist);
  printf("\n        Ice:\t");
  for(index=0;index<Nlayer;index++) printf("\t%f",layer[index].ice);
  printf("\n---------------------------------------------------------------------------\nLayer Moist:\t");
  sum_moist = 0.;
  for(index=0;index<Nlayer;index++) {
    layer_moist = layer[index].moist;
    sum_moist += layer_moist;
    printf("\t%f",layer_moist);
  }
  printf("\n\n-----> Total Moisture = %f\n\n",sum_moist);
}




