#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vicNl.h>

static char vcid[] = "$Id";

#if SAVE_STATE

void write_model_state(dist_prcp_struct    *prcp,
		       global_param_struct *gp,
		       int                  Nveg,
		       int                  cellnum,
		       outfiles_struct     *outfiles,
		       soil_con_struct     *soil_con) 
/*********************************************************************
  write_model_state      Keith Cherkauer           April 14, 2000

  This subroutine saves the model state at hour 0 of the date 
  defined in the global control file using STATEDAY, STATEMONTH,
  and STATEYEAR.  The saved files can then be used to initialize 
  the model to the same state as when the files were created.

  Soil moisture, soil thermal, and snowpack variables  are stored 
  for each vegetation type and snow band.  However moisture variables
  from the distributed precipitation model are averaged so that the
  model is restarted with mu = 1.

*********************************************************************/
{
  extern option_struct options;

  double tmpval;
  double Nsum;
  int    veg;
  int    band;
  int    lidx;
  int    nidx;
  int    dist;
  int    Ndist;
  int    Nbands;

  cell_data_struct     ***cell;
  snow_data_struct      **snow;
  energy_bal_struct     **energy;
  veg_var_struct       ***veg_var;

  if(options.DIST_PRCP) 
    Ndist = 2;
  else 
    Ndist = 1;
  Nbands = options.SNOW_BAND;

  cell    = prcp->cell;
  veg_var = prcp->veg_var;
  snow    = prcp->snow;
  energy  = prcp->energy;
  
  /* write cell information */
  fprintf(outfiles->statefile,"%i %i %i", cellnum, Nveg, Nbands);

  /* Write soil thermal node depths */
  Nsum = 0;
  for ( nidx = 0; nidx < options.Nnode; nidx++ ) {
    fprintf(outfiles->statefile," %f", Nsum);
    if ( nidx < options.Nnode - 1 )
      Nsum += (soil_con->dz_node[nidx] + soil_con->dz_node[nidx+1]) / 2.;
  }    
  fprintf(outfiles->statefile,"\n");

  /* Output for all vegetation types */
  for ( veg = 0; veg <= Nveg; veg++ ) {

    /* Output for all snow bands */
    for ( band = 0; band < Nbands; band++ ) {
      
      /* Write cell identification information */
      fprintf(outfiles->statefile,"%i %i", veg, band);
      
      /* Write average total soil moisture */
      for ( lidx = 0; lidx < options.Nlayer; lidx++ ) {
	tmpval = prcp->mu[veg] * (cell[WET][veg][band].layer[lidx].moist)
	  + (1. - prcp->mu[veg]) * (cell[DRY][veg][band].layer[lidx].moist);
	fprintf(outfiles->statefile," %f", tmpval);
      }
      
      /* Write average ice content */
      for ( lidx = 0; lidx < options.Nlayer; lidx++ ) {
	tmpval = prcp->mu[veg] * (cell[WET][veg][band].layer[lidx].ice)
	  + (1. - prcp->mu[veg]) * (cell[DRY][veg][band].layer[lidx].ice);
	fprintf(outfiles->statefile," %f", tmpval);
      }
      
      /* Write average dew storage */
      if ( veg < Nveg ) {
	tmpval = prcp->mu[veg] * (veg_var[WET][veg][band].Wdew)
	  + (1. - prcp->mu[veg]) * (veg_var[DRY][veg][band].Wdew);
	fprintf(outfiles->statefile," %f", tmpval);
      }
      
      /* Write snow data */
      fprintf(outfiles->statefile," %i %f %f %f %f %f", 
	      snow[veg][band].last_snow, snow[veg][band].swq, 
	      snow[veg][band].surf_temp, snow[veg][band].pack_temp, 
	      snow[veg][band].density, snow[veg][band].snow_canopy);
      
      /* Write soil thermal node temperatures */
      for ( nidx = 0; nidx < options.Nnode; nidx++ ) 
	fprintf(outfiles->statefile," %f", energy[veg][band].T[nidx]);
      
      fprintf(outfiles->statefile,"\n");
      
    }
  }

  /* Force file to be written */
  fflush(outfiles->statefile);

}

#endif
