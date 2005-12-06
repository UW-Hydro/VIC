#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vicNl.h>

static char vcid[] = "$Id$";


void write_model_state(dist_prcp_struct    *prcp,
		       global_param_struct *gp,
		       int                  Nveg,
		       int                  cellnum,
		       outfiles_struct     *outfiles,
		       soil_con_struct     *soil_con,
		       char                *STILL_STORM,
#if LAKE_MODEL
		       int                 *DRY_TIME,
		       lake_con_struct      lake_con)
#else
		       int                 *DRY_TIME)
#endif // LAKE_MODEL
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

  Modifications:
  04-10-03 Rewritten to handle updates to vicNl_def.h and to write
           the file as binary to minimize write time and differences
           with simulations started with the state file.         KAC
  04-10-03 Model is now restarted with the correct values for mu
           and LAST_STORM
  06-03-03 Modified to create ASCII as well as BINARY state file.  KAC
  06-06-03 It is not necessary to store the current ice content as
           it is recomputed in initialize_model_state.         KAC
  09-Oct-03 Removed initial space on veg/band info line in ASCII
	    file.                                                 TJB
  26-Oct-04 Changed calculation of Nbytes in binary state file to
	    account for bare soil values (extra veg class per grid
	    cell).  Without this fix, attempts to skip grid cells
	    fail.						TJB
  01-Nov-04 Added storage of state variables for SPATIAL_FROST and
	    LAKE_MODEL.						TJB
  02-Nov-04 Added a few more lake state variables.		TJB
  03-Nov-04 Now outputs extra_veg to aid other programs in parsing
	    state files.					TJB
  2005-11-29 SAVE_STATE is set in global param file not in user_def.h GCT
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
  int    byte, Nbytes;
  int    extra_veg;
#if SPATIAL_FROST
  int    frost_area;
#endif // SPATIAL_FROST

  cell_data_struct     ***cell;
  snow_data_struct      **snow;
  energy_bal_struct     **energy;
  veg_var_struct       ***veg_var;
#if LAKE_MODEL
  lake_var_struct         lake_var;
  int    node;
#endif // LAKE_MODEL

  if(options.DIST_PRCP) 
    Ndist = 2;
  else 
    Ndist = 1;
  Nbands = options.SNOW_BAND;

  cell    = prcp->cell;
  veg_var = prcp->veg_var;
  snow    = prcp->snow;
  energy  = prcp->energy;
#if LAKE_MODEL
  lake_var = prcp->lake_var;
#endif // LAKE_MODEL
 
  extra_veg = 0;
#if LAKE_MODEL
  if ( options.LAKES && lake_con.Cl[0] > 0 ) {
    extra_veg = 1; // add a veg type for the wetland
  }
#endif // LAKE_MODEL

  /* write cell information */
  if ( options.BINARY_STATE_FILE ) {
    fwrite( &cellnum, 1, sizeof(int), outfiles->statefile );
    fwrite( &Nveg, 1, sizeof(int), outfiles->statefile );
    fwrite( &extra_veg, 1, sizeof(int), outfiles->statefile );
    fwrite( &Nbands, 1, sizeof(int), outfiles->statefile );
  }
  else {
    fprintf( outfiles->statefile, "%i %i %i %i", cellnum, Nveg, extra_veg, Nbands );
  }
  // This stores the number of bytes from after this value to the end 
  // of the line.  DO NOT CHANGE unless you have changed the values
  // written to the state file.
  // IF YOU EDIT THIS FILE: UPDATE THIS VALUE!
  if ( options.BINARY_STATE_FILE ) {
    Nbytes = ( options.Nnode * sizeof(double) // dz_node
	       + (Nveg+1) * sizeof(double) // mu
	       + (Nveg+1) * sizeof(char) // STILL_STORM
	       + (Nveg+1) * sizeof(int) // DRY_TIME
	       + (Nveg+1) * Nbands * 2 * sizeof(int) // veg & band
	       + (Nveg+1) * Nbands * Ndist * options.Nlayer * sizeof(double) // soil moisture
#if SPATIAL_FROST
	       + (Nveg+1) * Nbands * Ndist * options.Nlayer * FROST_SUBAREAS * sizeof(double) // soil ice
#else
	       + (Nveg+1) * Nbands * Ndist * options.Nlayer * sizeof(double) // soil ice
#endif // SPATIAL_FROST
	       + Nveg * Nbands * Ndist * sizeof(double) // dew
	       + (Nveg+1) * Nbands * sizeof(int) // last_snow
	       + (Nveg+1) * Nbands * sizeof(char) // MELTING
	       + (Nveg+1) * Nbands * sizeof(double) * 9 // other snow parameters
	       + (Nveg+1) * Nbands * options.Nnode * sizeof(double) ); // soil temperatures
#if LAKE_MODEL
    if ( options.LAKES && lake_con.Cl[0] > 0 ) {
      Nbytes += sizeof(double) // wetland mu
		+ sizeof(char) // wetland STILL_STORM
		+ sizeof(int) // wetland DRY_TIME
		+ 2 * sizeof(int) // wetland veg and band (band = 0 for wetland)
		+ Ndist * options.Nlayer * sizeof(double) // wetland soil moisture
#if SPATIAL_FROST
	        + Ndist * options.Nlayer * FROST_SUBAREAS * sizeof(double) // wetland soil ice
#else
	        + Ndist * options.Nlayer * sizeof(double) // wetland soil ice
#endif // SPATIAL_FROST
		+ sizeof(double) // wetland dew
		+ sizeof(int) // wetland last_snow
		+ sizeof(char) // wetland MELTING
		+ 9 * sizeof(double) // wetland snow parameters
		+ options.Nnode * sizeof(double) // wetland soil temperatures
		+ sizeof(int) // numnod
		+ sizeof(double) // volume
		+ sizeof(double) // ldepth
		+ sizeof(double) // sarea
		+ sizeof(double) // dz
		+ sizeof(double) // surfdz
		+ lake_con.numnod * sizeof(double) // surface
		+ sizeof(double) // tp_in
		+ sizeof(double) // tempavg
		+ lake_con.numnod * sizeof(double) // temp
		+ lake_con.numnod * sizeof(double) // density
		+ sizeof(int) // mixmax
		+ sizeof(int) // activenod
		+ sizeof(double) // tempi
		+ sizeof(double) // hice
		+ sizeof(double) // fraci
		+ sizeof(double) // swe
		+ sizeof(double) // sdepth
      ;
    }
#endif // LAKE_MODEL
    fwrite( &Nbytes, 1, sizeof(int), outfiles->statefile );
  }

  /* Write soil thermal node depths */
  Nsum = 0;
  for ( nidx = 0; nidx < options.Nnode; nidx++ ) {
    if ( options.BINARY_STATE_FILE )
      fwrite( &soil_con->dz_node[nidx], 1, sizeof(double), 
	      outfiles->statefile );
    else
      fprintf( outfiles->statefile, " %f ", soil_con->dz_node[nidx] );
  }    
  if ( !options.BINARY_STATE_FILE )
    fprintf( outfiles->statefile, "\n" );

  /* Output for all vegetation types */
  for ( veg = 0; veg <= Nveg + extra_veg; veg++ ) {

    // Store distributed precipitation fraction
    if ( options.BINARY_STATE_FILE )
      fwrite( &prcp->mu[veg], 1, sizeof(double), outfiles->statefile );
    else
      fprintf( outfiles->statefile, "%f", prcp->mu[veg] );

    // Store distributed precipitation variables
    if ( options.BINARY_STATE_FILE ) {
      fwrite( &STILL_STORM[veg], 1, sizeof(char), outfiles->statefile );
      fwrite( &DRY_TIME[veg], 1, sizeof(int), outfiles->statefile );
    }
    else {
      fprintf( outfiles->statefile, " %i %i\n", (int)STILL_STORM[veg], 
	       DRY_TIME[veg] );
    }

#if LAKE_MODEL
  if ( options.LAKES && lake_con.Cl[0] > 0 && veg == Nveg + extra_veg ) {
    Nbands = 1; // wetland veg type only occurs in band 0
  }
#endif // LAKE_MODEL

    /* Output for all snow bands */
    for ( band = 0; band < Nbands; band++ ) {

      /* Write cell identification information */
      if ( options.BINARY_STATE_FILE ) {
	fwrite( &veg, 1, sizeof(int), outfiles->statefile );
	fwrite( &band, 1, sizeof(int), outfiles->statefile );
      }
      else {
	fprintf( outfiles->statefile, "%i %i", veg, band );
      }
      
      for ( dist = 0; dist < Ndist; dist ++ ) {
	// Store both wet and dry fractions if using distributed precipitation

	/* Write total soil moisture */
	for ( lidx = 0; lidx < options.Nlayer; lidx++ ) {
	  tmpval = cell[dist][veg][band].layer[lidx].moist;
	  if ( options.BINARY_STATE_FILE )
	    fwrite( &tmpval, 1, sizeof(double), outfiles->statefile );
	  else
	    fprintf( outfiles->statefile, " %f", tmpval );
	}

        /* Write average ice content */
        for ( lidx = 0; lidx < options.Nlayer; lidx++ ) {
#if SPATIAL_FROST
	  for ( frost_area = 0; frost_area < FROST_SUBAREAS; frost_area++ ) {
	    tmpval = cell[dist][veg][band].layer[lidx].ice[frost_area];
	    if ( options.BINARY_STATE_FILE ) {
	      fwrite( &tmpval, 1, sizeof(double), outfiles->statefile );
	    }
	    else {
	      fprintf( outfiles->statefile, " %f", tmpval );
	    }
	  }
#else
	  tmpval = cell[dist][veg][band].layer[lidx].ice;
	  if ( options.BINARY_STATE_FILE ) {
	    fwrite( &tmpval, 1, sizeof(double), outfiles->statefile );
	  }
	  else {
	    fprintf( outfiles->statefile, " %f", tmpval );
	  }
#endif // SPATIAL_FROST
        }

	/* Write dew storage */
	if ( veg < Nveg || ( veg == Nveg+extra_veg && extra_veg > 0 ) ) {
	  tmpval = veg_var[dist][veg][band].Wdew;
	  if ( options.BINARY_STATE_FILE )
	    fwrite( &tmpval, 1, sizeof(double), outfiles->statefile );
	  else
	    fprintf( outfiles->statefile, " %f", tmpval );
	}
      }
      
      /* Write snow data */
      if ( options.BINARY_STATE_FILE ) {
	fwrite( &snow[veg][band].last_snow, 1, sizeof(int), outfiles->statefile );
	fwrite( &snow[veg][band].MELTING, 1, sizeof(char), outfiles->statefile );
	fwrite( &snow[veg][band].coverage, 1, sizeof(double), outfiles->statefile );
	fwrite( &snow[veg][band].swq, 1, sizeof(double), outfiles->statefile );
	fwrite( &snow[veg][band].surf_temp, 1, sizeof(double), outfiles->statefile );
	fwrite( &snow[veg][band].surf_water, 1, sizeof(double), outfiles->statefile );
	fwrite( &snow[veg][band].pack_temp, 1, sizeof(double), outfiles->statefile );
	fwrite( &snow[veg][band].pack_water, 1, sizeof(double), outfiles->statefile );
	fwrite( &snow[veg][band].density, 1, sizeof(double), outfiles->statefile );
	fwrite( &snow[veg][band].coldcontent, 1, sizeof(double), outfiles->statefile );
	fwrite( &snow[veg][band].snow_canopy, 1, sizeof(double), outfiles->statefile );
      }
      else {
	fprintf( outfiles->statefile, " %i %i %f %f %f %f %f %f %f %f %f", 
		 snow[veg][band].last_snow, (int)snow[veg][band].MELTING, 
		 snow[veg][band].coverage, snow[veg][band].swq, 
		 snow[veg][band].surf_temp, snow[veg][band].surf_water, 
		 snow[veg][band].pack_temp, snow[veg][band].pack_water, 
		 snow[veg][band].density, snow[veg][band].coldcontent, 
		 snow[veg][band].snow_canopy );
      }
      
      /* Write soil thermal node temperatures */
      for ( nidx = 0; nidx < options.Nnode; nidx++ ) 
	if ( options.BINARY_STATE_FILE )
	  fwrite( &energy[veg][band].T[nidx], 1, sizeof(double), 
		  outfiles->statefile );
	else
	  fprintf( outfiles->statefile, " %f", energy[veg][band].T[nidx] );

      if ( !options.BINARY_STATE_FILE ) fprintf( outfiles->statefile, "\n" );
      
    }
  }

#if LAKE_MODEL
  if ( options.LAKES && lake_con.Cl[0] > 0 ) {
    if ( options.BINARY_STATE_FILE ) {
      fwrite( &lake_con.numnod, 1, sizeof(int), outfiles->statefile );
      fwrite( &lake_var.volume, 1, sizeof(double), outfiles->statefile );
      fwrite( &lake_var.ldepth, 1, sizeof(double), outfiles->statefile );
      fwrite( &lake_var.sarea, 1, sizeof(double), outfiles->statefile );
      fwrite( &lake_var.dz, 1, sizeof(double), outfiles->statefile );
      fwrite( &lake_var.surfdz, 1, sizeof(double), outfiles->statefile );
      for ( node = 0; node < lake_con.numnod; node++ ) {
        fwrite( &lake_var.surface[node], 1, sizeof(double), outfiles->statefile );
      }
      fwrite( &lake_var.tp_in, 1, sizeof(double), outfiles->statefile );
      fwrite( &lake_var.tempavg, 1, sizeof(double), outfiles->statefile );
      for ( node = 0; node < lake_con.numnod; node++ ) {
        fwrite( &lake_var.temp[node], 1, sizeof(double), outfiles->statefile );
      }
      for ( node = 0; node < lake_con.numnod; node++ ) {
        fwrite( &lake_var.density[node], 1, sizeof(double), outfiles->statefile );
      }
      fwrite( &lake_var.mixmax, 1, sizeof(int), outfiles->statefile );
      fwrite( &lake_var.activenod, 1, sizeof(int), outfiles->statefile );
      fwrite( &lake_var.tempi, 1, sizeof(double), outfiles->statefile );
      fwrite( &lake_var.hice, 1, sizeof(double), outfiles->statefile );
      fwrite( &lake_var.fraci, 1, sizeof(double), outfiles->statefile );
      fwrite( &lake_var.swe, 1, sizeof(double), outfiles->statefile );
      fwrite( &lake_var.sdepth, 1, sizeof(double), outfiles->statefile );
    }
    else {
      fprintf( outfiles->statefile, "%d", lake_con.numnod );
      fprintf( outfiles->statefile, " %f", lake_var.volume );
      fprintf( outfiles->statefile, " %f", lake_var.ldepth );
      fprintf( outfiles->statefile, " %f", lake_var.sarea );
      fprintf( outfiles->statefile, " %f", lake_var.dz );
      fprintf( outfiles->statefile, " %f", lake_var.surfdz );
      for ( node = 0; node < lake_con.numnod; node++ ) {
        fprintf( outfiles->statefile, " %f", lake_var.surface[node] );
      }
      fprintf( outfiles->statefile, " %f", lake_var.tp_in );
      fprintf( outfiles->statefile, " %f", lake_var.tempavg );
      for ( node = 0; node < lake_con.numnod; node++ ) {
        fprintf( outfiles->statefile, " %f", lake_var.temp[node] );
      }
      for ( node = 0; node < lake_con.numnod; node++ ) {
        fprintf( outfiles->statefile, " %f", lake_var.density[node] );
      }
      fprintf( outfiles->statefile, " %d", lake_var.mixmax );
      fprintf( outfiles->statefile, " %d", lake_var.activenod );
      fprintf( outfiles->statefile, " %f", lake_var.tempi );
      fprintf( outfiles->statefile, " %f", lake_var.hice );
      fprintf( outfiles->statefile, " %f", lake_var.fraci );
      fprintf( outfiles->statefile, " %f", lake_var.swe );
      fprintf( outfiles->statefile, " %f", lake_var.sdepth );
      fprintf( outfiles->statefile, "\n" );
    }
  }
#endif // LAKE_MODEL
  /* Force file to be written */
  fflush(outfiles->statefile);

}

