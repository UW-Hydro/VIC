#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vicNl.h>

static char vcid[] = "$Id$";

void read_initial_model_state(FILE                *statefile,
			      dist_prcp_struct    *prcp,
			      global_param_struct *gp,
			      int                  Nveg,
			      int                  Nbands,
			      int                  cellnum,
			      soil_con_struct     *soil_con,
			      int                  Ndist,
			      char                *init_STILL_STORM,
			      int                 *init_DRY_TIME) 
/*********************************************************************
  read_initial_model_state   Keith Cherkauer         April 14, 2000

  This subroutine initializes the model state at hour 0 of the date 
  defined in the given state file.  

  Soil moisture, soil thermal, and snowpack variables  are stored 
  for each vegetation type and snow band.  However moisture variables
  from the distributed precipitation model are averaged so that the
  model is restarted with mu = 1.

  Modifications:
  04-10-03 Rewritten to handle updates to vicNl_def.h and to write
           the file as binary to minimize write time and differences
           with simulations started with the state file.         KAC
  04-10-03 Modified to read storm parameters from the state file.  KAC
  06-03-03 Modified to read ASCII as well as BINARY state file.  KAC
  06-06-03 It is not necessary to initialize ice content from the
           model state file as the model recomutes it in subsequent
           steps.                                               KAC
  06-06-03 Added check to make sure that soil moisture obtained from
           the state file does not exceed the maximum moisture 
           content.                                             KAC
  06-07-03 Added check to verify that the sum of the defined nodes
           equals the damping depth.                            KAC
  03-Oct-03 Modified to loop over tmp_Nveg and tmp_Nband when searching
            for desired cellnum in ASCII file, rather than over Nveg
            and Nbands.  As we skip over other records in the state
            file while searching for the desired record, the loop
            must parse each undesired record differently, according
            to how many veg classes and snow bands exist in the
            record (tmp_Nveg and tmp_Nband, respectively), rather
            than the number of veg classes and snow bands in the
            desired record (Nveg and Nbands, respectively).       TJB

*********************************************************************/
{
  extern option_struct options;

  char   tmpstr[MAXSTRING];
  char   ErrStr[MAXSTRING];
  char   tmpchar;
  double tmpval;
  double Nsum;
  double depth_node[MAX_NODES];
  double sum;
  int    veg, iveg;
  int    band, iband;
  int    lidx;
  int    nidx;
  int    dist;
  int    tmp_cellnum;
  int    tmp_Nveg;
  int    tmp_Nband;
  int    tmp_char;
  int    byte, Nbytes;
#if SPATIAL_FROST
  int    frost_area;
#endif

  cell_data_struct     ***cell;
  snow_data_struct      **snow;
  energy_bal_struct     **energy;
  veg_var_struct       ***veg_var;

  cell    = prcp->cell;
  veg_var = prcp->veg_var;
  snow    = prcp->snow;
  energy  = prcp->energy;
  
#if !NO_REWIND 
  rewind(statefile);
  
  /* skip header */
  if ( options.BINARY_STATE_FILE ) 
    fread(&tmpstr, 1, sizeof(int)*5, statefile);
  else {
    fgets(tmpstr, MAXSTRING, statefile);
    fgets(tmpstr, MAXSTRING, statefile);
  }
#endif

  /* read cell information */
  if ( options.BINARY_STATE_FILE ) {
    fread( &tmp_cellnum, 1, sizeof(int), statefile );
    fread( &tmp_Nveg, 1, sizeof(int), statefile );
    fread( &tmp_Nband, 1, sizeof(int), statefile );
    fread( &Nbytes, 1, sizeof(int), statefile );
  }
  else {
    fscanf( statefile, "%i %i %i", &tmp_cellnum, &tmp_Nveg, &tmp_Nband );
  }
  
  // Skip over unused cell information
  while ( tmp_cellnum != cellnum && !feof(statefile) ) {
    if ( options.BINARY_STATE_FILE ) {
      // skip rest of current cells info
      for ( byte = 0; byte < Nbytes; byte++ ) 
	fread ( &tmpchar, 1, 1, statefile);
      // read info for next cell
      fread( &tmp_cellnum, 1, sizeof(int), statefile );
      fread( &tmp_Nveg, 1, sizeof(int), statefile );
      fread( &tmp_Nband, 1, sizeof(int), statefile );
      fread( &Nbytes, 1, sizeof(int), statefile );
    }
    else {
      // skip rest of current cells info
      fgets(tmpstr, MAXSTRING, statefile); // skip rest of general cell info
      for ( veg = 0; veg <= tmp_Nveg; veg++ ) {
	fgets(tmpstr, MAXSTRING, statefile); // skip dist precip info
	for ( band = 0; band < tmp_Nband; band++ )
	  fgets(tmpstr, MAXSTRING, statefile); // skip snowband info
      }
      // read info for next cell
      fscanf( statefile, "%i %i %i", &tmp_cellnum, &tmp_Nveg, &tmp_Nband );
    }
  }

  if ( feof(statefile) ) {
    sprintf(ErrStr, "Requested grid cell (%i) is not in the model state file.", 
	    cellnum);
    nrerror(ErrStr);
  }

  if ( tmp_Nveg != Nveg ) {
    sprintf(ErrStr,"The number of vegetation types in cell %i (%i) does not equal that defined in vegetation parameter file (%i).  Check your input files.", cellnum, tmp_Nveg, Nveg);
    nrerror(ErrStr);
  }
  if ( tmp_Nband != Nbands ) {
    sprintf(ErrStr,"The number of snow bands in cell %i (%i) does not equal that defined in the snow band file (%i).  Check your input files.", cellnum, tmp_Nband, Nbands);
    nrerror(ErrStr);
  }

  /* Read soil thermal node depths */
  sum = 0;
  for ( nidx = 0; nidx < options.Nnode; nidx++ ) {
    if ( options.BINARY_STATE_FILE ) 
      fread( &soil_con->dz_node[nidx], 1, sizeof(double), statefile );
    else 
      fscanf( statefile, "%lf", &soil_con->dz_node[nidx] );
    sum += soil_con->dz_node[nidx];
  }
  if ( options.Nnode == 1 ) soil_con->dz_node[0] = 0;
  sum -= 0.5 * soil_con->dz_node[0];
  sum -= 0.5 * soil_con->dz_node[options.Nnode-1];
  if ( abs( sum - soil_con->dp ) > SMALL ) {
    fprintf( stderr, "WARNING: Sum of soil nodes (%f) exceeds defined damping depth (%f).  Resetting damping depth.\n", sum, soil_con->dp );
    soil_con->dp = sum;
  }
    
  /* Input for all vegetation types */
  for ( veg = 0; veg <= Nveg; veg++ ) {

    // read distributed precipitation variables
    if ( options.BINARY_STATE_FILE ) {
      fread( &prcp->mu[veg], 1, sizeof(double), statefile );
      fread( &init_STILL_STORM[veg], 1, sizeof(char), statefile );
      fread( &init_DRY_TIME[veg], 1, sizeof(int), statefile );
    }
    else {
      fscanf( statefile, "%lf %i %i", &prcp->mu[veg], &tmp_char, 
	      &init_DRY_TIME[veg] );
      init_STILL_STORM[veg] = (char)tmp_char;
    }

    /* Input for all snow bands */
    for ( band = 0; band < Nbands; band++ ) {
      
      /* Read cell identification information */
      if ( options.BINARY_STATE_FILE ) {
	if ( fread( &iveg, 1, sizeof(int), statefile) != sizeof(int) ) 
	  nrerror("End of model state file found unexpectedly");
	if ( fread( &iband, 1, sizeof(int), statefile) != sizeof(int) ) 
	  nrerror("End of model state file found unexpectedly");
      }
      else {
	if ( fscanf(statefile,"%i %i", &iveg, &iband) == EOF ) 
	  nrerror("End of model state file found unexpectedly");
      }
      if ( iveg != veg || iband != band ) {
	fprintf(stderr,"The vegetation and snow band indices in the model state file (veg = %i, band = %i) do not match those currently requested (veg = %i , band = %i).  Model state file must be stored with variables for all vegetation indexed by variables for all snow bands.", iveg, iband, veg, band);
	nrerror(ErrStr);
      }
      
      // Read both wet and dry fractions if using distributed precipitation
      for ( dist = 0; dist < Ndist; dist ++ ) {

	/* Read total soil moisture */
	for ( lidx = 0; lidx < options.Nlayer; lidx++ ) {
	  if ( options.BINARY_STATE_FILE ) {
	    if ( fread( &cell[dist][veg][band].layer[lidx].moist, 1, 
			sizeof(double), statefile ) != sizeof(double) )
	      nrerror("End of model state file found unexpectedly");
	  }
	  else {
	    if ( fscanf(statefile," %lf", 
			&cell[WET][veg][band].layer[lidx].moist) == EOF ) 
	      nrerror("End of model state file found unexpectedly");
	  }
	  if ( cell[WET][veg][band].layer[lidx].moist > soil_con->max_moist[lidx] ) {
	    // Check that soil moisture does not exceed maximum allowed
	    fprintf( stderr, "WARNING: Maximum soil moisture exceeded in layer %i for veg type %i and snow band %i.  Value reset to maximum (%f mm).\n", lidx, veg, band, soil_con->max_moist[lidx] );
	    cell[WET][veg][band].layer[lidx].moist = soil_con->max_moist[lidx];
	  }
	}
      
	/* Read dew storage */
	if ( veg < Nveg ) {
	  if ( options.BINARY_STATE_FILE ) {
	    if ( fread( &veg_var[dist][veg][band].Wdew, 1, sizeof(double), 
			statefile ) != sizeof(double) ) 
	      nrerror("End of model state file found unexpectedly");
	  }
	  else {
	    if ( fscanf(statefile," %lf", &veg_var[WET][veg][band].Wdew) == EOF ) 
	      nrerror("End of model state file found unexpectedly");
	  }
	}
      }
      
      /* Read snow data */
      if ( options.BINARY_STATE_FILE ) {
	if ( fread( &snow[veg][band].last_snow, 1, sizeof(int), 
		    statefile ) != sizeof(int) )
	  nrerror("End of model state file found unexpectedly");
	if ( fread( &snow[veg][band].MELTING, 1, sizeof(char), 
		    statefile ) != sizeof(char) )
	  nrerror("End of model state file found unexpectedly");
	if ( fread( &snow[veg][band].coverage, 1, sizeof(double), 
		    statefile ) != sizeof(double) )
	  nrerror("End of model state file found unexpectedly");
	if ( fread( &snow[veg][band].swq, 1, sizeof(double), 
		    statefile ) != sizeof(double) )
	  nrerror("End of model state file found unexpectedly");
	if ( fread( &snow[veg][band].surf_temp, 1, sizeof(double), 
		    statefile ) != sizeof(double) )
	  nrerror("End of model state file found unexpectedly");
	if ( fread( &snow[veg][band].surf_water, 1, sizeof(double), 
		    statefile ) != sizeof(double) )
	  nrerror("End of model state file found unexpectedly");
	if ( fread( &snow[veg][band].pack_temp, 1, sizeof(double), 
		    statefile ) != sizeof(double) )
	  nrerror("End of model state file found unexpectedly");
	if ( fread( &snow[veg][band].pack_water, 1, sizeof(double), 
		    statefile ) != sizeof(double) )
	  nrerror("End of model state file found unexpectedly");
	if ( fread( &snow[veg][band].density, 1, sizeof(double), 
		    statefile ) != sizeof(double) )
	  nrerror("End of model state file found unexpectedly");
	if ( fread( &snow[veg][band].coldcontent, 1, sizeof(double), 
		    statefile ) != sizeof(double) )
	  nrerror("End of model state file found unexpectedly");
	if ( fread( &snow[veg][band].snow_canopy, 1, sizeof(double), 
		    statefile ) != sizeof(double) )
	  nrerror("End of model state file found unexpectedly");
      }
      else {
	if ( fscanf(statefile," %i %i %lf %lf %lf %lf %lf %lf %lf %lf %lf", 
		    &snow[veg][band].last_snow, &tmp_char,
		    &snow[veg][band].coverage, &snow[veg][band].swq, 
		    &snow[veg][band].surf_temp, &snow[veg][band].surf_water, 
		    &snow[veg][band].pack_temp, &snow[veg][band].pack_water, 
		    &snow[veg][band].density, &snow[veg][band].coldcontent, 
		    &snow[veg][band].snow_canopy) 
	     == EOF ) 
	  nrerror("End of model state file found unexpectedly");
	snow[veg][band].MELTING = (char)tmp_char;
      }
      if(snow[veg][band].density > 0.) 
	snow[veg][band].depth = 1000. * snow[veg][band].swq 
	  / snow[veg][band].density;
      
      /* Read soil thermal node temperatures */
      for ( nidx = 0; nidx < options.Nnode; nidx++ ) {
	if ( options.BINARY_STATE_FILE ) {
	  if ( fread( &energy[veg][band].T[nidx], 1, sizeof(double), 
		      statefile ) != sizeof(double) )
	    nrerror("End of model state file found unexpectedly");
	}
	else {
	  if ( fscanf(statefile," %lf", &energy[veg][band].T[nidx]) == EOF )
	    nrerror("End of model state file found unexpectedly");
	}
      }
    }
  }
}
