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

*********************************************************************/
{
  extern option_struct options;

  char   tmpstr[MAXSTRING];
  char   ErrStr[MAXSTRING];
  char   tmpchar;
  double tmpval;
  double Nsum;
  double depth_node[MAX_NODES];
  int    veg, iveg;
  int    band, iband;
  int    lidx;
  int    nidx;
  int    dist;
  int    tmp_cellnum;
  int    tmp_Nveg;
  int    tmp_Nband;
  int    byte, Nbytes;

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
  //fgets(tmpstr, MAXSTRING, statefile);
  //fgets(tmpstr, MAXSTRING, statefile);
  fread(&tmpstr, 1, sizeof(int)*5, statefile);
#endif

  /* read cell information */
  //fscanf(statefile, "%i %i %i", &tmp_cellnum, &tmp_Nveg, &tmp_Nband);
  fread( &tmp_cellnum, 1, sizeof(int), statefile );
  fread( &tmp_Nveg, 1, sizeof(int), statefile );
  fread( &tmp_Nband, 1, sizeof(int), statefile );
  fread( &Nbytes, 1, sizeof(int), statefile );
  
  while ( tmp_cellnum != cellnum && !feof(statefile) ) {
    for ( byte = 0; byte < Nbytes; byte++ ) 
      fread ( &tmpchar, 1, 1, statefile);
    /*fscanf(statefile, "%i %i %i", &tmp_cellnum, &tmp_Nveg, 
	   &tmp_Nband);*/
    fread( &tmp_cellnum, 1, sizeof(int), statefile );
    fread( &tmp_Nveg, 1, sizeof(int), statefile );
    fread( &tmp_Nband, 1, sizeof(int), statefile );
    fread( &Nbytes, 1, sizeof(int), statefile );
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
  for ( nidx = 0; nidx < options.Nnode; nidx++ ) 
    //fscanf(statefile," %lf", &depth_node[nidx]);
    //fread( &depth_node[nidx], 1, sizeof(double), statefile );
    fread( &soil_con->dz_node[nidx], 1, sizeof(double), statefile );

  // recompute distance between soil nodes
  //if ( options.Nnode > 1 )
    //compute_dz(soil_con->dz_node, depth_node, options.Nnode, soil_con->dp);
  //else soil_con->dz_node[0] = 0;
  if ( options.Nnode == 1 ) soil_con->dz_node[0] = 0;
    
  // read storm parameters
  fread( init_STILL_STORM, 1, sizeof(char), statefile );
  fread( init_DRY_TIME, 1, sizeof(int), statefile );

  /* Input for all vegetation types */
  for ( veg = 0; veg <= Nveg; veg++ ) {
    
    // read distributed precipitation fraction
    fread( &prcp->mu[veg], 1, sizeof(double), statefile );

    /* Input for all snow bands */
    for ( band = 0; band < Nbands; band++ ) {
      
      /* Read cell identification information */
      //if ( fscanf(statefile,"%i %i", &iveg, &iband) == EOF ) 
      if ( fread( &iveg, 1, sizeof(int), statefile) != sizeof(int) ) 
	nrerror("End of model state file found unexpectedly");
      if ( fread( &iband, 1, sizeof(int), statefile) != sizeof(int) ) 
	nrerror("End of model state file found unexpectedly");
      if ( iveg != veg || iband != band ) {
	fprintf(stderr,"The vegetation and snow band indexs in the model state file (veg = %i, band = %i) do not match those currently requested (veg = %i , band = %i).  Model state file must be stored with variables for all vegetation indexed by variables for all snow bands.", iveg, iband, veg, band);
	nrerror(ErrStr);
      }
      
      for ( dist = 0; dist < Ndist; dist ++ ) {
	// Store both wet and dry fractions if using distributed precipitation

	/* Read average total soil moisture */
	for ( lidx = 0; lidx < options.Nlayer; lidx++ ) {
	  //if ( fscanf(statefile," %lf", 
	  //&cell[WET][veg][band].layer[lidx].moist) == EOF ) 
	  if ( fread( &cell[dist][veg][band].layer[lidx].moist, 1, 
		      sizeof(double), statefile ) != sizeof(double) )
	    nrerror("End of model state file found unexpectedly");
	}
      
	/* Read average ice content */
	for ( lidx = 0; lidx < options.Nlayer; lidx++ ) {
	  //if ( fscanf(statefile," %lf", 
	  //&cell[WET][veg][band].layer[lidx].ice) == EOF ) 
	  if ( fread( &cell[dist][veg][band].layer[lidx].ice, 1, 
		      sizeof(double), statefile ) != sizeof(double) ) 
	    nrerror("End of model state file found unexpectedly");
	}
      
	/* Read average dew storage */
	if ( veg < Nveg ) {
	  //if ( fscanf(statefile," %lf", &veg_var[WET][veg][band].Wdew) == EOF ) 
	  //nrerror("End of model state file found unexpectedly");
	  if ( fread( &veg_var[dist][veg][band].Wdew, 1, sizeof(double), 
		      statefile ) != sizeof(double) ) 
	    nrerror("End of model state file found unexpectedly");
	}
      }
      
      /* Read snow data */
      /*if ( fscanf(statefile," %i %lf %lf %lf %lf %lf", 
		  &snow[veg][band].last_snow, &snow[veg][band].swq, 
		  &snow[veg][band].surf_temp, &snow[veg][band].pack_temp, 
		  &snow[veg][band].density, &snow[veg][band].snow_canopy) 
	   == EOF ) 
	nrerror("End of model state file found unexpectedly");*/
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
      if(snow[veg][band].density > 0.) 
	snow[veg][band].depth = 1000. * snow[veg][band].swq 
	  / snow[veg][band].density;
      
      /* Read soil thermal node temperatures */
      for ( nidx = 0; nidx < options.Nnode; nidx++ ) 
	//if ( fscanf(statefile," %lf", &energy[veg][band].T[nidx]) == EOF ) 
  	if ( fread( &energy[veg][band].T[nidx], 1, sizeof(double), 
		    statefile ) != sizeof(double) )
	  nrerror("End of model state file found unexpectedly");
    
    }
  }
}
