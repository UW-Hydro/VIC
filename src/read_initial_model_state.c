#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vicNl.h>

static char vcid[] = "$Id$";

void read_initial_model_state(FILE                *init_state,
			      all_vars_struct     *all_vars,
			      global_param_struct *gp,
			      int                  Nveg,
			      int                  Nbands,
			      int                  cellnum,
			      soil_con_struct     *soil_con,
			      lake_con_struct      lake_con)
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
            desired record (Nveg and Nbands, respectively).			TJB
  01-Nov-04 Modified to read state files containing SPATIAL_FROST
	    and LAKE_MODEL state variables.					TJB
  02-Nov-04 Added a few more lake state variables.				TJB
  03-Nov-04 Now reads extra_veg from state file.				TJB
  2005-Apr-10 Fixed incorrect check on soil node depths.			TJB
  2005-Jan-10 modified to read lake nodal variables for each of the
	      active nodes.							JCA
  2006-Jun-16 Skip reading if areafract < 0.					GCT
  2006-Aug-23 Changed order of fread/fwrite statements from ...1, sizeof...
	      to ...sizeof, 1,...						GCT
  2006-Sep-07 Changed "Skip reading if areafract < 0" to "<=0".			GCT
  2006-Oct-16 Merged infiles and outfiles structs into filep_struct;
	      This included moving infiles.statefile to filep.init_state.	TJB
  2006-Nov-07 Removed LAKE_MODEL option.					TJB
  2007-Apr-28 modified to read Zsum_node.					JCA
  2007-May-07 Fixed fread checks to make sure correct number of items were
	      read in rather than the size of the item read in.			JCA
  2007-May-07 Nsum and sum removed from declaration.				JCA
  2007-Aug-24 Added features for EXCESS_ICE option.				JCA
  2007-Sep-14 Fixed bug for read-in during EXCESS_ICE option.			JCA
  2007-Sep-18 Check for soil moist exceeding max moist moved from
	      here to initialize_model_state.					JCA
  2007-Nov-06 New list of lake state variables.					LCB via TJB
  2009-Jul-31 Removed extra lake/wetland veg tile; updated set of lake state
	      variables.							TJB
  2009-Aug-27 Now once again expects to read data for all bands, regardless of
	      whether they have area > 0.  This makes it much easier to ensure
	      that the value of Nbands stored in the state file matches the number
	      of bands actually stored in the state file.			TJB
  2009-Sep-28 Now stores soil, snow, and energy states from lake separately
	      from wetland.							TJB
  2010-Jan-10 Corrected typo in condition for checking Wdew.			TJB
  2012-Jan-01 Removed lake area condition from logic determining whether to read
	      lake state data.  Now, if options.LAKES is TRUE, every grid cell
	      will save lake state data.  If no lake is present, default NULL
	      values will be stored.						TJB
  2013-Jul-25 Added soil carbon terms.						TJB
  2013-Dec-26 Removed EXCESS_ICE option.				TJB
  2013-Dec-27 Moved SPATIAL_FROST to options_struct.			TJB
  2013-Dec-28 Removed NO_REWIND option.					TJB
  2014-Mar-28 Removed DIST_PRCP option.					TJB
*********************************************************************/
{
  extern option_struct options;

  char   tmpstr[MAXSTRING];
  char   ErrStr[MAXSTRING];
  char   tmpchar;
  double tmpval;
  double depth_node[MAX_NODES];
  int    veg, iveg;
  int    band, iband;
  int    lidx;
  int    nidx;
  int    tmp_cellnum;
  int    tmp_Nveg;
  int    tmp_Nband;
  int    tmp_char;
  int    byte, Nbytes;
  int    tmp_int, node;
  int    frost_area;

  cell_data_struct      **cell;
  snow_data_struct      **snow;
  energy_bal_struct     **energy;
  veg_var_struct        **veg_var;
  lake_var_struct        *lake_var;
  
  cell    = all_vars->cell;
  veg_var = all_vars->veg_var;
  snow    = all_vars->snow;
  energy  = all_vars->energy;
  lake_var = &all_vars->lake_var;

  /* read cell information */
  if ( options.BINARY_STATE_FILE ) {
    fread( &tmp_cellnum, sizeof(int), 1, init_state );
    fread( &tmp_Nveg, sizeof(int), 1, init_state );
    fread( &tmp_Nband, sizeof(int), 1, init_state );
    fread( &Nbytes, sizeof(int), 1, init_state );
  }
  else 
    fscanf( init_state, "%d %d %d", &tmp_cellnum, &tmp_Nveg, &tmp_Nband );
  // Skip over unused cell information
  while ( tmp_cellnum != cellnum && !feof(init_state) ) {
    if ( options.BINARY_STATE_FILE ) {
      // skip rest of current cells info
      for ( byte = 0; byte < Nbytes; byte++ ) 
	fread ( &tmpchar, 1, 1, init_state);
      // read info for next cell
      fread( &tmp_cellnum, sizeof(int), 1, init_state );
      fread( &tmp_Nveg, sizeof(int), 1, init_state );
      fread( &tmp_Nband, sizeof(int), 1, init_state );
      fread( &Nbytes, sizeof(int), 1, init_state );
    }
    else {
      // skip rest of current cells info
      fgets(tmpstr, MAXSTRING, init_state); // skip rest of general cell info
      for ( veg = 0; veg <= tmp_Nveg; veg++ ) {
	for ( band = 0; band < tmp_Nband; band++ )
	  fgets(tmpstr, MAXSTRING, init_state); // skip snowband info
      }
      if ( options.LAKES ) {
        fgets(tmpstr, MAXSTRING, init_state); // skip lake info
      }
      // read info for next cell
      fscanf( init_state, "%d %d %d", &tmp_cellnum, &tmp_Nveg, &tmp_Nband );
    }//end if
  }//end while
  
  if ( feof(init_state) ) {
    sprintf(ErrStr, "Requested grid cell (%d) is not in the model state file.", 
	    cellnum);
    nrerror(ErrStr);
  }
  
  if ( tmp_Nveg != Nveg ) {
    sprintf(ErrStr,"The number of vegetation types in cell %d (%d) does not equal that defined in vegetation parameter file (%d).  Check your input files.", cellnum, tmp_Nveg, Nveg);
    nrerror(ErrStr);
  }
  if ( tmp_Nband != Nbands ) {
    sprintf(ErrStr,"The number of snow bands in cell %d (%d) does not equal that defined in the snow band file (%d).  Check your input files.", cellnum, tmp_Nband, Nbands);
    nrerror(ErrStr);
  }
 
  /* Read soil thermal node deltas */
  for ( nidx = 0; nidx < options.Nnode; nidx++ ) {
    if ( options.BINARY_STATE_FILE ) 
      fread( &soil_con->dz_node[nidx], sizeof(double), 1, init_state );
    else 
      fscanf( init_state, "%lf", &soil_con->dz_node[nidx] );
  }
  if ( options.Nnode == 1 ) soil_con->dz_node[0] = 0;
  
  /* Read soil thermal node depths */
  for ( nidx = 0; nidx < options.Nnode; nidx++ ) {
    if ( options.BINARY_STATE_FILE ) 
      fread( &soil_con->Zsum_node[nidx], sizeof(double), 1, init_state );
    else 
      fscanf( init_state, "%lf", &soil_con->Zsum_node[nidx] );
  }
  if ( options.Nnode == 1 ) soil_con->Zsum_node[0] = 0;
  if ( soil_con->Zsum_node[options.Nnode-1] - soil_con->dp > SMALL) {
    fprintf( stderr, "WARNING: Sum of soil nodes (%f) exceeds defined damping depth (%f).  Resetting damping depth.\n", soil_con->Zsum_node[options.Nnode-1], soil_con->dp );
    soil_con->dp = soil_con->Zsum_node[options.Nnode-1];
  }

  /* Input for all vegetation types */
  for ( veg = 0; veg <= Nveg; veg++ ) {
    
    /* Input for all snow bands */
    for ( band = 0; band < Nbands; band++ ) {
      /* Read cell identification information */
      if ( options.BINARY_STATE_FILE ) {
	if ( fread( &iveg, sizeof(int), 1, init_state) != 1 ) 
	  nrerror("End of model state file found unexpectedly");
	if ( fread( &iband, sizeof(int), 1, init_state) != 1 ) 
	  nrerror("End of model state file found unexpectedly");
      }
      else {
	if ( fscanf(init_state,"%d %d", &iveg, &iband) == EOF ) 
	  nrerror("End of model state file found unexpectedly");
      }
      if ( iveg != veg || iband != band ) {
	fprintf(stderr,"The vegetation and snow band indices in the model state file (veg = %d, band = %d) do not match those currently requested (veg = %d , band = %d).  Model state file must be stored with variables for all vegetation indexed by variables for all snow bands.\n", iveg, iband, veg, band);
	nrerror(ErrStr);
      }
      
      /* Read total soil moisture */
      for ( lidx = 0; lidx < options.Nlayer; lidx++ ) {
	if ( options.BINARY_STATE_FILE ) {
	  if ( fread( &cell[veg][band].layer[lidx].moist,
			sizeof(double), 1, init_state ) != 1 )
	    nrerror("End of model state file found unexpectedly");
	}
	else {
	  if ( fscanf(init_state," %lf", 
			&cell[veg][band].layer[lidx].moist) == EOF ) 
	    nrerror("End of model state file found unexpectedly");
	}
      }
	
      /* Read average ice content */
      for ( lidx = 0; lidx < options.Nlayer; lidx++ ) {
	for ( frost_area = 0; frost_area < options.Nfrost; frost_area++ ) {
	  if ( options.BINARY_STATE_FILE ) {
	    if ( fread( &cell[veg][band].layer[lidx].ice[frost_area],
			  sizeof(double), 1, init_state ) != 1 )
	      nrerror("End of model state file found unexpectedly");
	  }
	  else {
	    if ( fscanf(init_state," %lf", 
			  &cell[veg][band].layer[lidx].ice[frost_area]) == EOF ) 
	      nrerror("End of model state file found unexpectedly");
	  }
	}
      }
	
      if ( veg < Nveg ) {
	/* Read dew storage */
	if ( options.BINARY_STATE_FILE ) {
	  if ( fread( &veg_var[veg][band].Wdew, sizeof(double), 1, 
			init_state ) != 1 ) 
	    nrerror("End of model state file found unexpectedly");
	}
	else {
	  if ( fscanf(init_state," %lf", &veg_var[veg][band].Wdew) == EOF ) 
	    nrerror("End of model state file found unexpectedly");
	}

        if ( options.CARBON ) {
          if ( options.BINARY_STATE_FILE ) {
	    /* Read cumulative annual NPP */
	    if ( fread( &(veg_var[veg][band].AnnualNPP), sizeof(double), 1, init_state ) != 1 )
	      nrerror("End of model state file found unexpectedly");
	    if ( fread( &(veg_var[veg][band].AnnualNPPPrev), sizeof(double), 1, init_state ) != 1 )
	      nrerror("End of model state file found unexpectedly");
	    /* Read Soil Carbon Storage */
	    if ( fread( &(cell[veg][band].CLitter), sizeof(double), 1, init_state ) != 1 )
	      nrerror("End of model state file found unexpectedly");
	    if ( fread( &(cell[veg][band].CInter), sizeof(double), 1, init_state ) != 1 )
	      nrerror("End of model state file found unexpectedly");
	    if ( fread( &(cell[veg][band].CSlow), sizeof(double), 1, init_state ) != 1 )
	      nrerror("End of model state file found unexpectedly");
          }
          else {
	    /* Read cumulative annual NPP */
	    if ( fscanf(init_state," %lf", &veg_var[veg][band].AnnualNPP) == EOF ) 
	      nrerror("End of model state file found unexpectedly");
	    if ( fscanf(init_state," %lf", &veg_var[veg][band].AnnualNPPPrev) == EOF ) 
	      nrerror("End of model state file found unexpectedly");
	    /* Read Soil Carbon Storage */
	    if ( fscanf(init_state," %lf", &cell[veg][band].CLitter) == EOF ) 
	      nrerror("End of model state file found unexpectedly");
	    if ( fscanf(init_state," %lf", &cell[veg][band].CInter) == EOF ) 
	      nrerror("End of model state file found unexpectedly");
	    if ( fscanf(init_state," %lf", &cell[veg][band].CSlow) == EOF ) 
	      nrerror("End of model state file found unexpectedly");
          }
        }

      }
      
      /* Read snow data */
      if ( options.BINARY_STATE_FILE ) {
	if ( fread( &snow[veg][band].last_snow, sizeof(int), 1, 
		    init_state ) != 1 )
	  nrerror("End of model state file found unexpectedly");
	if ( fread( &snow[veg][band].MELTING, sizeof(char), 1, 
		    init_state ) != 1 )
	  nrerror("End of model state file found unexpectedly");
	if ( fread( &snow[veg][band].coverage, sizeof(double), 1, 
		    init_state ) != 1 )
	  nrerror("End of model state file found unexpectedly");
	if ( fread( &snow[veg][band].swq, sizeof(double), 1, 
		    init_state ) != 1 )
	  nrerror("End of model state file found unexpectedly");
	if ( fread( &snow[veg][band].surf_temp, sizeof(double), 1, 
		    init_state ) != 1 )
	  nrerror("End of model state file found unexpectedly");
	if ( fread( &snow[veg][band].surf_water, sizeof(double), 1, 
		    init_state ) != 1 )
	  nrerror("End of model state file found unexpectedly");
	if ( fread( &snow[veg][band].pack_temp, sizeof(double), 1, 
		    init_state ) != 1 ) 
	  nrerror("End of model state file found unexpectedly");
	if ( fread( &snow[veg][band].pack_water, sizeof(double), 1, 
		    init_state ) != 1 )
	  nrerror("End of model state file found unexpectedly");
	if ( fread( &snow[veg][band].density, sizeof(double), 1, 
		    init_state ) != 1 )
	  nrerror("End of model state file found unexpectedly");
	if ( fread( &snow[veg][band].coldcontent, sizeof(double), 1, 
		    init_state ) != 1 )
	  nrerror("End of model state file found unexpectedly");
	if ( fread( &snow[veg][band].snow_canopy, sizeof(double), 1, 
		    init_state ) != 1 )
	  nrerror("End of model state file found unexpectedly");
      }
      else {
	if ( fscanf(init_state," %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf", 
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
	  if ( fread( &energy[veg][band].T[nidx], sizeof(double), 1, 
		      init_state ) != 1 )
	    nrerror("End of model state file found unexpectedly");
	}
	else {
	  if ( fscanf(init_state," %lf", &energy[veg][band].T[nidx]) == EOF )
	    nrerror("End of model state file found unexpectedly");
	}
      }
    }
  }
  if ( options.LAKES ) {
    if ( options.BINARY_STATE_FILE ) {

      /* Read total soil moisture */
      for ( lidx = 0; lidx < options.Nlayer; lidx++ ) {
	if ( fread( &lake_var->soil.layer[lidx].moist, sizeof(double), 1, init_state ) != 1 )
	  nrerror("End of model state file found unexpectedly");
      }
	
      /* Read average ice content */
      for ( lidx = 0; lidx < options.Nlayer; lidx++ ) {
        for ( frost_area = 0; frost_area < options.Nfrost; frost_area++ ) {
          if ( fread( &lake_var->soil.layer[lidx].ice[frost_area], sizeof(double), 1, init_state ) != 1 )
            nrerror("End of model state file found unexpectedly");
	}
      }

      if ( options.CARBON ) {
        /* Read Soil Carbon Storage */
        if ( fread( &(lake_var->soil.CLitter), sizeof(double), 1, init_state ) != 1 )
	  nrerror("End of model state file found unexpectedly");
	if ( fread( &(lake_var->soil.CInter), sizeof(double), 1, init_state ) != 1 )
	  nrerror("End of model state file found unexpectedly");
	if ( fread( &(lake_var->soil.CSlow), sizeof(double), 1, init_state ) != 1 )
	  nrerror("End of model state file found unexpectedly");
      }

      /* Read snow data */
      if ( fread( &lake_var->snow.last_snow, sizeof(int), 1, init_state ) != 1 )
	nrerror("End of model state file found unexpectedly");
      if ( fread( &lake_var->snow.MELTING, sizeof(char), 1, init_state ) != 1 )
	nrerror("End of model state file found unexpectedly");
      if ( fread( &lake_var->snow.coverage, sizeof(double), 1, init_state ) != 1 )
	nrerror("End of model state file found unexpectedly");
      if ( fread( &lake_var->snow.swq, sizeof(double), 1, init_state ) != 1 )
	nrerror("End of model state file found unexpectedly");
      if ( fread( &lake_var->snow.surf_temp, sizeof(double), 1, init_state ) != 1 )
	nrerror("End of model state file found unexpectedly");
      if ( fread( &lake_var->snow.surf_water, sizeof(double), 1, init_state ) != 1 )
	nrerror("End of model state file found unexpectedly");
      if ( fread( &lake_var->snow.pack_temp, sizeof(double), 1, init_state ) != 1 ) 
	nrerror("End of model state file found unexpectedly");
      if ( fread( &lake_var->snow.pack_water, sizeof(double), 1, init_state ) != 1 )
	nrerror("End of model state file found unexpectedly");
      if ( fread( &lake_var->snow.density, sizeof(double), 1, init_state ) != 1 )
	nrerror("End of model state file found unexpectedly");
      if ( fread( &lake_var->snow.coldcontent, sizeof(double), 1, init_state ) != 1 )
	nrerror("End of model state file found unexpectedly");
      if ( fread( &lake_var->snow.snow_canopy, sizeof(double), 1, init_state ) != 1 )
	nrerror("End of model state file found unexpectedly");
      if(lake_var->snow.density > 0.) 
	lake_var->snow.depth = 1000. * lake_var->snow.swq / lake_var->snow.density;
      
      /* Read soil thermal node temperatures */
      for ( nidx = 0; nidx < options.Nnode; nidx++ ) {
	if ( fread( &lake_var->energy.T[nidx], sizeof(double), 1, init_state ) != 1 )
	  nrerror("End of model state file found unexpectedly");
      }

      /* Read lake-specific variables */
      if ( fread( &lake_var->activenod, sizeof(int), 1, init_state ) != 1 )
	nrerror("End of model state file found unexpectedly");
      if ( fread( &lake_var->dz, sizeof(double), 1, init_state ) != 1 )
	nrerror("End of model state file found unexpectedly");
      if ( fread( &lake_var->surfdz, sizeof(double), 1, init_state ) != 1 )
	nrerror("End of model state file found unexpectedly");
      if ( fread( &lake_var->ldepth, sizeof(double), 1, init_state ) != 1 )
	nrerror("End of model state file found unexpectedly");
      for ( node = 0; node <= lake_var->activenod; node++ ) {
        if ( fread( &lake_var->surface[node], sizeof(double), 1, init_state ) != 1 )
	  nrerror("End of model state file found unexpectedly");
      }
      if ( fread( &lake_var->sarea, sizeof(double), 1, init_state ) != 1 )
	nrerror("End of model state file found unexpectedly");
      if ( fread( &lake_var->volume, sizeof(double), 1, init_state ) != 1 )
	nrerror("End of model state file found unexpectedly");
      for ( node = 0; node < lake_var->activenod; node++ ) {
        if ( fread( &lake_var->temp[node], sizeof(double), 1, init_state ) != 1 )
	  nrerror("End of model state file found unexpectedly");
      }
      if ( fread( &lake_var->tempavg, sizeof(double), 1, init_state ) != 1 )
	nrerror("End of model state file found unexpectedly");
      if ( fread( &lake_var->areai, sizeof(double), 1, init_state ) != 1 )
	nrerror("End of model state file found unexpectedly");
      if ( fread( &lake_var->new_ice_area, sizeof(double), 1, init_state ) != 1 )
	nrerror("End of model state file found unexpectedly");
      if ( fread( &lake_var->ice_water_eq, sizeof(double), 1, init_state ) != 1 )
	nrerror("End of model state file found unexpectedly");
      if ( fread( &lake_var->hice, sizeof(double), 1, init_state ) != 1 )
	nrerror("End of model state file found unexpectedly");
      if ( fread( &lake_var->tempi, sizeof(double), 1, init_state ) != 1 )
	nrerror("End of model state file found unexpectedly");
      if ( fread( &lake_var->swe, sizeof(double), 1, init_state ) != 1 )
	nrerror("End of model state file found unexpectedly");
      if ( fread( &lake_var->surf_temp, sizeof(double), 1, init_state ) != 1 )
	nrerror("End of model state file found unexpectedly");
      if ( fread( &lake_var->pack_temp, sizeof(double), 1, init_state ) != 1 )
	nrerror("End of model state file found unexpectedly");
      if ( fread( &lake_var->coldcontent, sizeof(double), 1, init_state ) != 1 )
	nrerror("End of model state file found unexpectedly");
      if ( fread( &lake_var->surf_water, sizeof(double), 1, init_state ) != 1 )
	nrerror("End of model state file found unexpectedly");
      if ( fread( &lake_var->pack_water, sizeof(double), 1, init_state ) != 1 )
	nrerror("End of model state file found unexpectedly");
      if ( fread( &lake_var->SAlbedo, sizeof(double), 1, init_state ) != 1 )
	nrerror("End of model state file found unexpectedly");
      if ( fread( &lake_var->sdepth, sizeof(double), 1, init_state ) != 1 )
	nrerror("End of model state file found unexpectedly");
      
    }
    else {

      /* Read total soil moisture */
      for ( lidx = 0; lidx < options.Nlayer; lidx++ ) {
	if ( fscanf(init_state," %lf", &lake_var->soil.layer[lidx].moist) == EOF ) 
          nrerror("End of model state file found unexpectedly");
      }
	
      /* Read average ice content */
      for ( lidx = 0; lidx < options.Nlayer; lidx++ ) {
        for ( frost_area = 0; frost_area < options.Nfrost; frost_area++ ) {
	  if ( fscanf(init_state," %lf", &lake_var->soil.layer[lidx].ice[frost_area]) == EOF ) 
	    nrerror("End of model state file found unexpectedly");
        }
      }
	
      if (options.CARBON) {
	/* Read Soil Carbon Storage */
	if ( fscanf(init_state," %lf", &lake_var->soil.CLitter) == EOF ) 
	  nrerror("End of model state file found unexpectedly");
	if ( fscanf(init_state," %lf", &lake_var->soil.CInter) == EOF ) 
	  nrerror("End of model state file found unexpectedly");
	if ( fscanf(init_state," %lf", &lake_var->soil.CSlow) == EOF ) 
	  nrerror("End of model state file found unexpectedly");
      }

      /* Read snow data */
      if ( fscanf(init_state," %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf", 
		  &lake_var->snow.last_snow, &tmp_char,
		  &lake_var->snow.coverage, &lake_var->snow.swq, 
		  &lake_var->snow.surf_temp, &lake_var->snow.surf_water, 
		  &lake_var->snow.pack_temp, &lake_var->snow.pack_water, 
		  &lake_var->snow.density, &lake_var->snow.coldcontent, 
		  &lake_var->snow.snow_canopy)
	   == EOF ) 
        nrerror("End of model state file found unexpectedly");
      lake_var->snow.MELTING = (char)tmp_char;
      if(lake_var->snow.density > 0.) 
	lake_var->snow.depth = 1000. * lake_var->snow.swq / lake_var->snow.density;
      
      /* Read soil thermal node temperatures */
      for ( nidx = 0; nidx < options.Nnode; nidx++ ) {
	if ( fscanf(init_state," %lf", &lake_var->energy.T[nidx]) == EOF )
	  nrerror("End of model state file found unexpectedly");
      }

      /* Read lake-specific variables */
      if ( fscanf(init_state," %d", &lake_var->activenod) == EOF )
	nrerror("End of model state file found unexpectedly");
      if ( fscanf(init_state," %lf", &lake_var->dz) == EOF )
	nrerror("End of model state file found unexpectedly");
      if ( fscanf(init_state," %lf", &lake_var->surfdz) == EOF )
	nrerror("End of model state file found unexpectedly");
      if ( fscanf(init_state," %lf", &lake_var->ldepth) == EOF )
	nrerror("End of model state file found unexpectedly");
      for ( node = 0; node <= lake_var->activenod; node++ ) {
        if ( fscanf(init_state," %lf", &lake_var->surface[node]) == EOF )
	  nrerror("End of model state file found unexpectedly");
      }
      if ( fscanf(init_state," %lf", &lake_var->sarea) == EOF )
	nrerror("End of model state file found unexpectedly");
      if ( fscanf(init_state," %lf", &lake_var->volume) == EOF )
	nrerror("End of model state file found unexpectedly");
      for ( node = 0; node < lake_var->activenod; node++ ) {
        if ( fscanf(init_state," %lf", &lake_var->temp[node]) == EOF )
	  nrerror("End of model state file found unexpectedly");
      }
      if ( fscanf(init_state," %lf", &lake_var->tempavg) == EOF )
	nrerror("End of model state file found unexpectedly");
      if ( fscanf(init_state," %lf", &lake_var->areai) == EOF )
	nrerror("End of model state file found unexpectedly");
      if ( fscanf(init_state," %lf", &lake_var->new_ice_area) == EOF )
	nrerror("End of model state file found unexpectedly");
      if ( fscanf(init_state," %lf", &lake_var->ice_water_eq) == EOF )
	nrerror("End of model state file found unexpectedly");
      if ( fscanf(init_state," %lf", &lake_var->hice) == EOF )
	nrerror("End of model state file found unexpectedly");
      if ( fscanf(init_state," %lf", &lake_var->tempi) == EOF )
	nrerror("End of model state file found unexpectedly");
      if ( fscanf(init_state," %lf", &lake_var->swe) == EOF )
	nrerror("End of model state file found unexpectedly");
      if ( fscanf(init_state," %lf", &lake_var->surf_temp) == EOF )
	nrerror("End of model state file found unexpectedly");
      if ( fscanf(init_state," %lf", &lake_var->pack_temp) == EOF )
	nrerror("End of model state file found unexpectedly");
      if ( fscanf(init_state," %lf", &lake_var->coldcontent) == EOF )
	nrerror("End of model state file found unexpectedly");
      if ( fscanf(init_state," %lf", &lake_var->surf_water) == EOF )
	nrerror("End of model state file found unexpectedly");
      if ( fscanf(init_state," %lf", &lake_var->pack_water) == EOF )
	nrerror("End of model state file found unexpectedly");
      if ( fscanf(init_state," %lf", &lake_var->SAlbedo) == EOF )
	nrerror("End of model state file found unexpectedly");
      if ( fscanf(init_state," %lf", &lake_var->sdepth) == EOF )
	nrerror("End of model state file found unexpectedly");
      
    }
  }

}
