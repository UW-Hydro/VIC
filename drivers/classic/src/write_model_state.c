#include <vic_def.h>
#include <vic_run.h>
#include <vic_driver_classic.h>
#define MAX_PFT 20

void write_model_state(all_vars_struct    *all_vars,
		       global_param_struct *gp,
		       int                  Nveg,
		       int                  cellnum,
		       filep_struct        *filep,
		       soil_con_struct     *soil_con,
		       lake_con_struct      lake_con,
		       cn_data_struct      *cn)
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
  09-Oct-03 Removed initial space on veg/band info line in ASCII file.		TJB
  26-Oct-04 Changed calculation of Nbytes in binary state file to
	    account for bare soil values (extra veg class per grid
	    cell).  Without this fix, attempts to skip grid cells
	    fail.								TJB
  01-Nov-04 Added storage of state variables for SPATIAL_FROST and
	    LAKE_MODEL.								TJB
  02-Nov-04 Added a few more lake state variables.				TJB
  03-Nov-04 Now outputs extra_veg to aid other programs in parsing
	    state files.							TJB
  2005-Dec-07 STATE_FILE option is set in global file.				GCT
  2005-Jan-10 writes temp[0] instead of tp_in for lake skin surface
	      temperature.							JCA
  2005-Jan-10 modified to write lake nodal variables for each of the
	      active nodes.							JCA
  2006-Jun-16 Skip writing snow band if areafract < 0.				GCT
  2006-Aug-23 Changed order of fread/fwrite statements from ...1, sizeof...
	      to ...sizeof, 1,...						GCT
  2006-Sep-07 Changed "Skip writing snow band if areafract < 0" to "<=0".	GCT
  2006-Oct-16 Merged infiles and outfiles structs into filep_struct.		TJB
  2006-Nov-07 Removed LAKE_MODEL option.					TJB
  2007-Apr-24 Modified to write Zsum_node.					JCA
  2007-Apr-25 Removed variable Nsum.						JCA
  2007-Aug-24 Added features for EXCESS_ICE option.				JCA
  2007-Nov-06 New list of lake state variables.					LCB via TJB
  2009-Jul-31 Removed extra lake/wetland veg tile; updated set of lake state
	      variables.							TJB
  2009-Aug-27 Now once again writes data for all bands, regardless of
	      whether they have area > 0.  This makes it much easier to ensure
	      that the value of Nbands stored in the state file matches the number
	      of bands actually stored in the state file.			TJB
  2009-Sep-28 Now stores soil, snow, and energy states from lake separately
	      from wetland.							TJB
  2010-Mar-05 Fixed typo in writing of number of lake active nodes.		TJB
  2012-Jan-01 Removed lake area condition from logic determining whether to write
	      lake state data.  Now, if options.LAKES is TRUE, every grid cell
	      will save lake state data.  If no lake is present, default NULL
	      values will be stored.						TJB
  2013-Jul-25 Added soil carbon terms.						TJB
  2013-Dec-26 Removed EXCESS_ICE option.				TJB
  2013-Dec-27 Moved SPATIAL_FROST to options_struct.			TJB
  2014-Mar-28 Removed DIST_PRCP option.					TJB
  2014-Nov-25 Added CN quantities.                                      MAB
*********************************************************************/
{
  extern option_struct options;

  double tmpval;
  int    veg;
  int    band;
  int    lidx;
  int    nidx;
  int    Nbands;
  int    byte, Nbytes;
  int    frost_area;

  cell_data_struct      **cell;
  snow_data_struct      **snow;
  energy_bal_struct     **energy;
  veg_var_struct        **veg_var;
  lake_var_struct         lake_var;
  int    node;

  Nbands = options.SNOW_BAND;

  cell    = all_vars->cell;
  veg_var = all_vars->veg_var;
  snow    = all_vars->snow;
  energy  = all_vars->energy;
  lake_var = all_vars->lake_var;
 
  /* write cell information */
  if ( options.BINARY_STATE_FILE ) {
    fwrite( &cellnum, sizeof(int), 1, filep->statefile );
    fwrite( &Nveg, sizeof(int), 1, filep->statefile );
    fwrite( &Nbands, sizeof(int), 1, filep->statefile );
  }
  else {
    fprintf( filep->statefile, "%i %i %i", cellnum, Nveg, Nbands );
  }
  // This stores the number of bytes from after this value to the end 
  // of the line.  DO NOT CHANGE unless you have changed the values
  // written to the state file.
  // IF YOU EDIT THIS FILE: UPDATE THIS VALUE!
  if ( options.BINARY_STATE_FILE ) {
    Nbytes =   options.Nnode * sizeof(double) // dz_node
	       + options.Nnode * sizeof(double) // Zsum_node
	       + (Nveg+1) * Nbands * 2 * sizeof(int) // veg & band
	       + (Nveg+1) * Nbands * options.Nlayer * sizeof(double) // soil moisture
	       + (Nveg+1) * Nbands * options.Nlayer * options.Nfrost * sizeof(double) // soil ice
	       + Nveg * Nbands * sizeof(double); // dew
    if ( options.CARBON ) {
      /* Carbon-specific state vars */
      Nbytes += Nveg * Nbands * 5 * sizeof(double); // AnnualNPP, AnnualNPPPrev, and 3 soil carbon storages
    }
    Nbytes += (Nveg+1) * Nbands * sizeof(int) // last_snow
	       + (Nveg+1) * Nbands * sizeof(char) // MELTING
	       + (Nveg+1) * Nbands * sizeof(double) * 9 // other snow parameters
	       + (Nveg+1) * Nbands * options.Nnode * sizeof(double); // soil temperatures
    if ( options.LAKES ) {
      /* Lake/wetland tiles have lake-specific state vars */
      Nbytes += sizeof(int) // activenod
	+ sizeof(double) // dz
	+ sizeof(double) // surfdz
	+ sizeof(double) // ldepth
	+ (lake_var.activenod+1) * sizeof(double) // surface
	+ sizeof(double) // sarea
	+ sizeof(double) // volume
	+ lake_var.activenod * sizeof(double) // temp
	+ sizeof(double) // tempavg
	+ sizeof(double) // areai
	+ sizeof(double) // new_ice_area
	+ sizeof(double) // ice_water_eq
	+ sizeof(double) // hice
	+ sizeof(double) // tempi
	+ sizeof(double) // swe
	+ sizeof(double) // surf_temp
	+ sizeof(double) // pack_temp
	+ sizeof(double) // coldcontent
	+ sizeof(double) // surf_water
	+ sizeof(double) // pack_water
	+ sizeof(double) // SAlbedo
	+ sizeof(double) // sdepth
	+ options.Nlayer * sizeof(double) // soil moisture
	+ options.Nlayer * options.Nfrost * sizeof(double) // soil ice
	+ sizeof(int) // last_snow
	+ sizeof(char) // MELTING
	+ sizeof(double) * 9 // other snow parameters
	+ options.Nnode * sizeof(double) // soil temperatures
	;
      if ( options.CARBON ) {
        /* Carbon-specific state vars */
        Nbytes += 3 * sizeof(double); // 3 soil carbon storages
      }
      if ( options.CARBON == CN_NORMAL || options.CARBON == CN_ADECOMP ) {
      /* CN-specific state vars */
	Nbytes += Nbands * (MAX_PFT + 1) * sizeof(double) // dormant_flag
	  + Nbands * (MAX_PFT + 1) * sizeof(double) // days_active
	  + Nbands * (MAX_PFT + 1) * sizeof(double) // onset_flag
	  + Nbands * (MAX_PFT + 1) * sizeof(double) // onset_counter
	  + Nbands * (MAX_PFT + 1) * sizeof(double) // onset_gddflag
	  + Nbands * (MAX_PFT + 1) * sizeof(double) // onset_fdd
	  + Nbands * (MAX_PFT + 1) * sizeof(double) // onset_gdd
	  + Nbands * (MAX_PFT + 1) * sizeof(double) // onset_swi
	  + Nbands * (MAX_PFT + 1) * sizeof(double) // offset_flag
	  + Nbands * (MAX_PFT + 1) * sizeof(double) // offset_counter
	  + Nbands * (MAX_PFT + 1) * sizeof(double) // offset_fdd
	  + Nbands * (MAX_PFT + 1) * sizeof(double) // offset_swi
	  + Nbands * (MAX_PFT + 1) * sizeof(double) // lgsf
	  + Nbands * (MAX_PFT + 1) * sizeof(double) // bglfr
	  + Nbands * (MAX_PFT + 1) * sizeof(double) // bgtr
	  + Nbands * (MAX_PFT + 1) * sizeof(double) // dayl
	  + Nbands * (MAX_PFT + 1) * sizeof(double) // prev_dayl
	  + Nbands * (MAX_PFT + 1) * sizeof(double) // annavg_t2m
	  + Nbands * (MAX_PFT + 1) * sizeof(double) // tempavg_t2m
	  + Nbands * (MAX_PFT + 1) * sizeof(double) // gpp
	  + Nbands * (MAX_PFT + 1) * sizeof(double) // availc
	  + Nbands * (MAX_PFT + 1) * sizeof(double) // xsmrpool_recover
	  + Nbands * (MAX_PFT + 1) * sizeof(double) // alloc_pnow
	  + Nbands * (MAX_PFT + 1) * sizeof(double) // c_allometry
	  + Nbands * (MAX_PFT + 1) * sizeof(double) // n_allometry
	  + Nbands * (MAX_PFT + 1) * sizeof(double) // tempsum_potential_gpp
	  + Nbands * (MAX_PFT + 1) * sizeof(double) // annsum_potential_gpp
	  + Nbands * (MAX_PFT + 1) * sizeof(double) // tempmax_retransn
	  + Nbands * (MAX_PFT + 1) * sizeof(double) // annmax_retransn
	  + Nbands * (MAX_PFT + 1) * sizeof(double) // avail_retransn
	  + Nbands * (MAX_PFT + 1) * sizeof(double) // plant_nalloc
	  + Nbands * (MAX_PFT + 1) * sizeof(double) // plant_calloc
	  + Nbands * (MAX_PFT + 1) * sizeof(double) // excess_cflux
	  + Nbands * (MAX_PFT + 1) * sizeof(double) // downreg
	  + Nbands * (MAX_PFT + 1) * sizeof(double) // prev_leafc_to_litter
	  + Nbands * (MAX_PFT + 1) * sizeof(double) // prev_frootc_to_litter
	  + Nbands * (MAX_PFT + 1) * sizeof(double) // tempsum_npp
	  + Nbands * (MAX_PFT + 1) * sizeof(double) // annsum_npp
	  + Nbands * (MAX_PFT + 1) * sizeof(double) // leafc
	  + Nbands * (MAX_PFT + 1) * sizeof(double) // leafc_storage
	  + Nbands * (MAX_PFT + 1) * sizeof(double) // leafc_xfer
	  + Nbands * (MAX_PFT + 1) * sizeof(double) // frootc
	  + Nbands * (MAX_PFT + 1) * sizeof(double) // frootc_storage
	  + Nbands * (MAX_PFT + 1) * sizeof(double) // frootc_xfer
	  + Nbands * (MAX_PFT + 1) * sizeof(double) // livestemc
	  + Nbands * (MAX_PFT + 1) * sizeof(double) // livestemc_storage
	  + Nbands * (MAX_PFT + 1) * sizeof(double) // livestemc_xfer
	  + Nbands * (MAX_PFT + 1) * sizeof(double) // deadstemc
	  + Nbands * (MAX_PFT + 1) * sizeof(double) // deadstemc_storage
	  + Nbands * (MAX_PFT + 1) * sizeof(double) // deadstemc_xfer
	  + Nbands * (MAX_PFT + 1) * sizeof(double) // livecrootc
	  + Nbands * (MAX_PFT + 1) * sizeof(double) // livecrootc_storage
	  + Nbands * (MAX_PFT + 1) * sizeof(double) // livecrootc_xfer
	  + Nbands * (MAX_PFT + 1) * sizeof(double) // deadcrootc
	  + Nbands * (MAX_PFT + 1) * sizeof(double) // deadcrootc_storage
	  + Nbands * (MAX_PFT + 1) * sizeof(double) // deadcrootc_xfer
	  + Nbands * (MAX_PFT + 1) * sizeof(double) // gresp_storage
	  + Nbands * (MAX_PFT + 1) * sizeof(double) // gresp_xfer
	  + Nbands * (MAX_PFT + 1) * sizeof(double) // cpool
	  + Nbands * (MAX_PFT + 1) * sizeof(double) // xsmrpool
	  + Nbands * (MAX_PFT + 1) * sizeof(double) // pft_ctrunc
	  + Nbands * (MAX_PFT + 1) * sizeof(double) // totvegc
	  + Nbands * (MAX_PFT + 1) * sizeof(double) // leafn
	  + Nbands * (MAX_PFT + 1) * sizeof(double) // leafn_storage
	  + Nbands * (MAX_PFT + 1) * sizeof(double) // leafn_xfer
	  + Nbands * (MAX_PFT + 1) * sizeof(double) // frootn
	  + Nbands * (MAX_PFT + 1) * sizeof(double) // frootn_storage
	  + Nbands * (MAX_PFT + 1) * sizeof(double) // frootn_xfer
	  + Nbands * (MAX_PFT + 1) * sizeof(double) // livestemn
	  + Nbands * (MAX_PFT + 1) * sizeof(double) // livestemn_storage
	  + Nbands * (MAX_PFT + 1) * sizeof(double) // livestemn_xfer
	  + Nbands * (MAX_PFT + 1) * sizeof(double) // deadstemn
	  + Nbands * (MAX_PFT + 1) * sizeof(double) // deadstemn_storage
	  + Nbands * (MAX_PFT + 1) * sizeof(double) // deadstemn_xfer
	  + Nbands * (MAX_PFT + 1) * sizeof(double) // livecrootn
	  + Nbands * (MAX_PFT + 1) * sizeof(double) // livecrootn_storage
	  + Nbands * (MAX_PFT + 1) * sizeof(double) // livecrootn_xfer
	  + Nbands * (MAX_PFT + 1) * sizeof(double) // deadcrootn
	  + Nbands * (MAX_PFT + 1) * sizeof(double) // deadcrootn_storage
	  + Nbands * (MAX_PFT + 1) * sizeof(double) // deadcrootn_xfer
	  + Nbands * (MAX_PFT + 1) * sizeof(double) // retransn
	  + Nbands * (MAX_PFT + 1) * sizeof(double) // npool
	  + Nbands * (MAX_PFT + 1) * sizeof(double) // pft_ntrunc
	  + Nbands * sizeof(double) // decl
	  + Nbands * sizeof(double) // fpi
	  + Nbands * sizeof(double) // fpg
	  + Nbands * sizeof(double) // annsum_counter
	  + Nbands * sizeof(double) // cannsum_npp
	  + Nbands * sizeof(double) // cannavg_t2m
	  + Nbands * options.Nnode * sizeof(double) // watfc
	  + Nbands * sizeof(double) // me
	  + Nbands * sizeof(double) // fire_prob
	  + Nbands * sizeof(double) // mean_fire_prob
	  + Nbands * sizeof(double) // fireseasonl
	  + Nbands * sizeof(double) // ann_farea_burned
	  + Nbands * sizeof(double) // cwdc
	  + Nbands * sizeof(double) // litr1c
	  + Nbands * sizeof(double) // litr2c
	  + Nbands * sizeof(double) // litr3c
	  + Nbands * sizeof(double) // soil1c
	  + Nbands * sizeof(double) // soil2c
	  + Nbands * sizeof(double) // soil3c
	  + Nbands * sizeof(double) // soil4c
	  + Nbands * sizeof(double) // seedc
	  + Nbands * sizeof(double) // col_ctrunc
	  + Nbands * sizeof(double) // totlitc
	  + Nbands * sizeof(double) // totcolc
	  + Nbands * sizeof(double) // prod10c
	  + Nbands * sizeof(double) // prod100c
	  + Nbands * sizeof(double) // cwdn
	  + Nbands * sizeof(double) // litr1n
	  + Nbands * sizeof(double) // litr2n
	  + Nbands * sizeof(double) // litr3n
	  + Nbands * sizeof(double) // soil1n
	  + Nbands * sizeof(double) // soil2n
	  + Nbands * sizeof(double) // soil3n
	  + Nbands * sizeof(double) // soil4n
	  + Nbands * sizeof(double) // sminn
	  + Nbands * sizeof(double) // col_ntrunc
	  + Nbands * sizeof(double) // seedn
	  + Nbands * sizeof(double) // totcoln
	  + Nbands * sizeof(double) // prod10n
	  + Nbands * sizeof(double) // prod100n
	  ;
      }
    }
    fwrite( &Nbytes, sizeof(int), 1, filep->statefile );
  }
  
  /* Write soil thermal node deltas */
  for ( nidx = 0; nidx < options.Nnode; nidx++ ) {
    if ( options.BINARY_STATE_FILE )
      fwrite( &soil_con->dz_node[nidx], sizeof(double), 1,
	      filep->statefile );
    else
      fprintf( filep->statefile, " %f ", soil_con->dz_node[nidx] );
  } 
  /* Write soil thermal node depths */
  for ( nidx = 0; nidx < options.Nnode; nidx++ ) {
    if ( options.BINARY_STATE_FILE )
      fwrite( &soil_con->Zsum_node[nidx], sizeof(double), 1, 
	      filep->statefile );
    else
      fprintf( filep->statefile, " %f ", soil_con->Zsum_node[nidx] );
  }    
  if ( !options.BINARY_STATE_FILE )
    fprintf( filep->statefile, "\n" );
 
  /* Output for all vegetation types */
  for ( veg = 0; veg <= Nveg; veg++ ) {
    
    /* Output for all snow bands */
    for ( band = 0; band < Nbands; band++ ) {
      /* Write cell identification information */
      if ( options.BINARY_STATE_FILE ) {
	fwrite( &veg, sizeof(int), 1, filep->statefile );
	fwrite( &band, sizeof(int), 1, filep->statefile );
      }
      else {
	fprintf( filep->statefile, "%i %i", veg, band );
      }

      /* Write total soil moisture */
      for ( lidx = 0; lidx < options.Nlayer; lidx++ ) {
	tmpval = cell[veg][band].layer[lidx].moist;
	if ( options.BINARY_STATE_FILE )
	  fwrite( &tmpval, sizeof(double), 1, filep->statefile );
	else
	  fprintf( filep->statefile, " %f", tmpval );
      }

      /* Write average ice content */
      for ( lidx = 0; lidx < options.Nlayer; lidx++ ) {
	for ( frost_area = 0; frost_area < options.Nfrost; frost_area++ ) {
	  tmpval = cell[veg][band].layer[lidx].ice[frost_area];
	  if ( options.BINARY_STATE_FILE ) {
	    fwrite( &tmpval, sizeof(double), 1, filep->statefile );
	  }
	  else {
	    fprintf( filep->statefile, " %f", tmpval );
	  }
	}
      }

      if ( veg < Nveg ) {
	/* Write dew storage */
	tmpval = veg_var[veg][band].Wdew;
	if ( options.BINARY_STATE_FILE )
	  fwrite( &tmpval, sizeof(double), 1, filep->statefile );
	else
	  fprintf( filep->statefile, " %f", tmpval );
        if (options.CARBON) {
	  /* Write cumulative NPP */
	  tmpval = veg_var[veg][band].AnnualNPP;
	  if ( options.BINARY_STATE_FILE )
	    fwrite( &tmpval, sizeof(double), 1, filep->statefile );
	  else
	    fprintf( filep->statefile, " %f", tmpval );
	  tmpval = veg_var[veg][band].AnnualNPPPrev;
	  if ( options.BINARY_STATE_FILE )
	    fwrite( &tmpval, sizeof(double), 1, filep->statefile );
	  else
	    fprintf( filep->statefile, " %f", tmpval );
	  /* Write soil carbon storages */
	  tmpval = cell[veg][band].CLitter;
	  if ( options.BINARY_STATE_FILE )
	    fwrite( &tmpval, sizeof(double), 1, filep->statefile );
	  else
	    fprintf( filep->statefile, " %f", tmpval );
	  tmpval = cell[veg][band].CInter;
	  if ( options.BINARY_STATE_FILE )
	    fwrite( &tmpval, sizeof(double), 1, filep->statefile );
	  else
	    fprintf( filep->statefile, " %f", tmpval );
	  tmpval = cell[veg][band].CSlow;
	  if ( options.BINARY_STATE_FILE )
	    fwrite( &tmpval, sizeof(double), 1, filep->statefile );
	  else
	    fprintf( filep->statefile, " %f", tmpval );
        }
      }
      
      /* Write snow data */
      if ( options.BINARY_STATE_FILE ) {
	fwrite( &snow[veg][band].last_snow, sizeof(int), 1, filep->statefile );
	fwrite( &snow[veg][band].MELTING, sizeof(char), 1, filep->statefile );
	fwrite( &snow[veg][band].coverage, sizeof(double), 1, filep->statefile );
	fwrite( &snow[veg][band].swq, sizeof(double), 1, filep->statefile );
	fwrite( &snow[veg][band].surf_temp, sizeof(double), 1, filep->statefile );
	fwrite( &snow[veg][band].surf_water, sizeof(double), 1, filep->statefile );
	fwrite( &snow[veg][band].pack_temp, sizeof(double), 1, filep->statefile );
	fwrite( &snow[veg][band].pack_water, sizeof(double), 1, filep->statefile );
	fwrite( &snow[veg][band].density, sizeof(double), 1, filep->statefile );
	fwrite( &snow[veg][band].coldcontent, sizeof(double), 1, filep->statefile );
	fwrite( &snow[veg][band].snow_canopy, sizeof(double), 1, filep->statefile );
      }
      else {
	fprintf( filep->statefile, " %i %i %f %f %f %f %f %f %f %f %f", 
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
	  fwrite( &energy[veg][band].T[nidx], sizeof(double), 1, 
		  filep->statefile );
	else
	  fprintf( filep->statefile, " %f", energy[veg][band].T[nidx] );

      if ( !options.BINARY_STATE_FILE ) fprintf( filep->statefile, "\n" );
      
    }
  }

  if ( options.LAKES ) {
    if ( options.BINARY_STATE_FILE ) {

      /* Write total soil moisture */
      for ( lidx = 0; lidx < options.Nlayer; lidx++ ) {
	fwrite( &lake_var.soil.layer[lidx].moist, sizeof(double), 1, filep->statefile );
      }

      /* Write average ice content */
      for ( lidx = 0; lidx < options.Nlayer; lidx++ ) {
        for ( frost_area = 0; frost_area < options.Nfrost; frost_area++ ) {
          fwrite( &lake_var.soil.layer[lidx].ice[frost_area], sizeof(double), 1, filep->statefile );
        }
      }
      if (options.CARBON) {
	/* Write soil carbon storages */
	tmpval = lake_var.soil.CLitter;
	if ( options.BINARY_STATE_FILE )
	  fwrite( &tmpval, sizeof(double), 1, filep->statefile );
	else
	  fprintf( filep->statefile, " %f", tmpval );
	tmpval = lake_var.soil.CInter;
	if ( options.BINARY_STATE_FILE )
	  fwrite( &tmpval, sizeof(double), 1, filep->statefile );
	else
	  fprintf( filep->statefile, " %f", tmpval );
	tmpval = lake_var.soil.CSlow;
	if ( options.BINARY_STATE_FILE )
	  fwrite( &tmpval, sizeof(double), 1, filep->statefile );
	else
	  fprintf( filep->statefile, " %f", tmpval );
      }

      /* Write snow data */
      fwrite( &lake_var.snow.last_snow, sizeof(int), 1, filep->statefile );
      fwrite( &lake_var.snow.MELTING, sizeof(char), 1, filep->statefile );
      fwrite( &lake_var.snow.coverage, sizeof(double), 1, filep->statefile );
      fwrite( &lake_var.snow.swq, sizeof(double), 1, filep->statefile );
      fwrite( &lake_var.snow.surf_temp, sizeof(double), 1, filep->statefile );
      fwrite( &lake_var.snow.surf_water, sizeof(double), 1, filep->statefile );
      fwrite( &lake_var.snow.pack_temp, sizeof(double), 1, filep->statefile );
      fwrite( &lake_var.snow.pack_water, sizeof(double), 1, filep->statefile );
      fwrite( &lake_var.snow.density, sizeof(double), 1, filep->statefile );
      fwrite( &lake_var.snow.coldcontent, sizeof(double), 1, filep->statefile );
      fwrite( &lake_var.snow.snow_canopy, sizeof(double), 1, filep->statefile );
      
      /* Write soil thermal node temperatures */
      for ( nidx = 0; nidx < options.Nnode; nidx++ ) 
	fwrite( &lake_var.energy.T[nidx], sizeof(double), 1, filep->statefile );

      /* Write lake-specific variables */
      fwrite( &lake_var.activenod, sizeof(int), 1, filep->statefile );
      fwrite( &lake_var.dz, sizeof(double), 1, filep->statefile );
      fwrite( &lake_var.surfdz, sizeof(double), 1, filep->statefile );
      fwrite( &lake_var.ldepth, sizeof(double), 1, filep->statefile );
      for ( node = 0; node <= lake_var.activenod; node++ ) {
        fwrite( &lake_var.surface[node], sizeof(double), 1, filep->statefile );
      }
      fwrite( &lake_var.sarea, sizeof(double), 1, filep->statefile );
      fwrite( &lake_var.volume, sizeof(double), 1, filep->statefile );
      for ( node = 0; node < lake_var.activenod; node++ ) {
        fwrite( &lake_var.temp[node], sizeof(double), 1, filep->statefile );
      }
      fwrite( &lake_var.tempavg, sizeof(double), 1, filep->statefile );
      fwrite( &lake_var.areai, sizeof(double), 1, filep->statefile );
      fwrite( &lake_var.new_ice_area, sizeof(double), 1, filep->statefile );
      fwrite( &lake_var.ice_water_eq, sizeof(double), 1, filep->statefile );
      fwrite( &lake_var.hice, sizeof(double), 1, filep->statefile );
      fwrite( &lake_var.tempi, sizeof(double), 1, filep->statefile );
      fwrite( &lake_var.swe, sizeof(double), 1, filep->statefile );
      fwrite( &lake_var.surf_temp, sizeof(double), 1, filep->statefile );
      fwrite( &lake_var.pack_temp, sizeof(double), 1, filep->statefile );
      fwrite( &lake_var.coldcontent, sizeof(double), 1, filep->statefile );
      fwrite( &lake_var.surf_water, sizeof(double), 1, filep->statefile );
      fwrite( &lake_var.pack_water, sizeof(double), 1, filep->statefile );
      fwrite( &lake_var.SAlbedo, sizeof(double), 1, filep->statefile );
      fwrite( &lake_var.sdepth, sizeof(double), 1, filep->statefile );
    }
    else {

      /* Write total soil moisture */
      for ( lidx = 0; lidx < options.Nlayer; lidx++ ) {
	fprintf( filep->statefile, " %f", lake_var.soil.layer[lidx].moist );
      }

      /* Write average ice content */
      for ( lidx = 0; lidx < options.Nlayer; lidx++ ) {
        for ( frost_area = 0; frost_area < options.Nfrost; frost_area++ ) {
          fprintf( filep->statefile, " %f", lake_var.soil.layer[lidx].ice[frost_area] );
        }
      }

      /* Write snow data */
      fprintf( filep->statefile, " %i %i %f %f %f %f %f %f %f %f %f", 
		 lake_var.snow.last_snow, (int)lake_var.snow.MELTING, 
		 lake_var.snow.coverage, lake_var.snow.swq, 
		 lake_var.snow.surf_temp, lake_var.snow.surf_water, 
		 lake_var.snow.pack_temp, lake_var.snow.pack_water, 
		 lake_var.snow.density, lake_var.snow.coldcontent, 
		 lake_var.snow.snow_canopy );
      
      /* Write soil thermal node temperatures */
      for ( nidx = 0; nidx < options.Nnode; nidx++ ) 
	fprintf( filep->statefile, " %f", lake_var.energy.T[nidx] );
      
      /* Write lake-specific variables */
      fprintf( filep->statefile, " %d", lake_var.activenod );
      fprintf( filep->statefile, " %f", lake_var.dz );
      fprintf( filep->statefile, " %f", lake_var.surfdz );
      fprintf( filep->statefile, " %f", lake_var.ldepth );
      for ( node = 0; node <= lake_var.activenod; node++ ) {
        fprintf( filep->statefile, " %f", lake_var.surface[node] );
      }
      fprintf( filep->statefile, " %f", lake_var.sarea );
      fprintf( filep->statefile, " %f", lake_var.volume );
      for ( node = 0; node < lake_var.activenod; node++ ) {
        fprintf( filep->statefile, " %f", lake_var.temp[node] );
      }
      fprintf( filep->statefile, " %f", lake_var.tempavg );
      fprintf( filep->statefile, " %f", lake_var.areai );
      fprintf( filep->statefile, " %f", lake_var.new_ice_area );
      fprintf( filep->statefile, " %f", lake_var.ice_water_eq );
      fprintf( filep->statefile, " %f", lake_var.hice );
      fprintf( filep->statefile, " %f", lake_var.tempi );
      fprintf( filep->statefile, " %f", lake_var.swe );
      fprintf( filep->statefile, " %f", lake_var.surf_temp );
      fprintf( filep->statefile, " %f", lake_var.pack_temp );
      fprintf( filep->statefile, " %f", lake_var.coldcontent );
      fprintf( filep->statefile, " %f", lake_var.surf_water );
      fprintf( filep->statefile, " %f", lake_var.pack_water );
      fprintf( filep->statefile, " %f", lake_var.SAlbedo );
      fprintf( filep->statefile, " %f", lake_var.sdepth );

      fprintf( filep->statefile, "\n" );
    }
  }

  if ( options.CARBON == CN_NORMAL || options.CARBON == CN_ADECOMP ) {
    for (band = 0; band < Nbands; band++ )
      {
	for(veg = 0; veg < MAX_PFT + 1; veg++)
	  {
	    if ( options.BINARY_STATE_FILE ) {
	      fwrite(&cn[band].dormant_flag[veg], sizeof(double), 1, filep->statefile);
	      fwrite(&cn[band].days_active[veg], sizeof(double), 1, filep->statefile);
	      fwrite(&cn[band].onset_flag[veg], sizeof(double), 1, filep->statefile);
	      fwrite(&cn[band].onset_counter[veg], sizeof(double), 1, filep->statefile);
	      fwrite(&cn[band].onset_gddflag[veg], sizeof(double), 1, filep->statefile);
	      fwrite(&cn[band].onset_fdd[veg], sizeof(double), 1, filep->statefile);
	      fwrite(&cn[band].onset_gdd[veg], sizeof(double), 1, filep->statefile);
	      fwrite(&cn[band].onset_swi[veg], sizeof(double), 1, filep->statefile);
	      fwrite(&cn[band].offset_flag[veg], sizeof(double), 1, filep->statefile);
	      fwrite(&cn[band].offset_counter[veg], sizeof(double), 1, filep->statefile);
	      fwrite(&cn[band].offset_fdd[veg], sizeof(double), 1, filep->statefile);
	      fwrite(&cn[band].offset_swi[veg], sizeof(double), 1, filep->statefile);
	      fwrite(&cn[band].lgsf[veg], sizeof(double), 1, filep->statefile);
	      fwrite(&cn[band].bglfr[veg], sizeof(double), 1, filep->statefile);
	      fwrite(&cn[band].bgtr[veg], sizeof(double), 1, filep->statefile);
	      fwrite(&cn[band].dayl[veg], sizeof(double), 1, filep->statefile);
	      fwrite(&cn[band].prev_dayl[veg], sizeof(double), 1, filep->statefile);
	      fwrite(&cn[band].annavg_t2m[veg], sizeof(double), 1, filep->statefile);
	      fwrite(&cn[band].tempavg_t2m[veg], sizeof(double), 1, filep->statefile);
	      fwrite(&cn[band].gpp[veg], sizeof(double), 1, filep->statefile);
	      fwrite(&cn[band].availc[veg], sizeof(double), 1, filep->statefile);
	      fwrite(&cn[band].xsmrpool_recover[veg], sizeof(double), 1, filep->statefile);
	      fwrite(&cn[band].alloc_pnow[veg], sizeof(double), 1, filep->statefile);
	      fwrite(&cn[band].c_allometry[veg], sizeof(double), 1, filep->statefile);
	      fwrite(&cn[band].n_allometry[veg], sizeof(double), 1, filep->statefile);
	      fwrite(&cn[band].plant_ndemand[veg], sizeof(double), 1, filep->statefile);
	      fwrite(&cn[band].tempsum_potential_gpp[veg], sizeof(double), 1, filep->statefile);
	      fwrite(&cn[band].annsum_potential_gpp[veg], sizeof(double), 1, filep->statefile);
	      fwrite(&cn[band].tempmax_retransn[veg], sizeof(double), 1, filep->statefile);
	      fwrite(&cn[band].annmax_retransn[veg], sizeof(double), 1, filep->statefile);
	      fwrite(&cn[band].avail_retransn[veg], sizeof(double), 1, filep->statefile);
	      fwrite(&cn[band].plant_nalloc[veg], sizeof(double), 1, filep->statefile);
	      fwrite(&cn[band].plant_calloc[veg], sizeof(double), 1, filep->statefile);
	      fwrite(&cn[band].excess_cflux[veg], sizeof(double), 1, filep->statefile);
	      fwrite(&cn[band].downreg[veg], sizeof(double), 1, filep->statefile);
	      fwrite(&cn[band].prev_leafc_to_litter[veg], sizeof(double), 1, filep->statefile);
	      fwrite(&cn[band].prev_frootc_to_litter[veg], sizeof(double), 1, filep->statefile);
	      fwrite(&cn[band].tempsum_npp[veg], sizeof(double), 1, filep->statefile);
	      fwrite(&cn[band].annsum_npp[veg], sizeof(double), 1, filep->statefile);
	      fwrite(&cn[band].leafc[veg], sizeof(double), 1, filep->statefile);
	      fwrite(&cn[band].leafc_storage[veg], sizeof(double), 1, filep->statefile);
	      fwrite(&cn[band].leafc_xfer[veg], sizeof(double), 1, filep->statefile);
	      fwrite(&cn[band].frootc[veg], sizeof(double), 1, filep->statefile);
	      fwrite(&cn[band].frootc_storage[veg], sizeof(double), 1, filep->statefile);
	      fwrite(&cn[band].frootc_xfer[veg], sizeof(double), 1, filep->statefile);
	      fwrite(&cn[band].livestemc[veg], sizeof(double), 1, filep->statefile);
	      fwrite(&cn[band].livestemc_storage[veg], sizeof(double), 1, filep->statefile);
	      fwrite(&cn[band].livestemc_xfer[veg], sizeof(double), 1, filep->statefile);
	      fwrite(&cn[band].deadstemc[veg], sizeof(double), 1, filep->statefile);
	      fwrite(&cn[band].deadstemc_storage[veg], sizeof(double), 1, filep->statefile);
	      fwrite(&cn[band].deadstemc_xfer[veg], sizeof(double), 1, filep->statefile);
	      fwrite(&cn[band].livecrootc[veg], sizeof(double), 1, filep->statefile);
	      fwrite(&cn[band].livecrootc_storage[veg], sizeof(double), 1, filep->statefile);
	      fwrite(&cn[band].livecrootc_xfer[veg], sizeof(double), 1, filep->statefile);
	      fwrite(&cn[band].deadcrootc[veg], sizeof(double), 1, filep->statefile);
	      fwrite(&cn[band].deadcrootc_storage[veg], sizeof(double), 1, filep->statefile);
	      fwrite(&cn[band].deadcrootc_xfer[veg], sizeof(double), 1, filep->statefile);
	      fwrite(&cn[band].gresp_storage[veg], sizeof(double), 1, filep->statefile);
	      fwrite(&cn[band].gresp_xfer[veg], sizeof(double), 1, filep->statefile);
	      fwrite(&cn[band].cpool[veg], sizeof(double), 1, filep->statefile);
	      fwrite(&cn[band].xsmrpool[veg], sizeof(double), 1, filep->statefile);
	      fwrite(&cn[band].pft_ctrunc[veg], sizeof(double), 1, filep->statefile);
	      fwrite(&cn[band].totvegc[veg], sizeof(double), 1, filep->statefile);
	      fwrite(&cn[band].leafn[veg], sizeof(double), 1, filep->statefile);
	      fwrite(&cn[band].leafn_storage[veg], sizeof(double), 1, filep->statefile);
	      fwrite(&cn[band].leafn_xfer[veg], sizeof(double), 1, filep->statefile);
	      fwrite(&cn[band].frootn[veg], sizeof(double), 1, filep->statefile);
	      fwrite(&cn[band].frootn_storage[veg], sizeof(double), 1, filep->statefile);
	      fwrite(&cn[band].frootn_xfer[veg], sizeof(double), 1, filep->statefile);
	      fwrite(&cn[band].livestemn[veg], sizeof(double), 1, filep->statefile);
	      fwrite(&cn[band].livestemn_storage[veg], sizeof(double), 1, filep->statefile);
	      fwrite(&cn[band].livestemn_xfer[veg], sizeof(double), 1, filep->statefile);
	      fwrite(&cn[band].deadstemn[veg], sizeof(double), 1, filep->statefile);
	      fwrite(&cn[band].deadstemn_storage[veg], sizeof(double), 1, filep->statefile);
	      fwrite(&cn[band].deadstemn_xfer[veg], sizeof(double), 1, filep->statefile);
	      fwrite(&cn[band].livecrootn[veg], sizeof(double), 1, filep->statefile);
	      fwrite(&cn[band].livecrootn_storage[veg], sizeof(double), 1, filep->statefile);
	      fwrite(&cn[band].livecrootn_xfer[veg], sizeof(double), 1, filep->statefile);
	      fwrite(&cn[band].deadcrootn[veg], sizeof(double), 1, filep->statefile);
	      fwrite(&cn[band].deadcrootn_storage[veg], sizeof(double), 1, filep->statefile);
	      fwrite(&cn[band].deadcrootn_xfer[veg], sizeof(double), 1, filep->statefile);
	      fwrite(&cn[band].retransn[veg], sizeof(double), 1, filep->statefile);
	      fwrite(&cn[band].npool[veg], sizeof(double), 1, filep->statefile);
	      fwrite(&cn[band].pft_ntrunc[veg], sizeof(double), 1, filep->statefile);
	    }

	    else

	      {

		fprintf(filep->statefile, "%d", veg);
		fprintf(filep->statefile, " %f", cn[band].dormant_flag[veg]);
		fprintf(filep->statefile, " %f", cn[band].days_active[veg]);
		fprintf(filep->statefile, " %f", cn[band].onset_flag[veg]);
		fprintf(filep->statefile, " %f", cn[band].onset_counter[veg]);
		fprintf(filep->statefile, " %f", cn[band].onset_gddflag[veg]);
		fprintf(filep->statefile, " %f", cn[band].onset_fdd[veg]);
		fprintf(filep->statefile, " %f", cn[band].onset_gdd[veg]);
		fprintf(filep->statefile, " %f", cn[band].onset_swi[veg]);
		fprintf(filep->statefile, " %f", cn[band].offset_flag[veg]);
		fprintf(filep->statefile, " %f", cn[band].offset_counter[veg]);
		fprintf(filep->statefile, " %f", cn[band].offset_fdd[veg]);
		fprintf(filep->statefile, " %f", cn[band].offset_swi[veg]);
		fprintf(filep->statefile, " %f", cn[band].lgsf[veg]);
		fprintf(filep->statefile, " %f", cn[band].bglfr[veg]);
		fprintf(filep->statefile, " %f", cn[band].bgtr[veg]);
		fprintf(filep->statefile, " %f", cn[band].dayl[veg]);
		fprintf(filep->statefile, " %f", cn[band].prev_dayl[veg]);
		fprintf(filep->statefile, " %f", cn[band].annavg_t2m[veg]);
		fprintf(filep->statefile, " %f", cn[band].tempavg_t2m[veg]);
		fprintf(filep->statefile, " %f", cn[band].gpp[veg]);
		fprintf(filep->statefile, " %f", cn[band].availc[veg]);
		fprintf(filep->statefile, " %f", cn[band].xsmrpool_recover[veg]);
		fprintf(filep->statefile, " %f", cn[band].alloc_pnow[veg]);
		fprintf(filep->statefile, " %f", cn[band].c_allometry[veg]);
		fprintf(filep->statefile, " %f", cn[band].n_allometry[veg]);
		fprintf(filep->statefile, " %f", cn[band].plant_ndemand[veg]);
		fprintf(filep->statefile, " %f", cn[band].tempsum_potential_gpp[veg]);
		fprintf(filep->statefile, " %f", cn[band].annsum_potential_gpp[veg]);
		fprintf(filep->statefile, " %f", cn[band].tempmax_retransn[veg]);
		fprintf(filep->statefile, " %f", cn[band].annmax_retransn[veg]);
		fprintf(filep->statefile, " %f", cn[band].avail_retransn[veg]);
		fprintf(filep->statefile, " %f", cn[band].plant_nalloc[veg]);
		fprintf(filep->statefile, " %f", cn[band].plant_calloc[veg]);
		fprintf(filep->statefile, " %f", cn[band].excess_cflux[veg]);
		fprintf(filep->statefile, " %f", cn[band].downreg[veg]);
		fprintf(filep->statefile, " %f", cn[band].prev_leafc_to_litter[veg]);
		fprintf(filep->statefile, " %f", cn[band].prev_frootc_to_litter[veg]);
		fprintf(filep->statefile, " %f", cn[band].tempsum_npp[veg]);
		fprintf(filep->statefile, " %f", cn[band].annsum_npp[veg]);
		fprintf(filep->statefile, " %f", cn[band].leafc[veg]);
		fprintf(filep->statefile, " %f", cn[band].leafc_storage[veg]);
		fprintf(filep->statefile, " %f", cn[band].leafc_xfer[veg]);
		fprintf(filep->statefile, " %f", cn[band].frootc[veg]);
		fprintf(filep->statefile, " %f", cn[band].frootc_storage[veg]);
		fprintf(filep->statefile, " %f", cn[band].frootc_xfer[veg]);
		fprintf(filep->statefile, " %f", cn[band].livestemc[veg]);
		fprintf(filep->statefile, " %f", cn[band].livestemc_storage[veg]);
		fprintf(filep->statefile, " %f", cn[band].livestemc_xfer[veg]);
		fprintf(filep->statefile, " %f", cn[band].deadstemc[veg]);
		fprintf(filep->statefile, " %f", cn[band].deadstemc_storage[veg]);
		fprintf(filep->statefile, " %f", cn[band].deadstemc_xfer[veg]);
		fprintf(filep->statefile, " %f", cn[band].livecrootc[veg]);
		fprintf(filep->statefile, " %f", cn[band].livecrootc_storage[veg]);
		fprintf(filep->statefile, " %f", cn[band].livecrootc_xfer[veg]);
		fprintf(filep->statefile, " %f", cn[band].deadcrootc[veg]);
		fprintf(filep->statefile, " %f", cn[band].deadcrootc_storage[veg]);
		fprintf(filep->statefile, " %f", cn[band].deadcrootc_xfer[veg]);
		fprintf(filep->statefile, " %f", cn[band].gresp_storage[veg]);
		fprintf(filep->statefile, " %f", cn[band].gresp_xfer[veg]);
		fprintf(filep->statefile, " %f", cn[band].cpool[veg]);
		fprintf(filep->statefile, " %f", cn[band].xsmrpool[veg]);
		fprintf(filep->statefile, " %f", cn[band].pft_ctrunc[veg]);
		fprintf(filep->statefile, " %f", cn[band].totvegc[veg]);
		fprintf(filep->statefile, " %f", cn[band].leafn[veg]);
		fprintf(filep->statefile, " %f", cn[band].leafn_storage[veg]);
		fprintf(filep->statefile, " %f", cn[band].leafn_xfer[veg]);
		fprintf(filep->statefile, " %f", cn[band].frootn[veg]);
		fprintf(filep->statefile, " %f", cn[band].frootn_storage[veg]);
		fprintf(filep->statefile, " %f", cn[band].frootn_xfer[veg]);
		fprintf(filep->statefile, " %f", cn[band].livestemn[veg]);
		fprintf(filep->statefile, " %f", cn[band].livestemn_storage[veg]);
		fprintf(filep->statefile, " %f", cn[band].livestemn_xfer[veg]);
		fprintf(filep->statefile, " %f", cn[band].deadstemn[veg]);
		fprintf(filep->statefile, " %f", cn[band].deadstemn_storage[veg]);
		fprintf(filep->statefile, " %f", cn[band].deadstemn_xfer[veg]);
		fprintf(filep->statefile, " %f", cn[band].livecrootn[veg]);
		fprintf(filep->statefile, " %f", cn[band].livecrootn_storage[veg]);
		fprintf(filep->statefile, " %f", cn[band].livecrootn_xfer[veg]);
		fprintf(filep->statefile, " %f", cn[band].deadcrootn[veg]);
		fprintf(filep->statefile, " %f", cn[band].deadcrootn_storage[veg]);
		fprintf(filep->statefile, " %f", cn[band].deadcrootn_xfer[veg]);
		fprintf(filep->statefile, " %f", cn[band].retransn[veg]);
		fprintf(filep->statefile, " %f", cn[band].npool[veg]);
		fprintf(filep->statefile, " %f", cn[band].pft_ntrunc[veg]);

	      }
	    
	    if ( !options.BINARY_STATE_FILE )
	      fprintf( filep->statefile, "\n" );
	  }

	if ( options.BINARY_STATE_FILE )
	  {
	    fwrite(&cn[band].decl, sizeof(double), 1, filep->statefile);
	    fwrite(&cn[band].fpi, sizeof(double), 1, filep->statefile);
	    fwrite(&cn[band].fpg, sizeof(double), 1, filep->statefile);
	    fwrite(&cn[band].annsum_counter, sizeof(double), 1, \
		   filep->statefile);
	    fwrite(&cn[band].cannsum_npp, sizeof(double), 1, \
		   filep->statefile);
	    fwrite(&cn[band].cannavg_t2m, sizeof(double), 1, \
		   filep->statefile);
	    for(nidx = 0; nidx < options.Nnode; nidx++)
	      {
		fwrite(&cn[band].watfc[nidx], sizeof(double), 1,	\
		   filep->statefile);
	      }
	    fwrite(&cn[band].me, sizeof(double), 1, filep->statefile);
	    fwrite(&cn[band].fire_prob, sizeof(double), 1, filep->statefile);
	    fwrite(&cn[band].mean_fire_prob, sizeof(double), 1, \
		   filep->statefile);
	    fwrite(&cn[band].fireseasonl, sizeof(double), 1, filep->statefile);
	    fwrite(&cn[band].ann_farea_burned, sizeof(double), 1, \
		   filep->statefile);
	    fwrite(&cn[band].cwdc, sizeof(double), 1, filep->statefile);
	    fwrite(&cn[band].litr1c, sizeof(double), 1, filep->statefile);
	    fwrite(&cn[band].litr2c, sizeof(double), 1, filep->statefile);
	    fwrite(&cn[band].litr3c, sizeof(double), 1, filep->statefile);
	    fwrite(&cn[band].soil1c, sizeof(double), 1, filep->statefile);
	    fwrite(&cn[band].soil2c, sizeof(double), 1, filep->statefile);
	    fwrite(&cn[band].soil3c, sizeof(double), 1, filep->statefile);
	    fwrite(&cn[band].soil4c, sizeof(double), 1, filep->statefile);
	    fwrite(&cn[band].seedc, sizeof(double), 1, filep->statefile);
	    fwrite(&cn[band].col_ctrunc, sizeof(double), 1, filep->statefile);
	    fwrite(&cn[band].totlitc, sizeof(double), 1, filep->statefile);
	    fwrite(&cn[band].totcolc, sizeof(double), 1, filep->statefile);
	    fwrite(&cn[band].prod10c, sizeof(double), 1, filep->statefile);
	    fwrite(&cn[band].prod100c, sizeof(double), 1, filep->statefile);
	    fwrite(&cn[band].cwdn, sizeof(double), 1, filep->statefile);
	    fwrite(&cn[band].litr1n, sizeof(double), 1, filep->statefile);
	    fwrite(&cn[band].litr2n, sizeof(double), 1, filep->statefile);
	    fwrite(&cn[band].litr3n, sizeof(double), 1, filep->statefile);
	    fwrite(&cn[band].soil1n, sizeof(double), 1, filep->statefile);
	    fwrite(&cn[band].soil2n, sizeof(double), 1, filep->statefile);
	    fwrite(&cn[band].soil3n, sizeof(double), 1, filep->statefile);
	    fwrite(&cn[band].soil4n, sizeof(double), 1, filep->statefile);
	    fwrite(&cn[band].sminn, sizeof(double), 1, filep->statefile);
	    fwrite(&cn[band].col_ntrunc, sizeof(double), 1, filep->statefile);
	    fwrite(&cn[band].seedn, sizeof(double), 1, filep->statefile);
	    fwrite(&cn[band].totcoln, sizeof(double), 1, filep->statefile);
	    fwrite(&cn[band].prod10n, sizeof(double), 1, filep->statefile);
	    fwrite(&cn[band].prod100n, sizeof(double), 1, filep->statefile);
	  }

	else

	  {
	    fprintf(filep->statefile, "%f", cn[band].decl);
	    fprintf(filep->statefile, " %f", cn[band].fpi);
	    fprintf(filep->statefile, " %f", cn[band].fpg);
	    fprintf(filep->statefile, " %f", cn[band].annsum_counter);
	    fprintf(filep->statefile, " %f", cn[band].cannsum_npp);
	    fprintf(filep->statefile, " %f", cn[band].cannavg_t2m);
	    for(nidx = 0; nidx < options.Nnode; nidx++)
	      {
		fprintf(filep->statefile, " %f", cn[band].watfc[nidx]);
	      }
	    fprintf(filep->statefile, " %f", cn[band].me);
	    fprintf(filep->statefile, " %f", cn[band].fire_prob);
	    fprintf(filep->statefile, " %f", cn[band].mean_fire_prob);
	    fprintf(filep->statefile, " %f", cn[band].fireseasonl);
	    fprintf(filep->statefile, " %f", cn[band].ann_farea_burned);
	    fprintf(filep->statefile, " %f", cn[band].cwdc);
	    fprintf(filep->statefile, " %f", cn[band].litr1c);
	    fprintf(filep->statefile, " %f", cn[band].litr2c);
	    fprintf(filep->statefile, " %f", cn[band].litr3c);
	    fprintf(filep->statefile, " %f", cn[band].soil1c);
	    fprintf(filep->statefile, " %f", cn[band].soil2c);
	    fprintf(filep->statefile, " %f", cn[band].soil3c);
	    fprintf(filep->statefile, " %f", cn[band].soil4c);
	    fprintf(filep->statefile, " %f", cn[band].seedc);
	    fprintf(filep->statefile, " %f", cn[band].col_ctrunc);
	    fprintf(filep->statefile, " %f", cn[band].totlitc);
	    fprintf(filep->statefile, " %f", cn[band].totcolc);
	    fprintf(filep->statefile, " %f", cn[band].prod10c);
	    fprintf(filep->statefile, " %f", cn[band].prod100c);
	    fprintf(filep->statefile, " %f", cn[band].cwdn);
	    fprintf(filep->statefile, " %f", cn[band].litr1n);
	    fprintf(filep->statefile, " %f", cn[band].litr2n);
	    fprintf(filep->statefile, " %f", cn[band].litr3n);
	    fprintf(filep->statefile, " %f", cn[band].soil1n);
	    fprintf(filep->statefile, " %f", cn[band].soil2n);
	    fprintf(filep->statefile, " %f", cn[band].soil3n);
	    fprintf(filep->statefile, " %f", cn[band].soil4n);
	    fprintf(filep->statefile, " %f", cn[band].sminn);
	    fprintf(filep->statefile, " %f", cn[band].col_ntrunc);
	    fprintf(filep->statefile, " %f", cn[band].seedn);
	    fprintf(filep->statefile, " %f", cn[band].totcoln);
	    fprintf(filep->statefile, " %f", cn[band].prod10n);
	    fprintf(filep->statefile, " %f", cn[band].prod100n);
	  }

	if ( !options.BINARY_STATE_FILE )
	  fprintf( filep->statefile, "\n" );

     }

  }

  /* Force file to be written */
  fflush(filep->statefile);

}

