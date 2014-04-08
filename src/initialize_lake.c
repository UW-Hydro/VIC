#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

static char vcid[] = "$Id$";

int initialize_lake (lake_var_struct   *lake, 
		      lake_con_struct   lake_con,
		      soil_con_struct  *soil_con,
		      cell_data_struct *cell,
		      double            airtemp,
		      int               skip_hydro)

/**********************************************************************
	initialize_lake		Laura Bowling		March 8, 2000

  This routine initializes the lake variables for each new
  grid cell.

  VARIABLES INITIALIZED:
  lake.temp[MAXNOD]        Water temperature at each node.
  lake.tempi[MAXNOD]       Water temperature under ice at each node.
  lake.hice                Depth of lake ice.
  lake.areai               Area of lake ice. 
  lake.volume
  lake.sarea
  
  modifications:
  04-Oct-04 Merged with Laura Bowling's updated lake model code.	TJB
  23-Feb-05 Merged with Laura Bowling's second update to lake model code.	TJB
  2005-03-24 Added check for negative lake volumes.			TJB
  2006-Oct-16 Added RCS ID string.					TJB
  2006-Nov-07 Initialized aero_resist, aero_resist_used, and MELTING.	TJB
  2006-Nov-07 Removed LAKE_MODEL option.				TJB
  2007-Apr-23 Added initialization of lake->surface, lake->swe, and
	      lake->sdepth.						TJB
  2007-Oct-24 Changed get_sarea, get_volume, and get_depth to return exit
	      status so that errors can be trapped and communicated up the
	      chain of function calls.					KAC via TJB
  2007-Oct-24 Changed the error conditions so that get_depth does not
	      exit when depth == 0.0 (as long as volume == 0.0 when
	      depth == 0.0).						KAC via TJB
  2007-Nov-06 Replaced lake.fraci with lake.areai.  Added ice_depth()
	      function.							LCB via TJB
  2008-Jan-23 Added initialization of lake_snow->surf_temp, pack_water,
	      and pack_temp in conjunction with 2-layer snow pack over
	      lake ice.							LCB via TJB
  2008-Sep-09 Deleted initial volume print statement.			LCB via TJB
  2009-Jun-09 Lake_var data structure now only stores final (corrected)
	      values of aero_resist.					TJB
  2009-Jul-31 Removed references to lake_snow structure, which doesn't
	      exist outside of full_energy().				TJB
  2009-Sep-28 Added initialization of the new lake->snow, lake->soil,
	      and lake->energy structures.				TJB
  2009-Sep-30 Miscellaneous fixes for lake model.			TJB
  2009-Oct-08 Extended T fallback scheme to snow and ice T.		TJB
  2009-Dec-11 Removed min_liq and options.MIN_LIQ.			TJB
  2010-Sep-24 Added channel_in to store channel inflow separately from
	      incoming runoff from the catchment.			TJB
  2010-Nov-02 Added initialization of several lake_var variables that 
	      hadn't been initialized.					TJB
  2010-Nov-11 Added conditional skipping of initialization of moisture
	      fluxes, so that this function can be used to initialize
	      a lake that appears mid-timestep but whose moisture fluxes
	      have already been calculated.				TJB
  2010-Nov-21 Added lake->swe_save and lake->volume_save.		TJB
  2010-Nov-26 Added initialization of snow-related terms that are stored
	      in the lake_var structure.				TJB
  2011-Mar-01 Lake->soil state terms are now initialized to match those
	      of the cell data structure for the lake/wetland tile.	TJB
  2012-Jan-16 Removed LINK_DEBUG code					BN
  2012-Feb-07 Removed OUT_ZWT2 and OUT_ZWTL; renamed OUT_ZWT3 to
	      OUT_ZWT_LUMPED.						TJB
  2012-Feb-08 Renamed depth_full_snow_cover to max_snow_distrib_slope
	      and clarified the descriptions of the SPATIAL_SNOW
	      option.							TJB
  2013-Jul-25 Added soil carbon terms.					TJB
  2013-Dec-27 Moved SPATIAL_FROST to options_struct.			TJB
**********************************************************************/
{
  extern option_struct options;
  int i, k;
  int status;
  double depth;
  double remain;
  double in;
  double tmp_volume;

  /*  Assume no ice present, lake completely equilibrated with atmosphere. */

  for ( i = 0 ; i < MAX_LAKE_NODES; i++ ) {      
    lake->temp[i] = max(airtemp,0.0);
  }

  lake->areai = 0.0;
  lake->coldcontent = 0.0;
  lake->hice = 0.0;
  lake->ice_water_eq = 0.0;
  lake->new_ice_area = 0.0;
  lake->pack_temp = 0.0;
  lake->pack_water = 0.0;
  lake->SAlbedo = 0.0;
  lake->sdepth = 0.0;
  lake->surf_temp = 0.0;
  lake->surf_water = 0.0;
  lake->swe = 0.0;
  lake->swe_save = 0.0;
  lake->tempi = 0.0;

  /********************************************************************/
  /* Initialize lake physical parameters.                             */
  /********************************************************************/

  if (!skip_hydro) {

    lake->ldepth = lake_con.depth_in;

    if(lake->ldepth > MAX_SURFACE_LAKE && lake->ldepth < 2*MAX_SURFACE_LAKE) {
      /* Not quite enough for two full layers. */
      lake->surfdz = lake->ldepth/2.;
      lake->dz = lake->ldepth/2.;
      lake->activenod = 2;
    }
    else if(lake->ldepth >= 2* MAX_SURFACE_LAKE) {
      /* More than two layers. */	
      lake->surfdz = MAX_SURFACE_LAKE;
      lake->activenod = (int) (lake->ldepth/MAX_SURFACE_LAKE);
      if(lake->activenod > MAX_LAKE_NODES)
        lake->activenod = MAX_LAKE_NODES;
      lake->dz = (lake->ldepth-lake->surfdz)/((float)(lake->activenod-1));
    }
    else if(lake->ldepth > 0.0) {
      lake->surfdz = lake->ldepth;
      lake->dz = 0.0;
      lake->activenod = 1;
    }
    else {
      lake->surfdz = 0.0;
      lake->dz = 0.0;
      lake->activenod = 0;
      lake->ldepth = 0.0;
    }

    // lake_con.basin equals the surface area at specific depths as input by
    // the user in the lake parameter file or calculated in read_lakeparam(), 
    // lake->surface equals the area at the top of each dynamic solution layer 
 
    for(k=0; k<= lake->activenod; k++) {
      if(k==0)
        depth = lake->ldepth;
      else
        depth = lake->dz*(lake->activenod - k);
      status = get_sarea(lake_con, depth, &(lake->surface[k]));
      if (status < 0) {
        fprintf(stderr, "Error in get_sarea: record = %d, depth = %f, sarea = %e\n",0,depth,lake->surface[k]);
        return(status);
      }
    }

    lake->sarea = lake->surface[0];
    lake->sarea_save = lake->sarea;
    status = get_volume(lake_con, lake->ldepth, &tmp_volume);
    if (status < 0) {
      fprintf(stderr, "Error in get_volume: record = %d, depth = %f, volume = %e\n",0,depth,tmp_volume);
      return(status);
    }
    else if (status > 0) {
      fprintf(stderr, "Warning in get_volume: lake depth exceeds maximum; setting to maximum; record = %d\n",0);
    }
    lake->volume = tmp_volume+lake->ice_water_eq;
    lake->volume_save = lake->volume;
 
    // Initialize lake moisture fluxes to 0
    lake->baseflow_in=0.0;
    lake->baseflow_out=0.0;
    lake->channel_in=0.0;
    lake->evapw=0.0;
    lake->prec=0.0;
    lake->recharge=0.0;
    lake->runoff_in=0.0;
    lake->runoff_out=0.0;
    lake->snowmlt=0.0;
    lake->vapor_flux=0.0;

  } // if (!skip_hydro)

  // Initialize other miscellaneous lake properties
  lake->aero_resist = 0;
  for(k=0; k < lake->activenod; k++) {
    lake->density[k] = RHO_W;
  }

  // Initialize the snow, energy, and soil components of lake structure
  // If we implement heat flux between lake and underlying soil, we will need to initialize these more correctly
  // Snow state vars
  lake->snow.albedo            = 0.0;
  lake->snow.canopy_albedo     = 0.0;
  lake->snow.coldcontent       = 0.0;
  lake->snow.coverage          = 0.0;
  lake->snow.density           = 0.0;
  lake->snow.depth             = 0.0;
  lake->snow.last_snow         = MISSING;
  lake->snow.max_snow_depth    = 0.0;
  lake->snow.MELTING           = FALSE;
  lake->snow.pack_temp         = 0.0;
  lake->snow.pack_water        = 0.0;
  lake->snow.snow              = FALSE;
  lake->snow.snow_canopy       = 0.0;
  lake->snow.store_coverage    = 0.0;
  lake->snow.store_snow        = FALSE;
  lake->snow.store_swq         = 0.0;
  lake->snow.surf_temp         = 0.0;
  lake->snow.surf_temp_fbflag  = 0;
  lake->snow.surf_temp_fbcount = 0;
  lake->snow.surf_water        = 0.0;
  lake->snow.swq               = 0.0;
  lake->snow.snow_distrib_slope= 0.0;
  lake->snow.tmp_int_storage   = 0.0;
  // Snow fluxes
  lake->snow.blowing_flux      = 0.0;
  lake->snow.canopy_vapor_flux = 0.0;
  lake->snow.mass_error        = 0.0;
  lake->snow.melt              = 0.0;
  lake->snow.Qnet              = 0.0;
  lake->snow.surface_flux      = 0.0;
  lake->snow.transport         = 0.0;
  lake->snow.vapor_flux        = 0.0;
  // Energy state vars
  lake->energy.AlbedoLake       = 0.0;
  lake->energy.AlbedoOver       = 0.0;
  lake->energy.AlbedoUnder      = 0.0;
  lake->energy.frozen           = 0.0;
  lake->energy.Nfrost           = 0;
  lake->energy.Nthaw            = 0;
  lake->energy.T1_index         = 0;
  lake->energy.Tcanopy          = 0.0;
  lake->energy.Tcanopy_fbflag   = 0;
  lake->energy.Tcanopy_fbcount  = 0;
  lake->energy.Tfoliage         = 0.0;
  lake->energy.Tfoliage_fbflag  = 0;
  lake->energy.Tfoliage_fbcount = 0;
  lake->energy.Tsurf            = lake->temp[0];
  lake->energy.Tsurf_fbflag     = 0;
  lake->energy.Tsurf_fbcount    = 0;
  lake->energy.unfrozen         = 0.0;
  for (i=0; i<MAX_FRONTS; i++) {
    lake->energy.fdepth[i]      = 0.0;
    lake->energy.tdepth[i]      = 0.0;
  }
  for (i=0; i<2; i++) {
    lake->energy.Cs[i]          = 0.0;
    lake->energy.kappa[i]       = 0.0;
  }
  for (i=0; i<MAX_NODES; i++) {
    lake->energy.Cs_node[i]     = 0.0;
    lake->energy.ice[i]         = 0.0;
    lake->energy.kappa_node[i]  = 0.0;
    lake->energy.moist[i]       = 0.0;
    lake->energy.T[i]           = lake->temp[0];
    lake->energy.T_fbflag[i]    = 0;
    lake->energy.T_fbcount[i]   = 0;
  }
  // Energy fluxes
  lake->energy.advected_sensible = 0.0;
  lake->energy.advection         = 0.0;
  lake->energy.AtmosError        = 0.0;
  lake->energy.AtmosLatent       = 0.0;
  lake->energy.AtmosLatentSub    = 0.0;
  lake->energy.AtmosSensible     = 0.0;
  lake->energy.canopy_advection  = 0.0;
  lake->energy.canopy_latent     = 0.0;
  lake->energy.canopy_latent_sub = 0.0;
  lake->energy.canopy_refreeze   = 0.0;
  lake->energy.canopy_sensible   = 0.0;
  lake->energy.deltaCC           = 0.0;
  lake->energy.deltaH            = 0.0;
  lake->energy.error             = 0.0;
  lake->energy.fusion            = 0.0;
  lake->energy.grnd_flux         = 0.0;
  lake->energy.latent            = 0.0;
  lake->energy.latent_sub        = 0.0;
  lake->energy.longwave          = 0.0;
  lake->energy.LongOverIn        = 0.0;
  lake->energy.LongUnderIn       = 0.0;
  lake->energy.LongUnderOut      = 0.0;
  lake->energy.melt_energy       = 0.0;
  lake->energy.NetLongAtmos      = 0.0;
  lake->energy.NetLongOver       = 0.0;
  lake->energy.NetLongUnder      = 0.0;
  lake->energy.NetShortAtmos     = 0.0;
  lake->energy.NetShortGrnd      = 0.0;
  lake->energy.NetShortOver      = 0.0;
  lake->energy.NetShortUnder     = 0.0;
  lake->energy.out_long_canopy   = 0.0;
  lake->energy.out_long_surface  = 0.0;
  lake->energy.refreeze_energy   = 0.0;
  lake->energy.sensible          = 0.0;
  lake->energy.shortwave         = 0.0;
  lake->energy.ShortOverIn       = 0.0;
  lake->energy.ShortUnderIn      = 0.0;
  lake->energy.snow_flux         = 0.0;
  // Soil states and fluxes
  lake->soil.asat                = 1.0;
  if (!skip_hydro) {
    lake->soil.baseflow            = 0.0;
    lake->soil.inflow              = 0.0;
    lake->soil.runoff              = 0.0;
  }
  lake->soil.rootmoist           = 0.0;
  lake->soil.wetness             = 1.0;
  for (i=0; i<2; i++) {
    lake->soil.aero_resist[i]    = 0.0;
  }
  for (i=0; i<MAX_LAYERS; i++) {
    lake->soil.layer[i].Cs       = cell->layer[i].Cs;
    lake->soil.layer[i].T        = lake->temp[0];
    lake->soil.layer[i].evap     = 0.0;
    lake->soil.layer[i].kappa    = cell->layer[i].kappa;
    lake->soil.layer[i].moist    = soil_con->porosity[i]*soil_con->depth[i]*1000.;
    lake->soil.layer[i].phi      = cell->layer[i].phi;
    for (k=0; k<options.Nfrost; k++) {
      lake->soil.layer[i].ice[k]     = 0.0;
    }
  }
  lake->soil.zwt = 0.0;
  lake->soil.zwt_lumped = 0.0;
  if (!skip_hydro) {
    for (i=0; i<N_PET_TYPES; i++) {
      lake->soil.pot_evap[i]       = 0.0;
    }
  }
  if (options.CARBON) {
    lake->soil.RhLitter = 0.0;
    lake->soil.RhLitter2Atm = 0.0;
    lake->soil.RhInter = 0.0;
    lake->soil.RhSlow = 0.0;
    lake->soil.RhTot = 0.0;
    lake->soil.CLitter = 0.0;
    lake->soil.CInter = 0.0;
    lake->soil.CSlow = 0.0;
  }

  return(0);

}


int get_sarea(lake_con_struct lake_con, double depth, double *sarea)
/******************************************************************************
  Function to compute surface area of liquid water in the lake, given the
  current depth of liquid water.

  Modifications:
  2007-Oct-24 Added exit status.						TJB
    Exit status values:
          0: No errors
      ERROR: Error: area cannot be reconciled with given lake depth and nodes
******************************************************************************/
{
  int i;
  int status;

  status = 0;
  *sarea = 0.0;

  if (depth > lake_con.z[0]) {
    *sarea = lake_con.basin[0];
  }
  else {	
    for (i=0; i< lake_con.numnod; i++) {
      if (depth <= lake_con.z[i] && depth > lake_con.z[i+1]) 
	*sarea = lake_con.basin[i+1] + (depth-lake_con.z[i+1])*(lake_con.basin[i] - lake_con.basin[i+1])/(lake_con.z[i] - lake_con.z[i+1]);
    }
    if (*sarea == 0.0 && depth != 0.0) {
      status = ERROR;
    }
  }

  return status;

}

int get_volume(lake_con_struct lake_con, double depth, double *volume)
/******************************************************************************
  Function to compute liquid water volume stored within the lake basin, given
  the current depth of liquid water.

  Modifications:
  2007-Oct-24 Added exit status.						TJB
    Exit status values:
          0: No errors
          1: Warning: lake depth exceeds maximum; setting to maximum
      ERROR: Error: volume cannot be reconciled with given lake depth and nodes
******************************************************************************/
{
  int i;
  int status;
  double m, b;

  status = 0;
  *volume = 0.0;

  if (depth > lake_con.z[0]) {
    status = 1;
    *volume = lake_con.maxvolume;
  }

  for (i=lake_con.numnod-1; i>= 0; i--) {
    if (depth >= lake_con.z[i]) 
      *volume += (lake_con.basin[i] + lake_con.basin[i+1]) * (lake_con.z[i] - lake_con.z[i+1])/2.;
    else if (depth < lake_con.z[i] && depth >= lake_con.z[i+1]) {
      m = (lake_con.basin[i]-lake_con.basin[i+1])/(lake_con.z[i]-lake_con.z[i+1]);
      *volume += (depth - lake_con.z[i+1])*(m*(depth - lake_con.z[i+1])/2. + lake_con.basin[i+1]);
    }
  }

  if (*volume == 0.0  && depth != 0.0) {
    status = ERROR;
  }

  return status;

}

int get_depth(lake_con_struct lake_con, double volume, double *depth)
/******************************************************************************
  Function to compute the depth of liquid water in the lake (distance between
  surface and deepest point), given volume of liquid water currently stored in
  lake.

  Modifications:
  2007-Oct-24 Added exit status.						TJB
    Exit status values:
          0: No errors
          1: Warning: lake volume negative; setting to 0
      ERROR: Error: depth cannot be reconciled with given lake volume and nodes
  2007-Oct-30 Initialized surface area for lake bottom.				LCB via TJB
******************************************************************************/
{
  int k;
  int status;
  double m;
  double tempvolume;	

  status = 0;

  if (volume < -1*SMALL) {
    volume = 0.0;
    status = 1;
  }

  if (volume >= lake_con.maxvolume) {
    *depth = lake_con.maxdepth;
    *depth += (volume - lake_con.maxvolume)/lake_con.basin[0];	
  }
  else if ( volume < SMALL ) {
    *depth = 0.0;
  }
  else { 	
    // Update lake depth
    *depth = 0.0;
    tempvolume = volume;
    for ( k = lake_con.numnod - 1 ; k >= 0; k-- ) {
      if ( tempvolume > ((lake_con.z[k]-lake_con.z[k+1])
			 *(lake_con.basin[k]+lake_con.basin[k+1])/2.)) {
	// current layer completely filled
	tempvolume -= (lake_con.z[k]-lake_con.z[k+1])*(lake_con.basin[k]+lake_con.basin[k+1])/2.;
	*depth += lake_con.z[k] - lake_con.z[k+1];
      }
      else if (tempvolume > 0.0 ) {
        if (lake_con.basin[k]==lake_con.basin[k+1]) {
          *depth += tempvolume/lake_con.basin[k+1];
          tempvolume = 0.0;
	}
	else {
	  m = (lake_con.basin[k]-lake_con.basin[k+1])/(lake_con.z[k] - lake_con.z[k+1]);
	  *depth += ((-1*lake_con.basin[k+1]) + sqrt(lake_con.basin[k+1]*lake_con.basin[k+1] + 2.*m*tempvolume))/m;
	  tempvolume = 0.0;
	}
      }
    } 
    if (tempvolume/lake_con.basin[0] > SMALL )   {                  
      status = ERROR;
    }
  }

  if (*depth < 0.0 || (*depth == 0.0 && volume >= SMALL) ) {
    status = ERROR;
  }
  	
  return status;

}
