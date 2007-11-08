#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

static char vcid[] = "$Id$";

int initialize_lake (lake_var_struct   *lake, 
		      lake_con_struct   lake_con,
		      snow_data_struct *lake_snow,
		      double            airtemp)

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
**********************************************************************/
{
  extern option_struct options;
#if LINK_DEBUG
  extern debug_struct debug;
#endif

  int i, k;
  int status;
  double depth;
  double remain;
  double in;

  /*  Assume no ice present, lake completely equilibrated with atmosphere. */

  for ( i = 0 ; i < MAX_LAKE_NODES; i++ ) {      
    lake->temp[i] = max(airtemp,0.0);
    lake->surface[i] = 0.0;
  }

  lake->tempi = 0.0;
  lake->hice = 0.0;
  lake->areai = .0;
  lake->new_ice_area = 0.0;
  lake->ice_water_eq = 0.0;
  lake->aero_resist = 0;
  lake->aero_resist_used = 0;
  lake_snow->swq = 0.0;
  lake_snow->depth = 0.0;
  lake_snow->surf_temp = 0.0;
  lake_snow->MELTING = FALSE;
  lake->swe = lake_snow->swq;
  lake->sdepth = lake_snow->depth;

  /********************************************************************/
  /* Initialize lake physical parameters.                             */
  /********************************************************************/

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
  status = get_volume(lake_con, lake->ldepth, &(lake->volume));
  if (status < 0) {
    fprintf(stderr, "Error in get_volume: record = %d, depth = %f, volume = %e\n",0,depth,lake->volume);
    return(status);
  }
  else if (status > 0) {
    fprintf(stderr, "Warning in get_volume: lake depth exceeds maximum; setting to maximum; record = %d\n",0);
  }
 
  printf("initial volume = %e km3, initial depth = %f m, initial area = %e km2\n",lake->volume/(1000.*1000.*1000.), lake->ldepth, lake->sarea/(1000.*1000.));
  lake->runoff_out=0.0;
  lake->baseflow_out=0.0;

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

int ice_depth(lake_con_struct lake_con, double volume, double ice_water_eq, double *hice)
/******************************************************************************
  Function to compute liquid water equivalent of lake ice (expressed in mm over
  the ice area), given the volume of liquid water currently stored in the lake
  and the current equivalent liquid water volume of lake ice.

  Modifications:
  2007-Oct-30 Created function so ice volume is state variable, not depth.	LCB via TJB
  2007-Oct-30 Added exit status.						TJB
    Exit status values:
          0: No errors
          1: Warning: ice volume negative; setting to 0
      ERROR: Error: ice depth cannot be reconciled with given ice volume
             and water equivalent
******************************************************************************/
{
  int k;
  double ldepth;
  double m;
  double tempvolume;
  double surfacearea;
  int status;

  status = 0;
  ldepth = 0.0;
  *hice = 0.0;

  if (ice_water_eq < 0.0) {
    ice_water_eq = 0.0;
    status = 1;
  }

  status = get_depth(lake_con, volume-ice_water_eq, &ldepth);
  if (status == ERROR) return(ERROR);
  status = get_sarea(lake_con, ldepth, &surfacearea);
  if (status == ERROR) return(ERROR);
  tempvolume = ice_water_eq* RHOICE/RHO_W;

  if(tempvolume >= lake_con.maxvolume) {
    *hice = lake_con.maxdepth;
    *hice += ((tempvolume - lake_con.maxvolume)/lake_con.basin[0]);
  }
  else {
  // Update ice depth

    for ( k = lake_con.numnod - 1 ; k >= 0; k-- ) {

      /* Start calculation at top of water layer */
      if(lake_con.z[k] > ldepth) {

        /* First layer on top  of water. */
        if(ldepth > lake_con.z[k+1]) {

          /* Ice volume fills up this layer. */
          if ( tempvolume >= ((lake_con.z[k]-ldepth) *(lake_con.basin[k]+surfacearea)/2.)) {
            tempvolume -= (lake_con.z[k]-ldepth)*(lake_con.basin[k]+surfacearea)/2.;
            *hice += lake_con.z[k] - ldepth;
          }
          /* Ice volume falls within this layer. */
          else {
            if(lake_con.basin[k] != lake_con.basin[k+1] ) {
              m = (lake_con.basin[k]-lake_con.basin[k+1])/(lake_con.z[k] - lake_con.z[k+1]);
              *hice += ((-1*surfacearea) + sqrt(surfacearea*surfacearea + 2.*m*tempvolume))/m;
            }
            else
              *hice += tempvolume/surfacearea;

            tempvolume = 0.0;
          }
        }
        /* Next layers, until all volume accounted for. */
        else if(tempvolume > 0.0) {

          /* Ice volume fills up this layer. */
          if ( tempvolume > ((lake_con.z[k]-lake_con.z[k+1]) *(lake_con.basin[k]+lake_con.basin[k+1])/2.)) {
            // current layer completely filled
            tempvolume -= (lake_con.z[k]-lake_con.z[k+1])*(lake_con.basin[k]+lake_con.basin[k+1])/2.;
            *hice += lake_con.z[k] - lake_con.z[k+1];
          }
          /* Ice volume falls within this layer. */
          else {
             if(lake_con.basin[k] != lake_con.basin[k+1] ) {
               m = (lake_con.basin[k]-lake_con.basin[k+1])/(lake_con.z[k] - lake_con.z[k+1]);
               *hice += ((-1*lake_con.basin[k+1]) + sqrt(lake_con.basin[k+1]*lake_con.basin[k+1] + 2.*m*tempvolume))/m;
             }
             else
               *hice += tempvolume/lake_con.basin[k];
            tempvolume = 0.0;
          }

        } // end if (ldepth > lake_con.z[k+1])

      } // end if (lake_con.z[k] > ldepth)

    } // end loop over nodes

  } // end if (tempvolume >= lake_con.maxvolume)

  if(tempvolume/lake_con.basin[0] > SMALL )   {
    status = ERROR;
  }
  else if(*hice <= 0.0  && ice_water_eq != 0.0) {
    status = ERROR;
  }
  else if(*hice < 0.0) *hice = 0.0;

  return(status);

}

