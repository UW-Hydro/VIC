#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

#if LAKE_MODEL

static char vcid[] = "$Id$";

void initialize_lake (lake_var_struct  *lake, 
		      lake_con_struct   lake_con,
		      snow_data_struct *lake_snow,
		      double            airtemp)

/**********************************************************************
	initialize_lake		Laura Bowling		March 8, 2000

  This routine initializes the lake variables for each new
  grid cell.

  VARIABLES INITIALIZED:
  lake.tp_in               Lake skin temperature (C). 
  lake.temp[MAXNOD]        Water temperature at each node.
  lake.tempi[MAXNOD]       Water temperature under ice at each node.
  lake.hice                Depth of lake ice.
  lake.fraci               Fractional coverage of lake ice. 
  lake.mixmax              Depth of local instability (node #).
  lake.volume
  lake.sarea
  
  modifications:
  04-Oct-04 Merged with Laura Bowling's updated lake model code.	TJB
  23-Feb-05 Merged with Laura Bowling's second update to lake model code.	TJB
  2005-03-24 Added check for negative lake volumes.			TJB
  2006-Oct-16 Added RCS ID string.					TJB
  2006-Nov-07 Initialized aero_resist, aero_resist_used, and MELTING.	TJB

**********************************************************************/
{
  extern option_struct options;
#if LINK_DEBUG
  extern debug_struct debug;
#endif

  int i, k;
  double depth;
  double remain;
  double in;

  /*  Assume no ice present, lake completely equilibrated with atmosphere. */

  for ( i = 0 ; i < MAX_LAKE_NODES; i++ )
    {      
      lake->temp[i] = max(airtemp,0.0);
    }

  lake->tempi = 0.0;
  lake->hice = 0.0;
  lake->fraci = .0;
  lake->mixmax = 0;
  lake->aero_resist = 0;
  lake->aero_resist_used = 0;
  lake_snow->swq = 0.0;
  lake_snow->depth = 0.0;
  lake_snow->surf_temp = 0.0;
  lake_snow->MELTING = FALSE;

  /********************************************************************/
  /* Initialize lake physical parameters.                             */
  /********************************************************************/

  lake->ldepth = lake_con.depth_in;

  if(lake->ldepth > MAX_SURFACE_LAKE && lake->ldepth < 2*MAX_SURFACE_LAKE) 		{
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
   else if(lake->ldepth > 0.0){
  	lake->surfdz = lake->ldepth;
	lake->dz = 0.0;
	lake->activenod = 1;
	}
  else {
	lake->surfdz = 0.0;
	lake->dz = 0.0;
	lake->activenod = 0;
	} 	

  // lake_con.basin equals the surface area at specific depths as input by
  // the user in the lake parameter file or calculated in read_lakeparam(), 
  // lake->surface equals the area at the top of each dynamic solution layer 
 
  for(k=0; k< lake->activenod; k++) {
    if(k==0)
      depth = lake->ldepth;
    else
      depth = lake->dz*(lake->activenod - k);

    lake->surface[k] = get_sarea(lake_con, depth);
  }

  lake->sarea = lake->surface[0];
  lake->volume = get_volume(lake_con, lake->ldepth);	  
 
  printf("initial volume = %e km3, initial depth = %f initial area = %e km2\n",lake->volume/(1000.*1000.*1000.), lake->ldepth, lake->sarea/(1000.*1000.));
  lake->runoff_out=0.0;
  lake->baseflow_out=0.0;
}


double get_sarea(lake_con_struct lake_con, double depth)
{
  int i;
  double sarea;

  sarea = 0.0;

  if(depth > lake_con.z[0])
    {
      sarea = lake_con.basin[0];
    }
  else {	
  for(i=0; i< lake_con.numnod; i++)
    {
     if (depth <= lake_con.z[i] && depth > lake_con.z[i+1]) 
	sarea = lake_con.basin[i+1] + (depth-lake_con.z[i+1])*(lake_con.basin[i] - lake_con.basin[i+1])/(lake_con.z[i] - lake_con.z[i+1]);
    }
  if(sarea == 0.0 && depth != 0.0) {
    fprintf(stderr, "Somthing went wrong in get_sarea\n");
    fprintf(stderr, "depth = %f, sarea = %e\n",depth,sarea);	
    exit(0);
  }
 }

  return sarea;
}

double get_volume(lake_con_struct lake_con, double depth)
{
  int i;
  double volume;	
  double m, b;

  volume = 0.0;

  if(depth > lake_con.z[0])
    {
      fprintf(stderr, "Depth exceeds maximum depth.\n");
      volume = lake_con.maxvolume;
    }
  for(i=lake_con.numnod-1; i>= 0; i--)
    {
      if (depth >= lake_con.z[i]) 
	  volume += (lake_con.basin[i] + lake_con.basin[i+1]) * (lake_con.z[i] - lake_con.z[i+1])/2.;
      else if (depth < lake_con.z[i] && depth >= lake_con.z[i+1]){
	m = (lake_con.basin[i]-lake_con.basin[i+1])/(lake_con.z[i]-lake_con.z[i+1]);
	volume += (depth - lake_con.z[i+1])*(m*(depth - lake_con.z[i+1])/2. + lake_con.basin[i+1]);
      }
    }
  if(volume == 0.0  && depth != 0.0) {
    fprintf(stderr, "Somthing went wrong in get_volume\n");
    exit(0);
  }

  return volume;
}

double get_depth(lake_con_struct lake_con, double volume)
{
  int k;
  double depth;
  double m;
  double tempvolume;	
  depth = 0.0;

  if (volume < 0.0) {
#if VERBOSE
    fprintf(stderr, "Warning: lake volume %e is negative; setting to 0.\n", volume);
#endif // VERBOSE
    volume = 0.0;
  }

  if(volume >= lake_con.maxvolume) {
    depth = lake_con.maxdepth;
    depth += (volume - lake_con.maxvolume)/lake_con.basin[0];	
  }
  else { 	
  // Update lake depth
    tempvolume = volume;	
    for ( k = lake_con.numnod - 1 ; k >= 0; k-- ) {
      if ( tempvolume > ((lake_con.z[k]-lake_con.z[k+1])
			 *(lake_con.basin[k]+lake_con.basin[k+1])/2.)) {
	// current layer completely filled
	tempvolume -= (lake_con.z[k]-lake_con.z[k+1])*(lake_con.basin[k]+lake_con.basin[k+1])/2.;
	depth += lake_con.z[k] - lake_con.z[k+1];
      }
      else if(tempvolume > 0.0 ) {
        if(lake_con.basin[k]==lake_con.basin[k+1]) {
           depth += tempvolume/lake_con.basin[k+1];
           tempvolume = 0.0;
	}
	else {
	  m = (lake_con.basin[k]-lake_con.basin[k+1])/(lake_con.z[k] - lake_con.z[k+1]);
	  depth += ((-1*lake_con.basin[k+1]) + sqrt(lake_con.basin[k+1]*lake_con.basin[k+1] + 2.*m*tempvolume))/m;
	  tempvolume = 0.0;
	}
      }
    } 

    if(tempvolume/lake_con.basin[0] > SMALL )   {                  
      fprintf(stderr, "Error in get depth tempvolume = %e, depth=%f\n", tempvolume,depth);
    }
  }

  if(depth <= 0.0  && volume != 0.0) {
    fprintf(stderr, "Somthing went wrong in get_depth\n");
    fprintf(stderr, "depth = %f, volume = %e\n",depth,volume);	
    exit(0);
  }
  	
  if(depth < 0.0) depth = 0.0;

  return depth;
}
#endif // LAKE_MODEL
