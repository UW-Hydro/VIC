#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

#if LAKE_MODEL
#define MINLAYERTHICKNESS .25

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
  lake.tempi[MAXNOD]               Lake ice temp.
  lake.hice                Depth of lake ice.
  lake.fraci               Fractional coverage of lake ice. 
  lake.mixmax              Depth of local instability (node #).
  lake.volume
  lake.sarea

  modifications:
  11-18-02 Improvments made to the initialization of lake variables.  LCB

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

  for ( i = 0 ; i < lake_con.numnod ; i++ )
    {      
      lake->temp[i] = 15.0;
    }

  lake->tempi = 0.0;
  lake->hice = 0.0;
  lake->fraci = 0.0;
  lake->mixmax = 0;
  lake_snow->swq = 0.0;
  lake_snow->depth = 0.0;
  lake_snow->surf_temp = 0.0;

  /********************************************************************/
  /* Initialize lake physical parameters.                             */
  /********************************************************************/

  lake->ldepth = lake_con.depth_in;
  lake->dz = (lake->ldepth-SURF)/((float)(lake_con.numnod-1.));

  if(lake->dz < MINLAYERTHICKNESS) {
    lake->dz = MINLAYERTHICKNESS;
    lake->activenod = 1 + (int)(lake->ldepth-SURF)/lake->dz;
  }
  else
    lake->activenod = lake_con.numnod;

  // lake_con.basin equals the surface area at specific depths, lake->surface equals the
  // area at the top of each dynamic solution layer 
 
  for(k=0; k<lake->activenod;k++) {
    if(k==0)
      depth = SURF+lake->dz*(lake->activenod - 1);
    else
      depth = lake->dz*(lake->activenod - k);

    lake->surface[k] = get_sarea(lake_con, depth);
  }

  
  lake->volume=0.0;
   for(k=0; k<lake->activenod;k++) {
     if(k==0)
       lake->volume += (lake->surface[0] + lake->surface[1]) * SURF/2.;
     else if(k < lake->activenod-1)
       lake->volume += (lake->surface[k] + lake->surface[k+1]) * lake->dz/2.;
     else 
       lake->volume += (lake->surface[k]) * lake->dz;
    //  printf("area layer %d = %f, cum volume = %f\n",k,lake->surface[k],lake->volume);
   }
   
 
   lake->sarea = lake->surface[0];
 
  //  printf("initial volume = %f, initial depth = %f\n",lake->volume, lake->ldepth);
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
      fprintf(stderr, "Depth exceeds maximum depth.\n");
      sarea = lake_con.basin[0];
    }
  for(i=0; i< lake_con.numnod; i++)
    {
      if(i < lake_con.numnod -1 ) {
	if (depth <= lake_con.z[i] && depth > lake_con.z[i+1]) 
	  sarea = lake_con.basin[i+1] + (depth-lake_con.z[i+1])*(lake_con.basin[i] - lake_con.basin[i+1])/(lake_con.z[i] - lake_con.z[i+1]);
      }
      else {
	if(depth <= lake_con.z[i])
	  sarea = lake_con.basin[i];
      }
    }
  if(sarea == 0.0) {
    fprintf(stderr, "Somthing went wrong in get_sarea\n");
    exit(0);
  }

  return sarea;
}




#endif // LAKE_MODEL
