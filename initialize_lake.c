#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

#if LAKE_MODEL

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

**********************************************************************/
{
  extern option_struct options;
#if LINK_DEBUG
  extern debug_struct debug;
#endif

  int i, k;
  double remain;
  double in;

  /*  Assume no ice present, lake completely equilibrated with atmosphere. */

  for ( i = 0 ; i < lake_con.numnod ; i++ )
    {      
      lake->temp[i] = 6.0;     
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
  remain = modf(lake->ldepth/lake_con.dz, &in);
  remain += lake_con.dz - SURF;
  lake->activenod = in;

  lake->volume=0.0;
  for(k=0; k<lake->activenod;k++) {
    lake->surface[k] = lake_con.basin[k + lake_con.numnod - lake->activenod];
    if(k==0)
      lake->volume += lake->surface[k] * SURF;
    else 
      lake->volume += lake->surface[k] * lake_con.dz;
  }
  lake->volume += remain*lake->surface[0];
  lake->sarea = lake->surface[0];

  printf("Lake initial volume = %f m3, depth = %f m, surface area = %f m\n",lake->volume, lake->ldepth, lake->sarea);

  lake->runoff_out=0.0;
  lake->baseflow_out=0.0;
}

#endif // LAKE_MODEL
