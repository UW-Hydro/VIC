#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
#include <string.h>

static char vcid[] = "$Id$";

void initialize_energy_bal (energy_bal_struct  **energy, 
                            cell_data_struct  ***cell,
                            soil_con_struct     *soil_con,
                            double               surf_temp,
			    double              *mu,
		            int                  veg_num,
                            int                  Nnodes,
			    int                  Ndist,
                            FILE                *fsoil)
/**********************************************************************
	initialize_energy_bal	Keith Cherkauer		July 30, 1996

  This routine initializes the energy balance data structure for use
  with the frozen soils model.

  Soil temperatures are initialized using a linear interpolation
  between the surface temperature (assumed = air temperature) and
  the average annual air temperature (assumed that bottom of 2nd soil
  layer is at const T = average annual air temp).

  UNITS: (m, s, kg, C, moisture in mm) unless otherwise specified

  variable modified:
	energy[veg].dz
	energy[veg].T
	energy[veg].fdepth
	cell[veg].layer[index].T
	cell[veg].layer[index].T_thaw
	cell[veg].layer[index].T_froz
	cell[veg].layer[index].moist
	cell[veg].layer[index].moist_thaw
	cell[veg].layer[index].moist_froz
	cell[veg].layer[index].ice
	cell[veg].layer[index].tdepth
	cell[veg].layer[index].fdepth
	cell[veg].layer[index].kappa
	cell[veg].layer[index].Cs

  Modifications:
  10-15-96  modified to reflect fixes in the frozen soils code	KAC
  07-10-97  modified to read water content in (m/m) instead of (mm) KAC
  09-18-97  modified to initialize with no files, and to initialize
		FULL_ENERGY case without FROZEN_SOIL		KAC
  06-24-98  modified to use distributed soil moisture only      KAC

**********************************************************************/
{
  extern option_struct options;
  extern debug_struct debug;

  char   tmpstr[MAXSTRING];
  char   ErrStr[MAXSTRING];
  int    i, j, veg, index, tmpdepth;
  int    tmpint;
  int    dry;
  int    band;
  int    zindex;
  double sum, Lsum, Zsum, dp, Ltotal;
  double *kappa, *Cs, *M;
  double moist[MAXlayer], ice[MAXlayer];
  double unfrozen, frozen;
  double *thermdepths;
  double **layer_ice;

  dp = soil_con[0].dp;
  Ltotal = 0;
  for(index=0;index<options.Nlayer;index++) Ltotal += soil_con[0].depth[index];

  /************************************************************************
    CASE 1: Frozen soils activated, and initial conditions files provided
  ************************************************************************/

  if(options.INIT_SOIL && options.FROZEN_SOIL) {

    thermdepths = (double *)calloc(Nnodes,sizeof(double));

    for ( veg = 0 ; veg <= veg_num ; veg++) {
      for(band=0;band<options.SNOW_BAND;band++) {
	rewind(fsoil);
	fgets(tmpstr,MAXSTRING,fsoil);
	
	fscanf(fsoil,"%*s %i\n",&tmpint);
	if(tmpint!=Nnodes) nrerror("Nnodes in soil initialization file, does not match that of the current model version");
	
	fscanf(fsoil,"%*s %lf %lf\n",&energy[veg][band].fdepth[0],
	       &energy[veg][band].fdepth[1]);
	fscanf(fsoil,"%*s");
	for(i=0;i<options.Nlayer;i++) fscanf(fsoil,"%lf",&moist[i]);
	fscanf(fsoil,"%*s");
	for(i=0;i<options.Nlayer;i++) fscanf(fsoil,"%lf",&ice[i]);
	
	fgets(tmpstr,MAXSTRING,fsoil);
	fgets(tmpstr,MAXSTRING,fsoil);
	for(j=0;j<Nnodes;j++) {
	  fscanf(fsoil,"%lf",&thermdepths[j]);
	  fscanf(fsoil,"%lf",&energy[veg][band].T[j]);
	}
	if(thermdepths[Nnodes-1] != dp) {
	  sprintf(ErrStr,"Thermal solution depth %i (Nnodes-1) must equal thermal damping depth %lf, but is equal to %lf",
		  Nnodes-1,dp,thermdepths[Nnodes-1]);
	  nrerror(ErrStr);
	}
	for(j=Nnodes-1;j>0;j--) {
	  thermdepths[j] -= thermdepths[j-1];
	  thermdepths[j] = (double)((int)(thermdepths[j] * 10000. + 0.5))
	    / 10000.;
	}
	for(j=1;j<Nnodes-1;j++) {
	  if((int)(thermdepths[j]*1000.+0.5) 
	     == (int)(thermdepths[j+1]*1000.+0.5))
	    energy[veg][band].dz[j] = (thermdepths[j] + thermdepths[j+1]) / 2.;
	  else {
	    energy[veg][band].dz[j] = thermdepths[j];
	    j++;
	    energy[veg][band].dz[j] = thermdepths[j+1];
	    if((int)((energy[veg][band].dz[j-1]
		      +energy[veg][band].dz[j])/2.*1000.+0.5)
	       != (int)(thermdepths[j]*1000.+0.5)) {
	      sprintf(ErrStr,"Check spacing between thermal layers %i and %i\n",
		      j-1,j);
	      nrerror(ErrStr);
	    }
	  }
	}
	energy[veg][band].dz[0] = thermdepths[1];
	energy[veg][band].dz[Nnodes-1] = thermdepths[Nnodes-1];

	sum = energy[veg][band].dz[0]/2. + energy[veg][band].dz[Nnodes-1]/2.;
	for(index=1;index<Nnodes-1;index++) sum += energy[veg][band].dz[index];

	for(dry=0;dry<Ndist;dry++) {
	  Lsum=0.;
	  for(index=0;index<options.Nlayer;index++) {
	    Lsum += soil_con[0].depth[index];
	    if(energy[veg][band].fdepth[1] 
	       > (Lsum - soil_con[0].depth[index])) {
	      if(energy[veg][band].fdepth[1] < Lsum) {
		cell[dry][veg][band].layer[index].tdepth 
		  = energy[veg][band].fdepth[1] 
		  - (Lsum - soil_con[0].depth[index]);
		cell[dry][veg][band].layer[index].moist_thaw = moist[index] 
		  * soil_con[0].depth[index] * 1000.;
		if(energy[veg][band].fdepth[0] < Lsum) {
		  cell[dry][veg][band].layer[index].fdepth 
		    = energy[veg][band].fdepth[0] 
		    - (Lsum - soil_con[0].depth[index]);
		  cell[dry][veg][band].layer[index].moist_froz 
		    = soil_con[0].depth[index]
		    * moist[index] * 1000.;
		  cell[dry][veg][band].layer[index].ice 
		    = soil_con[0].depth[index]
		    * ice[index] * 1000.;
		  cell[dry][veg][band].layer[index].moist 
		    = soil_con[0].depth[index] 
		    * moist[index] * 1000.;
		}
		else {
		  cell[dry][veg][band].layer[index].fdepth 
		    = soil_con[0].depth[index];
		  cell[dry][veg][band].layer[index].moist_froz 
		    = soil_con[0].depth[index]
		    * moist[index] * 1000.;
		  cell[dry][veg][band].layer[index].ice 
		    = soil_con[0].depth[index]
		    * ice[index] * 1000.;
		  cell[dry][veg][band].layer[index].moist = 0.;
		}
	      }
	      else {
		cell[dry][veg][band].layer[index].tdepth 
		  = soil_con[0].depth[index];
		cell[dry][veg][band].layer[index].fdepth 
		  = soil_con[0].depth[index];
		cell[dry][veg][band].layer[index].moist_thaw 
		  = soil_con[0].depth[index] * moist[index] * 1000.;
		cell[dry][veg][band].layer[index].moist_froz = 0.;
		cell[dry][veg][band].layer[index].ice = 0.;
		cell[dry][veg][band].layer[index].moist = 0.;
	      }
	    }
	    else if(energy[veg][band].fdepth[0] 
		    > (Lsum - soil_con[0].depth[index])) {
	      if(energy[veg][band].fdepth[0] < Lsum) {
		cell[dry][veg][band].layer[index].tdepth = 0.;
		cell[dry][veg][band].layer[index].fdepth 
		  = energy[veg][band].fdepth[0]
		  - (Lsum - soil_con[0].depth[index]);
		cell[dry][veg][band].layer[index].moist_thaw = 0.;
		cell[dry][veg][band].layer[index].moist_froz 
		  = soil_con[0].depth[index] * moist[index] * 1000.;
		cell[dry][veg][band].layer[index].ice 
		  = soil_con[0].depth[index]
		  * ice[index] * 1000.;
		cell[dry][veg][band].layer[index].moist 
		  = soil_con[0].depth[index] 
		  * moist[index] * 1000.;
	      }
	      else {
		cell[dry][veg][band].layer[index].tdepth = 0.;
		cell[dry][veg][band].layer[index].fdepth 
		  = soil_con[0].depth[index];
		cell[dry][veg][band].layer[index].moist_thaw = 0.;
		cell[dry][veg][band].layer[index].moist_froz 
		  = soil_con[0].depth[index] * moist[index] * 1000.;
		cell[dry][veg][band].layer[index].ice 
		  = soil_con[0].depth[index]
		  * ice[index] * 1000.;
		cell[dry][veg][band].layer[index].moist = 0.;
	      }
	    }
	    else {
	      cell[dry][veg][band].layer[index].tdepth = 0.;
	      cell[dry][veg][band].layer[index].fdepth = 0.;
	      cell[dry][veg][band].layer[index].moist_thaw = 0.;
	      cell[dry][veg][band].layer[index].moist_froz = 0.;
	      cell[dry][veg][band].layer[index].ice = 0.;
	      cell[dry][veg][band].layer[index].moist 
		= soil_con[0].depth[index] 
		* moist[index] * 1000.;
	    }
	  }
	  layer_ice = (double **)calloc(options.Nlayer,sizeof(double *));
	  for(i=0;i<options.Nlayer;i++) {
	    layer_ice[i] = (double *)calloc(3,sizeof(double));
	    layer_ice[i][0] = 0.;
	    layer_ice[i][1] = cell[dry][veg][band].layer[i].ice 
	      / (soil_con[0].depth[i] * 1000.);
	    layer_ice[i][2] = 0.;
	  }
	  distribute_soil_property(energy[veg][band].dz,
				   energy[veg][band].fdepth[0],
				   energy[veg][band].fdepth[1],
				   layer_ice,options.Nlayer,Nnodes,
				   soil_con[0].depth,energy[veg][band].ice);
	  for(i=0;i<options.Nlayer;i++) free((char *)layer_ice[i]);
	  free((char *)layer_ice);
	}
	
      }
    }

    free((char *)thermdepths);

  }

  /****************************************************************************
    CASE 2: Initialize soil if frozen soils not activated, and initialization
    file is provided.
  ****************************************************************************/

  else if(options.INIT_SOIL && options.FULL_ENERGY) {

    thermdepths = (double *)calloc(Nnodes,sizeof(double));

    for ( veg = 0 ; veg <= veg_num ; veg++) {
      for(band=0;band<options.SNOW_BAND;band++) {
	rewind(fsoil);
	fgets(tmpstr,MAXSTRING,fsoil);
	
	fscanf(fsoil,"%*s %i\n",&tmpint);
	if(tmpint!=Nnodes) nrerror("Nnodes in soil initialization file, does not match that of the current model version");
	
	energy[veg][band].fdepth[0]=energy[veg][band].fdepth[1]=0.;
	
	for(dry=0;dry<Ndist;dry++) {
	  for(i=0;i<options.Nlayer;i++) {
	    cell[dry][veg][band].layer[i].fdepth = 0.;
	    cell[dry][veg][band].layer[i].tdepth = 0.;
	  }
	}
	fscanf(fsoil,"%*s");
	for(i=0;i<options.Nlayer;i++) {
	  fscanf(fsoil,"%lf",&cell[WET][veg][band].layer[i].moist);
	  cell[WET][veg][band].layer[i].moist *= soil_con[0].depth[i]*1000.;
	  cell[WET][veg][band].layer[i].moist_froz = 0.;
	  cell[WET][veg][band].layer[i].moist_thaw = 0.;
	  if(options.DIST_PRCP) {
	    cell[DRY][veg][band].layer[i].moist 
	      = cell[WET][veg][band].layer[i].moist;
	    cell[DRY][veg][band].layer[i].moist_froz = 0.;
	    cell[DRY][veg][band].layer[i].moist_thaw = 0.;
	  }
	}
	fgets(tmpstr,MAXSTRING,fsoil);
	fgets(tmpstr,MAXSTRING,fsoil);
	for(j=0;j<Nnodes;j++) {
	  fscanf(fsoil,"%lf",&thermdepths[j]);
	  fscanf(fsoil,"%lf",&energy[veg][band].T[j]);
	}
	if(thermdepths[Nnodes-1] != dp) {
	  sprintf(ErrStr,"Thermal solution depth %i (Nnodes-1) must equal thermal damping depth %lf, but is equal to %lf",Nnodes-1,dp,thermdepths[Nnodes-1]);
	  nrerror(ErrStr);
	}
	for(j=Nnodes-1;j>0;j--) {
	  thermdepths[j] -= thermdepths[j-1];
	  thermdepths[j] = (double)((int)(thermdepths[j] * 10000. + 0.5))
	    / 10000.;
	}
	
	energy[veg][band].dz[0] = soil_con[0].depth[0];
	energy[veg][band].dz[1] = soil_con[0].depth[0];
	energy[veg][band].dz[2] = 2. * (dp - 1.5 * soil_con[0].depth[0]);
	energy[veg][band].dz[2] = 2. * (Ltotal - 1.5 * soil_con[0].depth[0]
					- energy[veg][band].dz[2]);
      }
    }
    
    free((char *)thermdepths);
    
  }

  /************************************************************************
    CASE 3: Initialize soil if frozen soils not activated, and no initial
    soil properties file given
  ************************************************************************/
    
  else if(options.FULL_ENERGY && !options.FROZEN_SOIL) {

    for ( veg = 0 ; veg <= veg_num ; veg++) {
      for(band=0;band<options.SNOW_BAND;band++) {
	energy[veg][band].fdepth[0] = 0.;
	energy[veg][band].fdepth[1] = 0.;
	for(dry=0;dry<Ndist;dry++) {
	  for(index=0;index<options.Nlayer;index++) {
	    cell[dry][veg][band].layer[index].fdepth = 0.;
	    cell[dry][veg][band].layer[index].tdepth = 0.;
	    cell[dry][veg][band].layer[index].T = 5.;
	    cell[dry][veg][band].layer[index].T_thaw = -999.;
	    cell[dry][veg][band].layer[index].T_froz = -999.;
	    cell[dry][veg][band].layer[index].moist_thaw = 0.;
	    cell[dry][veg][band].layer[index].moist_froz = 0.;
	    cell[dry][veg][band].layer[index].ice = 0.;
	  }
	}
	energy[veg][band].dz[0] = soil_con[0].depth[0];
	energy[veg][band].dz[1] = soil_con[0].depth[0];
	energy[veg][band].dz[2] = 2. * (dp - 1.5 * soil_con[0].depth[0]);
	energy[veg][band].T[0] = surf_temp;
	energy[veg][band].T[1] = surf_temp;
	energy[veg][band].T[2] = soil_con[0].avg_temp;
      }
    }
  }

  /*****************************************************************
    CASE 4: Initialize Energy Balance Variables if Frozen Soils 
    Activated, and no Initial Condition File Given 
  *****************************************************************/
  else if(options.FROZEN_SOIL) {
    for ( veg = 0 ; veg <= veg_num ; veg++) {
      for(band=0;band<options.SNOW_BAND;band++) {
	energy[veg][band].T[0] = surf_temp;
	energy[veg][band].dz[0] = soil_con[0].depth[0];
	energy[veg][band].dz[1] = soil_con[0].depth[0];
	energy[veg][band].T[Nnodes-1] = soil_con[0].avg_temp;
	energy[veg][band].T[1] = linear_interp(soil_con[0].depth[0],0.,
					       dp,surf_temp,
					       soil_con[0].avg_temp);
	
	/************************************************
	  Initialize Soil Layer Depths and Temperatures
	************************************************/
	Zsum = soil_con[0].depth[0];
	dp -= soil_con[0].depth[0] * 1.5;
	for(index=2;index<Nnodes;index++) {
	  energy[veg][band].dz[index] = dp/(((double)Nnodes-2.5));
	  Zsum += (energy[veg][band].dz[index]
		   +energy[veg][band].dz[index-1])/2.;
	  energy[veg][band].T[index] = linear_interp(Zsum,0.,soil_con[0].dp,
						     surf_temp,
						     soil_con[0].avg_temp);
	}
	dp += soil_con[0].depth[0] * 1.5;
	if((int)(Zsum*1000+0.5) != (int)(dp*1000+0.5)) {
	  sprintf(ErrStr,"Sum of thermal node thicknesses (%lf) in initialize_energy_bal do not equal dp (%lf), check initialization procedure",Zsum,dp);
	  nrerror(ErrStr);
	}
	
	for(dry=0;dry<Ndist;dry++) {
	  find_0_degree_fronts(&energy[veg][band],
			       cell[dry][veg][band].layer,Ltotal,
			       soil_con[0].depth,energy[veg][band].T,
			       options.Nlayer,
			       Nnodes);
	  
	  find_sublayer_temperatures(cell[dry][veg][band].layer,
				     energy[veg][band].T,
				     energy[veg][band].dz,
				     soil_con[0].depth,energy[veg][band].fdepth[0],
				     energy[veg][band].fdepth[1],
				     options.Nlayer,Nnodes);
	  
	  for(index=0;index<options.Nlayer;index++) {
	    if(cell[dry][veg][band].layer[index].tdepth > 0.)
	      cell[dry][veg][band].layer[index].moist_thaw 
		= cell[dry][veg][band].layer[index].moist;
	    if(cell[dry][veg][band].layer[index].tdepth 
	       == soil_con[0].depth[index])
	      cell[dry][veg][band].layer[index].moist = 0.;
	    else {
	      if(cell[dry][veg][band].layer[index].fdepth > 0.)
		cell[dry][veg][band].layer[index].moist_froz 
		  = cell[dry][veg][band].layer[index].moist;
	      if(cell[dry][veg][band].layer[index].fdepth 
		 == soil_con[0].depth[index])
		cell[dry][veg][band].layer[index].moist = 0.;
	    }
	    if(cell[dry][veg][band].layer[index].T_froz<0.
	       && cell[dry][veg][band].layer[index].T_froz != -999.) {
	      unfrozen 
		= maximum_unfrozen_water(cell[dry][veg][band].layer[index].T_froz,
					 soil_con[0].max_moist[index], 
					 soil_con[0].bubble,
					 soil_con[0].expt[index]);
	      if(unfrozen>soil_con[0].max_moist[index] || unfrozen<0.)
		unfrozen = soil_con[0].max_moist[index];
	      cell[dry][veg][band].layer[index].unfrozen = unfrozen;
	      
	      frozen = cell[dry][veg][band].layer[index].moist_froz - unfrozen;
	      if(frozen < 0.0) {
		frozen = 0.0;
		unfrozen = cell[dry][veg][band].layer[index].moist_froz;
	      }
	      cell[dry][veg][band].layer[index].ice = frozen;
	      cell[dry][veg][band].layer[index].moist_froz = unfrozen;
	    }
	    else if(cell[dry][veg][band].layer[index].T_froz == 0.) {
	      cell[dry][veg][band].layer[index].unfrozen 
		= soil_con[0].max_moist[index];
	      cell[dry][veg][band].layer[index].ice = 0.;
	    }
	    else if(cell[dry][veg][band].layer[index].T_froz != -999.)
	      nrerror("ERROR: Frozen Layer Temperature > 0C");
	    
	  }
	  
	  layer_ice = (double **)calloc(options.Nlayer,sizeof(double *));
	  for(i=0;i<options.Nlayer;i++) {
	    layer_ice[i] = (double *)calloc(3,sizeof(double));
	    layer_ice[i][0] = 0.;
	    layer_ice[i][1] = cell[dry][veg][band].layer[i].ice 
	      / (soil_con[0].depth[i] * 1000.);
	    layer_ice[i][2] = 0.;
	  }
	  distribute_soil_property(energy[veg][band].dz,
				   energy[veg][band].fdepth[0],
				   energy[veg][band].fdepth[1],
				   layer_ice,options.Nlayer,Nnodes,
				   soil_con[0].depth,energy[veg][band].ice);
	  for(i=0;i<options.Nlayer;i++) free((char *)layer_ice[i]);
	  free((char *)layer_ice);
	}
      }

    }
    
  }
  else {
    for ( veg = 0 ; veg <= veg_num ; veg++) {
      for(band=0;band<options.SNOW_BAND;band++) {
	for(index=0;index<options.Nlayer;index++) {
	  energy[veg][band].dz[index] = 1.;
	}
      }
    }
  }

  /** set initial soil properties for energy balance **/

  if(options.FULL_ENERGY || options.FROZEN_SOIL) {
    for ( veg = 0 ; veg <= veg_num ; veg++) {
      for(band=0;band<options.SNOW_BAND;band++) {
	if(soil_con[0].AreaFract[band]>0) {
	  index = 0;
	  Zsum = energy[veg][band].dz[index]/2.;
	  while(index<Nnodes 
		&& (int)((Zsum+energy[veg][band].dz[index+1]/2.)*1000.+0.5) 
		< (int)((soil_con[0].depth[0])*1000.+0.5)) {
	    index++;
	    Zsum += energy[veg][band].dz[index];
	  }
	  
	  if(index>=Nnodes 
	     || (int)((Zsum+energy[veg][band].dz[index+1]/2.)*1000.+0.5) 
	     > (int)((soil_con[0].depth[0])*1000.+0.5)) {
	    nrerror("One soil temperature depth must equal the depth of the top soil layer");
	  }
	  energy[veg][band].T1_index = index+1;
	  
	  for(dry=0;dry<Ndist;dry++) {
	    find_sublayer_temperatures(cell[dry][veg][band].layer,energy[veg][band].T,
				       energy[veg][band].dz,soil_con[0].depth,
				       energy[veg][band].fdepth[0],
				       energy[veg][band].fdepth[1],
				       options.Nlayer,Nnodes);
	    
	    kappa = NULL;
	    Cs = NULL;
	    M = NULL;
	    soil_thermal_calc(soil_con[0], cell[dry][veg][band].layer, 
			      energy[veg][band], kappa,
			      Cs, M, options.Nlayer, Nnodes);
	    
	    /** Save Thermal Conductivities for Energy Balance **/
	    if(dry==1) mu[veg] = 1. - mu[veg];
	    energy[veg][band].kappa[0] += cell[dry][veg][band].layer[0].kappa*mu[veg];
	    energy[veg][band].Cs[0] += cell[dry][veg][band].layer[0].Cs*mu[veg];
	    energy[veg][band].kappa[1] += cell[dry][veg][band].layer[1].kappa*mu[veg];
	    energy[veg][band].Cs[1] += cell[dry][veg][band].layer[1].Cs*mu[veg];
	  }
	  
	  if(energy[veg][band].fdepth[0]>0.) energy[veg][band].frozen=TRUE;
	}
      }  /** End loop through elevation bands **/
    }  /** End loop through vegetation types **/

    if(options.FROZEN_SOIL) {

      /***********************************************************
        Prepare soil constants for use in thermal flux solutions
      ***********************************************************/

      soil_con[0].dz_node        = (double *)calloc(Nnodes,sizeof(double));
      soil_con[0].expt_node      = (double *)calloc(Nnodes,sizeof(double));
      soil_con[0].max_moist_node = (double *)calloc(Nnodes,sizeof(double));
      soil_con[0].alpha          = (double *)calloc(Nnodes-2,sizeof(double));
      soil_con[0].beta           = (double *)calloc(Nnodes-2,sizeof(double));
      soil_con[0].gamma          = (double *)calloc(Nnodes-2,sizeof(double));

      Lsum = 0.;
      Zsum = 0.;
      index = 0;
      for(zindex=0;zindex<Nnodes;zindex++) {
	
	soil_con[0].dz_node[zindex] = energy[0][0].dz[zindex];
	if(Zsum+energy[0][0].dz[zindex]/2.>=Lsum+soil_con[0].depth[index]) {
	  soil_con[0].expt_node[zindex] 
	    = (soil_con[0].expt[index]*(Lsum+soil_con[0].depth[index]
				     - (Zsum-energy[0][0].dz[zindex]/2.))
	       + soil_con[0].expt[index+1]*(Zsum+energy[0][0].dz[zindex]/2.
					 - soil_con[0].depth[index]))
	    / energy[0][0].dz[zindex];
	  soil_con[0].max_moist_node[zindex] 
	    = (soil_con[0].max_moist[index]/soil_con[0].depth[index]/1000.
	       * (Lsum+soil_con[0].depth[index] 
		  - (Zsum-energy[0][0].dz[zindex]/2.))
	       + soil_con[0].max_moist[index+1]/soil_con[0].depth[index+1]/1000.
	       * (Zsum+energy[0][0].dz[zindex]/2. - soil_con[0].depth[index]))
	    / energy[0][0].dz[zindex];
	}
	else {
	  soil_con[0].expt_node[zindex] = soil_con[0].expt[index];
	  soil_con[0].max_moist_node[zindex] 
	    = soil_con[0].max_moist[index] / soil_con[0].depth[index] / 1000.;
	}
	if(zindex<Nnodes-1) {
	  Zsum += (energy[0][0].dz[zindex] + energy[0][0].dz[zindex+1]) / 2.;
	  if((int)(((Zsum - Lsum)+0.0005)*1000.) 
	     > (int)(soil_con[0].depth[index]*1000.)) {
	    Lsum += soil_con[0].depth[index];
	    index++;
	  }
	}
      }
      for(zindex=0;zindex<Nnodes-2;zindex++) {
	soil_con[0].alpha[zindex] = ((soil_con[0].dz_node[zindex+2] 
				   + soil_con[0].dz_node[zindex+1]) / 2.0 
				  + (soil_con[0].dz_node[zindex+1] 
				     + soil_con[0].dz_node[zindex]) / 2.0);
	soil_con[0].beta[zindex] = pow((soil_con[0].dz_node[zindex+2] 
				     + soil_con[0].dz_node[zindex+1])/2.0, 2.0) 
	  + pow((soil_con[0].dz_node[zindex+1]+soil_con[0].dz_node[zindex])/2.0, 
		2.0);
	soil_con[0].gamma[zindex] = ((soil_con[0].dz_node[zindex+2] 
				   + soil_con[0].dz_node[zindex+1]) / 2.0 
				  - (soil_con[0].dz_node[zindex+1] 
				     + soil_con[0].dz_node[zindex]) / 2.0);
      }
    }
  }  
}
