#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vicNl.h>
#include <global.h>

static char vcid[] = "$Id$";

/** Main Program **/

void main(int argc, char *argv[])
/**********************************************************************
	vicNl.c		Dag Lohmann		January 1996

  This program controls file I/O and variable initialization as well as
  being the primary driver for the model.

  For details about variables, input files and subroutines check:
	http://ce.washington.edu/~hydro/Lettenmaier/Models/VIC/VIC_home.html

  UNITS: unless otherwise marked:
         all water balance components are in mm
	 all energy balance components are in mks
	 depths, and lengths are in m

  modifications:
  1997-98 Model was updated from simple 2 layer water balance to 
          an extension of the full energy and water balance 3 layer
	  model.                                                  KAC

**********************************************************************/
{

  extern veg_lib_struct *veg_lib;
  extern option_struct options;
  extern debug_struct debug;
  extern Error_struct Error;

  /** Variable Declarations **/

  char    NEWCELL;
  char    LASTREC;
  char    MODEL_DONE;
  int     rec, i;
  int     Ndist;
  int     Nveg_type;
  int     cellnum;
  int     index;
  int     RUN_MODEL;
  int     Ncells;
  int     cell_cnt;
  int     force_dt[2];
  double *mu;
  double  storage;
  dmy_struct *dmy;
  atmos_data_struct *atmos;
  veg_con_struct *veg_con;
  soil_con_struct soil_con;
  global_param_struct global_param;
  dist_prcp_struct prcp;	/* stores information about distributed 
                                   precipitation */
  filenames_struct filenames;
  filenames_struct builtnames;
  infiles_struct infiles;
  outfiles_struct outfiles;

  fprintf(stderr,"Running Model Version: %s\n",vcid);

  /** Read Model Options **/
  initialize_global();
  filenames = cmd_proc(argc, argv);

  /** Print Options **/
  fprintf(stderr,"Distributed Precipitation (%i)\n",options.DIST_PRCP);
  if(options.DIST_PRCP) {
    fprintf(stderr,"Use Radar Precipitation (%i)\n",options.RADAR);
  }
  fprintf(stderr,"Frozen Soils (%i)\n",options.FROZEN_SOIL);


  /** Read Global Control File **/
  force_dt[0] = force_dt[1] = -99;
  infiles.globalparam = open_file(filenames.global,"r");
  global_param = get_global_param(&filenames,infiles.globalparam,force_dt);

  /** Check and Open Files **/
  check_files(&infiles, filenames);
  open_debug();

  /** Read Vegetation Library File **/
  veg_lib = read_veglib(infiles.veglib,&Nveg_type);

  /** Initialize Parameters **/
  if(options.DIST_PRCP) Ndist = 2;
  else Ndist = 1;
  cellnum = -1;

  /************************************
    Run Model for all Active Grid Cells
    ************************************/
  MODEL_DONE = FALSE;
  cell_cnt=0;
  while(!MODEL_DONE) {
    if(!options.ARC_SOIL) {
      if((fscanf(infiles.soilparam, "%d", &flag))!=EOF) {
	if(flag) RUN_MODEL=TRUE;
	else     RUN_MODEL=FALSE;
      }
      else {
	MODEL_DONE = TRUE;
	RUN_MODEL = FALSE;
      }
      if(!MODEL_DONE) soil_con = read_soilparam(infiles.soilparam);
    }
    else {
      soil_con = read_soilparam_arc(infiles.soilparam, 
				    filenames.soil_dir, &Ncells, 
				    &RUN_MODEL, cell_cnt);
      cell_cnt++;
      if(cell_cnt==Ncells) MODEL_DONE = TRUE;
    }
    if(RUN_MODEL) {
      if(debug.PRT_SOIL) write_soilparam(soil_con); 

      NEWCELL=TRUE;
      cellnum++;

      /** Read Grid Cell Vegetation Parameters **/
      veg_con = read_vegparam(infiles.vegparam, soil_con.gridcel,
                              Nveg_type);
      calc_root_fractions(veg_con,soil_con);
      if(debug.PRT_VEGE) write_vegparam(veg_con); 

      /** Build Gridded Filenames, and Open **/
      builtnames = make_in_and_outfiles(&infiles, filenames, soil_con,
                   &outfiles);

      /** Read Forcing Data **/
      atmos = read_forcing_data(infiles, global_param.starthour,
				&global_param.nrecs,
				global_param.dt, force_dt);
      if(debug.PRT_ATMOS) write_atmosdata(atmos, global_param.nrecs);

      /** Read Elevation Band Data if Used **/
      read_snowband(infiles.snowband,rec,soil_con.gridcel,
		    (double)soil_con.elevation,
		    &soil_con.Tfactor,&soil_con.Pfactor,&soil_con.AreaFract);

      /** Make Date Data Structure **/
      dmy      = make_dmy(global_param);

      /** Make Precipitation Distribution Control Structure **/
      prcp     = make_dist_prcp(veg_con[0].vegetat_type_num, 
				&global_param.Nnodes);

      /** Initialize Soil and Vegetation Variables **/
      fprintf(stderr,"Initializing Variables\n");
      for(i=0;i<Ndist;i++) {
        initialize_soil(prcp.cell[i],soil_con,
			veg_con[0].vegetat_type_num);
        initialize_veg(prcp.veg_var[i],veg_con,global_param);
      }

      /**************************************************
         Initialize Meteological Forcing Values That
         Have not Been Specifically Set
       **************************************************/
      initialize_atmos(atmos,dmy,(double)soil_con.time_zone_lng,
		       (double)soil_con.lng,(double)soil_con.lat,
		       (double)soil_con.elevation,
		       global_param.MAX_SNOW_TEMP,global_param.MIN_RAIN_TEMP,
		       soil_con.Tfactor,global_param.nrecs,global_param.dt);

      /**************************************************
        Initialize Energy Balance and Snow Variables 
      **************************************************/
      if(options.SNOW_MODEL) {
        fprintf(stderr,"Snow Model Initialization\n");
	initialize_snow(prcp.snow,veg_con[0].vegetat_type_num,
			infiles.init_snow);
      }
      if(options.FULL_ENERGY) {
	fprintf(stderr,"Energy Balance Initialization\n");
	initialize_energy_bal(prcp.energy,prcp.cell,&soil_con,
			      atmos[0].air_temp,prcp.mu,
			      veg_con[0].vegetat_type_num,
			      global_param.Nnodes,
			      Ndist,infiles.init_soil);
      }

      fprintf(stderr,"Running Model\n");

      /** Update Error Handling Structure **/
      Error.outfp = outfiles;
      Error.infp = infiles;

      /***************************************************
	Intialize Moisture and Energy Balance Error Checks
	***************************************************/
      storage = 0.;
      for(index=0;index<options.Nlayer;index++)
	storage += find_total_layer_moisture(prcp.cell[0][0][0].layer[index],
					     soil_con.depth[index]);
      if(options.SNOW_MODEL) storage += prcp.snow[0][0].swq * 1000.;
      calc_water_balance_error(-global_param.nrecs,0.,0.,storage);
      calc_energy_balance_error(-global_param.nrecs,0.,0.,0.,0.,0.);

      /******************************************
	Run Model in Grid Cell for all Time Steps
	******************************************/

      for ( rec = 0 ; rec < global_param.nrecs; rec++ ) {

        if(rec==global_param.nrecs-1) LASTREC = TRUE;
        else LASTREC = FALSE;

        dist_prec(&atmos[rec],&prcp,soil_con,veg_con,
                  dmy,global_param,outfiles,rec,cellnum,
                  NEWCELL,LASTREC);
        NEWCELL=FALSE;

      }	/* End Rec Loop */

      close_files(infiles,outfiles,builtnames); 

      free_dist_prcp(&prcp,veg_con[0].vegetat_type_num);
      free_vegcon(&veg_con);
      free((char *)atmos);  
      free((char *)dmy);
      free((char *)soil_con.AreaFract);
      free((char *)soil_con.Tfactor);
      free((char *)soil_con.Pfactor);
      if(options.FROZEN_SOIL) {
	free((char*)soil_con.dz_node);
	free((char*)soil_con.expt_node);
	free((char*)soil_con.max_moist_node);
	free((char*)soil_con.alpha);
	free((char*)soil_con.beta);
	free((char*)soil_con.gamma);
      }
    }	/* End Run Model Condition */
  } 	/* End Grid Loop */
}	/* End Main Program */
