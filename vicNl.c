#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vicNl.h>
#include <global.h>

/** Main Program **/

void main(int argc, char *argv[])
/**********************************************************************
	vic2l.c		Dag Lohmann		January 1996

  This program controls file I/O and model steps for the water 
  balance version of the VIC-2L.

  For details about variables, input files and subroutines check:
	http://ce.washington.edu/~hydro/Lettenmaier/Models/VIC/VIC%20Water%20Balance%20Model.html

  UNITS: (cm, s, g, C) unless otherwise marked {at least in my subroutines}

  modifications:
  7-30-96  ifdef statements added to allow non working code to
           be cut out of the compilation step when a less complete
           version of the code is desired			KAC
  1-8-97   Began adding code to solve energy balance over bare soil.
           will be followed with 2-layer snow code, which is
           eventually to be made aware of the soil.		KAC

**********************************************************************/
{

  extern veg_lib_struct *veg_lib;
  extern option_struct options;
  extern debug_struct debug;
  extern Error_struct Error;

  /** Variable Declarations **/

  char    NEWCELL;
  char    LASTREC;
  int     STILL_STORM;
  int     rec, i;
  int     Ndist;
  int     Nveg_type;
  int     cellnum;
  int     index;
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

  /** Initialize Variables **/

  STILL_STORM = 0;	/* TRUE = currently precipitation */

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
  infiles.globalparam = open_file(filenames.global,"r");
  global_param = get_global_param(&filenames,infiles.globalparam);

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
  while((fscanf(infiles.soilparam, "%d", &flag))!=EOF) {
    soil_con = read_soilparam(infiles.soilparam);
    if(debug.PRT_SOIL) write_soilparam(soil_con); 
    if (flag == 1) {

      NEWCELL=TRUE;
      cellnum++;

      /** Read Grid Cell Vegetation Parameters **/
      veg_con = read_vegparam(infiles.vegparam, soil_con.gridcel,
                              Nveg_type);
      if(debug.PRT_VEGE) write_vegparam(veg_con); 

      /** Build Gridded Filenames, and Open **/
      builtnames = make_in_and_outfiles(&infiles, filenames, soil_con,
                   &outfiles);

      /** Read Forcing Data **/
      atmos = read_forcing_data(infiles, &global_param.nrecs,
				global_param.dt);
      if(debug.PRT_ATMOS) write_atmosdata(atmos, global_param.nrecs);

      /** Make Date Data Structure **/
      dmy      = make_dmy(global_param);

      /** Make Precipitation Distribution Control Structure **/
      prcp     = make_dist_prcp(veg_con[0].vegetat_type_num);
      for(i=0;i<Ndist;i++) {
        prcp.dist[i].veg_var  = make_veg_var(veg_con[0].vegetat_type_num);
        prcp.dist[i].cell     = make_cell_data(veg_con[0].vegetat_type_num+1,
                                options.Nlayer);
        if(options.FULL_ENERGY || options.SNOW_MODEL) {
          prcp.dist[i].snow = make_snow_data(global_param.nrecs);
          prcp.dist[i].energy = make_energy_bal(veg_con[0].vegetat_type_num+1,
                                &global_param.Ulayer,&global_param.Llayer);
        }
      }

      /** Initialize Soil and Vegetation Variables **/
      fprintf(stderr,"Initializing Variables\n");
      for(i=0;i<Ndist;i++) {
        initialize_soil(prcp.dist[i].cell,soil_con,
            veg_con[0].vegetat_type_num);
        initialize_veg(prcp.dist[i].veg_var,veg_con,global_param);
      }

      /**************************************************
        Initialize Meteological Forcing Values That
        Have not Been Specifically Set
      **************************************************/
      initialize_atmos(atmos,dmy,(double)soil_con.time_zone_lng,
		       (double)soil_con.lng,(double)soil_con.lat,
		       global_param.MAX_SNOW_TEMP,global_param.MIN_RAIN_TEMP,
		       global_param.nrecs,global_param.dt);

      if(!options.FULL_ENERGY)
        rad_and_vpd(atmos,soil_con,global_param.nrecs,dmy);

      /**************************************************
        Initialize Energy Balance and Snow Variables 
      **************************************************/
      if(options.FULL_ENERGY || options.FROZEN_SOIL) {
        fprintf(stderr,"Energy Balance Initialization\n");
        for(i=0;i<Ndist;i++)
          initialize_snow(prcp.dist[i].snow,veg_con[0].vegetat_type_num,
              infiles.init_snow);
        for(i=0;i<Ndist;i++)
          initialize_energy_bal(prcp.dist[i].energy,prcp.dist[i].cell,soil_con,
              atmos[0].air_temp,veg_con[0].vegetat_type_num,global_param.Ulayer,
              global_param.Llayer,infiles.init_soil);
      }

      fprintf(stderr,"Running Model\n");

      /** Update Error Handling Structure **/
      Error.outfp = outfiles;
      Error.infp = infiles;

      /***************************************************
	Intialize Moisture and Energy Balance Error Checks
	***************************************************/
      storage = 0.;
      for(index=0;index<options.Nlayer;index++) {
	storage += (prcp.dist[0].cell[0].layer[index].moist_thaw 
	         * prcp.dist[0].cell[0].layer[index].tdepth
	         / soil_con.depth[index]);
	storage += (prcp.dist[0].cell[0].layer[index].moist_froz 
	         * (prcp.dist[0].cell[0].layer[index].fdepth 
	         - prcp.dist[0].cell[0].layer[index].tdepth)
	         / soil_con.depth[index]);
	storage += (prcp.dist[0].cell[0].layer[index].moist
      	         * (soil_con.depth[index]
                 - prcp.dist[0].cell[0].layer[index].fdepth)
                 / soil_con.depth[index]);
	storage += (prcp.dist[0].cell[0].layer[index].ice 
	         * (prcp.dist[0].cell[0].layer[index].fdepth
	         - prcp.dist[0].cell[0].layer[index].tdepth)
	         / soil_con.depth[index]);
      }
      calc_water_balance_error(-1,0.,0.,storage);
      calc_energy_balance_error(-global_param.nrecs,0.,0.,0.,0.,0.);

      /******************************************
	Run Model in Grid Cell for all Time Steps
	******************************************/

      for ( rec = 0 ; rec < global_param.nrecs; rec++ ) {

        if(rec==global_param.nrecs-1) LASTREC = TRUE;
        else LASTREC = FALSE;

        dist_prec(&atmos[rec],&prcp,soil_con,veg_con,
                  dmy,global_param,outfiles,rec,cellnum,&STILL_STORM,
                  NEWCELL,LASTREC);
        NEWCELL=FALSE;

      }	/* End Rec Loop */

      close_files(infiles,outfiles,builtnames); 

/*****
      free((char *)Error.veg_con);
      free((char *)Error.atmos);
*****/

      free_dist_prcp(&prcp,veg_con[0].vegetat_type_num);
      free((char *)veg_con);
      free((char *)atmos);  
      free((char *)dmy);
      if(options.RADAR) free((char *)mu);

    }	/* End Run Flag Condition */
  } 	/* End Grid Loop */
}	/* End Main Program */
