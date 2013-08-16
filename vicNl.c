#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vicNl.h>
#include <global.h>

static char vcid[] = "$Id: vicNl.c,v 4.2.2.15 2007/10/08 19:14:07 vicadmin Exp $";

/** Main Program **/

int main(int argc, char *argv[])
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
  03-12-03 Modifed to add AboveTreeLine to soil_con_struct so that
           the model can make use of the computed treeline.     KAC
  04-10-03 Modified to initialize storm parameters using the state
           file.                                                KAC
  04-10-03 Modified to start the model by skipping records until the
           state file date is found.  This replaces the previous method
           of modifying the global file start date, which can change 
           the interpolation of atmospheric forcing data.        KAC
  04-15-03 Modified to store wet and dry fractions when intializing 
           water balance storage.  This accounts for changes in model
           state initialization, which now stores wet and dry fractions
           rather than just averagedvalues.                      KAC
  29-Oct-03 Modified the version display banner to print the version
	    string defined in global.h.				TJB
  16-Jun-04 Modified to pass soil_con.avgJulyAirTemp to
	    initialize_atmos().					TJB
  2005-11-09 (Port from 4.1.0) Updated arglist for make_dist_prcp(), 
            as part of fix for QUICK_FLUX state file compatibility. GCT
  2006-Sep-01 (Port from 4.1.0) OUTPUT_FORCE option now calls close_files(). TJB
  2006-Sep-11 Implemented flexible output configuration; uses the new
              out_data and out_data_files structures. TJB
  2006-Sep-14 Implemented ALMA-compliant input and output; uses the
	      new save_data structure.  TJB
  2006-Oct-26 Merged infiles and outfiles structs into filep_struct;
	      This included merging builtnames into filenames.	TJB
  2007-Jan-15 Added PRT_HEADER option; added call to
	      write_header().					TJB
  2007-Apr-21 Added calls to free_dmy(), free_out_data_files(),
	      free_out_data(), and free_veglib().  Added closing of
	      all parameter files.				TJB
  2007-Sep-14 Excluded calls to free_veglib() and closing of parameter
	      files other than the soil param file for the case
	      when OUTPUT_FORCE=TRUE.				TJB
  2007-Oct-08 Fixed typo in call to check_state_file().  Was assigning
	      init_state file pointer to filep.statefile; now assigns
	      pointer to filep.init_state.			TJB

**********************************************************************/
{

  extern veg_lib_struct *veg_lib;
  extern option_struct options;
#if LINK_DEBUG
  extern debug_struct debug;
#endif
  extern Error_struct Error;
  extern global_param_struct global_param;

  /** Variable Declarations **/

  char                     NEWCELL;
  char                     LASTREC;
  char                     MODEL_DONE;
  char                     init_STILL_STORM;
  int                      rec, i, j;
  int                      veg;
  int                      dist;
  int                      band;
  int                      Ndist;
  int                      Nveg_type;
  int                      cellnum;
  int                      index;
  int                      init_DRY_TIME;
  int                      RUN_MODEL;
  int                      Ncells;
  int                      cell_cnt;
  int                      startrec;
  double                   storage;
  double                   veg_fract;
  double                   band_fract;
  int month,month_new,layer,veg_class,veg_class_tmp,iveg; //ingjerd dec 2008
  double irrig_current;
  double irrig_new;
  double noirrig_current;
  double noirrig_new;
  double soilmoist_irr_current;
  double soilmoist_noirr_current;
  double soilmoist_noirr_new;
  double soilmoist_irr_new;
  double cv_tot;
  double wdew_temp1,wdew_temp2;
  double snow_temp1,snow_temp2;
  double dummy1,dummy2;
  int iha;

  dmy_struct              *dmy;
  atmos_data_struct       *atmos;
  veg_con_struct          *veg_con;
  soil_con_struct          soil_con;
  dist_prcp_struct         prcp; /* stores information about distributed 
				    precipitation */
  filenames_struct         filenames;
  filep_struct             filep;
  out_data_file_struct     *out_data_files;
  out_data_struct          *out_data;
  save_data_struct         save_data;
  
  /** Read Model Options **/
  initialize_global();
  filenames = cmd_proc(argc, argv);

#if VERBOSE
  display_current_settings(DISP_VERSION,(filenames_struct*)NULL,(global_param_struct*)NULL);
#endif

  /** Read Global Control File **/
  filep.globalparam = open_file(filenames.global,"r");
  global_param = get_global_param(&filenames, filep.globalparam);

  /** Set up output data structures **/
  out_data = create_output_list();
  out_data_files = set_output_defaults(out_data);
  filep.globalparam = open_file(filenames.global,"r");
  parse_output_info(&filenames, filep.globalparam, &out_data_files, out_data);

  /** Check and Open Files **/
  check_files(&filep, &filenames);

#if !OUTPUT_FORCE

  /** Check and Open Debugging Files **/
#if LINK_DEBUG
  open_debug();
#endif

  /** Read Vegetation Library File **/
  veg_lib = read_veglib(filep.veglib,&Nveg_type);

#endif // !OUTPUT_FORCE

  /** Initialize Parameters **/
  if(options.DIST_PRCP) Ndist = 2;
  else Ndist = 1;
  cellnum = -1;

  /** Make Date Data Structure **/
  dmy      = make_dmy(&global_param);

  /** allocate memory for the atmos_data_struct **/
  alloc_atmos(global_param.nrecs, &atmos);

  startrec = 0;
#if !OUTPUT_FORCE
  if ( options.INIT_STATE ) 
    filep.init_state = check_state_file(filenames.init_state, dmy, 
					 &global_param, options.Nlayer, 
					 options.Nnode, &startrec);

  /** open state file if model state is to be saved **/
  if ( options.SAVE_STATE && strcmp( filenames.statefile, "NONE" ) != 0 ) 
    filep.statefile = open_state_file(&global_param, filenames, options.Nlayer, 
					 options.Nnode);
  else filep.statefile = NULL;
#endif // !OUTPUT_FORCE

  /************************************
    Run Model for all Active Grid Cells
    ************************************/
  MODEL_DONE = FALSE;
  cell_cnt=0;
  while(!MODEL_DONE) {
    if(!options.ARC_SOIL) {
      if((fscanf(filep.soilparam, "%d", &flag))!=EOF) {
	if(flag) RUN_MODEL=TRUE;
	else     RUN_MODEL=FALSE;
      }
      else {
	MODEL_DONE = TRUE;
	RUN_MODEL = FALSE;
      }
      if(!MODEL_DONE) soil_con = read_soilparam(filep.soilparam, RUN_MODEL);
    }
    else {
      soil_con = read_soilparam_arc(filep.soilparam, 
				    filenames.soil_dir, &Ncells, 
				    &RUN_MODEL, cell_cnt);
      cell_cnt++;
      if(cell_cnt==Ncells) MODEL_DONE = TRUE;
    }
    if(RUN_MODEL) {
#if LINK_DEBUG
      if(debug.PRT_SOIL) write_soilparam(&soil_con); 
#endif

#if QUICK_FS
      /** Allocate Unfrozen Water Content Table **/
      if(options.FROZEN_SOIL) {
	for(i=0;i<MAX_LAYERS;i++) {
	  soil_con.ufwc_table_layer[i] = (double **)malloc((QUICK_FS_TEMPS+1)*sizeof(double *));
	  for(j=0;j<QUICK_FS_TEMPS+1;j++) 
	    soil_con.ufwc_table_layer[i][j] = (double *)malloc(2*sizeof(double));
	}
	for(i=0;i<MAX_NODES;i++) {
	  soil_con.ufwc_table_node[i] = (double **)malloc((QUICK_FS_TEMPS+1)*sizeof(double *));

	  for(j=0;j<QUICK_FS_TEMPS+1;j++) 
	    soil_con.ufwc_table_node[i][j] = (double *)malloc(2*sizeof(double));
	}
      }
#endif /* QUICK_FS */

      NEWCELL=TRUE;
      cellnum++;

#if !OUTPUT_FORCE

      /** Read Grid Cell Vegetation Parameters **/
      veg_con = read_vegparam(filep.vegparam, soil_con.gridcel,
                              Nveg_type,soil_con.NRoots); /* soil_con.Nroots added ingjerd dec 2008 */ 
      calc_root_fractions(veg_con, &soil_con);
#if LINK_DEBUG
      if(debug.PRT_VEGE) write_vegparam(veg_con); 
#endif /* LINK_DEBUG*/
#endif /* !OUTPUT_FORCE */

      /** Build Gridded Filenames, and Open **/
      make_in_and_outfiles(&filep, &filenames, &soil_con, out_data_files);

      if (options.PRT_HEADER) {
        /** Write output file headers **/
        write_header(out_data_files, out_data, dmy, global_param);
      }

#if !OUTPUT_FORCE
      /** Read Elevation Band Data if Used **/
      read_snowband(filep.snowband,soil_con.gridcel,
		    (double)soil_con.elevation, &soil_con.Tfactor, 
		    &soil_con.Pfactor, &soil_con.AreaFract, 
		    &soil_con.AboveTreeLine);

      /** Make Precipitation Distribution Control Structure **/
      prcp     = make_dist_prcp(veg_con[0].vegetat_type_num);

#endif /* !OUTPUT_FORCE */


      /**************************************************
         Initialize Meteological Forcing Values That
         Have not Been Specifically Set
       **************************************************/

#if VERBOSE
      fprintf(stderr,"Initializing Forcing Data\n");
#endif /* VERBOSE */

      initialize_atmos(atmos, dmy, filep.forcing, 
		       (double)soil_con.time_zone_lng, (double)soil_con.lng,
		       (double)soil_con.lat, soil_con.elevation,
		       soil_con.annual_prec, global_param.wind_h, 
		       soil_con.rough, soil_con.avgJulyAirTemp, 
		       soil_con.Tfactor,
#if OUTPUT_FORCE
		       soil_con.AboveTreeLine, out_data_files, out_data); 
#else /* OUTPUT_FORCE */
                       soil_con.AboveTreeLine); 
#endif /* OUTPUT_FORCE */
      
#if !OUTPUT_FORCE
#if LINK_DEBUG
      if(debug.PRT_ATMOS) write_atmosdata(atmos, global_param.nrecs);
#endif

      /**************************************************
        Initialize Energy Balance and Snow Variables 
      **************************************************/

#if VERBOSE
      fprintf(stderr,"Model State Initialization\n");
#endif /* VERBOSE */
      initialize_model_state(&prcp, dmy[0], &global_param, filep,
			     soil_con.gridcel, veg_con[0].vegetat_type_num,
			     options.Nnode, Ndist, atmos[0].air_temp[NR], 
			     &soil_con, veg_con, &init_STILL_STORM,
			     &init_DRY_TIME, &save_data);


#if VERBOSE
      fprintf(stderr,"Running Model\n");
#endif /* VERBOSE */

      /** Update Error Handling Structure **/
      Error.filep = filep;
      Error.out_data_files = out_data_files;

      /***************************************************
	Intialize Moisture and Energy Balance Error Checks
        --- As of 4/15/03 this does not properly initialize
            storage from bands above treeline, when the model 
            state is restored from a file.  This can lead to 
            water balance errors in the initial time step but 
            does not impact the actual simulation.  It will
            be addressed in the next release version.  KAC
	***************************************************/
      storage = 0.;
      for ( veg = 0; veg <= veg_con[0].vegetat_type_num; veg++ ) {
	if ( veg < veg_con[0].vegetat_type_num ) veg_fract = veg_con[veg].Cv;
	else veg_fract = ( 1.0 - veg_con[0].Cv_sum );
	for ( band = 0; band < options.SNOW_BAND; band++ ) {
	  band_fract = soil_con.AreaFract[band];
	  if ( veg_fract > SMALL && band_fract > SMALL ) {
	    for(index=0;index<options.Nlayer;index++)
	      for ( dist = 0; dist < Ndist; dist ++ )
		storage += prcp.cell[dist][veg][band].layer[index].moist 
		  * veg_fract * band_fract;
	    storage += prcp.snow[veg][band].swq * 1000. * veg_fract 
	      * band_fract;
	    if ( veg != veg_con[0].vegetat_type_num ) {
	      for ( dist = 0; dist < Ndist; dist ++ ) 
		storage += prcp.veg_var[dist][veg][band].Wdew 
		  * veg_fract * band_fract;
	      storage += prcp.snow[veg][band].snow_canopy * 1000. 
		* veg_fract * band_fract;
	    }
	  }
	}
      }
      calc_water_balance_error(-global_param.nrecs,0.,0.,storage);
      calc_energy_balance_error(-global_param.nrecs,0.,0.,0.,0.,0.);

      /******************************************
	Run Model in Grid Cell for all Time Steps
	******************************************/

      for ( rec = startrec ; rec < global_param.nrecs; rec++ ) {

        if ( rec == global_param.nrecs - 1 ) LASTREC = TRUE;
        else LASTREC = FALSE;

        dist_prec( &atmos[rec], &prcp, &soil_con, veg_con, dmy, &global_param,
                   &filep, out_data_files, out_data, &save_data, rec, cellnum,
		   NEWCELL, LASTREC, init_STILL_STORM, init_DRY_TIME );
        NEWCELL=FALSE;
	init_DRY_TIME = -999;

	/* After time step simulation: adjust vegetation fractions and redistribute soil
	   moisture in case of irrigation and last day of month. ingjerd dec 2008. 
	   OBS! Works only for one snow band! Do only? if no snow on ground or in canopy.  */
	if( dmy[rec+1].day == 1 && dmy[rec+1].hour == 0 && options.IRRIGATION 
	    && veg_con[0].irrveg==1 && rec>1 && rec<global_param.nrecs-1) {
	    iveg=veg_con[0].vegetat_type_num;
	    veg_class_tmp = veg_con[iveg-1].veg_class;
	    veg_class = veg_lib[veg_class_tmp].veg_class;
	    cv_tot = veg_con[iveg-1].Cv+veg_con[iveg-2].Cv;
	    month = dmy[rec].month-1; //month starts at 0
	    month_new=month+1; if(month_new>11) month_new=0;
	    irrig_current = veg_lib[veg_class_tmp].irrpercent[month]*cv_tot; //fraction irrigated
	    irrig_new = veg_lib[veg_class_tmp].irrpercent[month_new]*cv_tot;
	    noirrig_current = (1-veg_lib[veg_class_tmp].irrpercent[month])*cv_tot; //fraction non irrigated
	    noirrig_new = (1-veg_lib[veg_class_tmp].irrpercent[month_new])*cv_tot;
	    //printf("month=%d cv=%f %f irrig_cv%f noirrig_cv%f irrig_new_cv%f noirr_new_cv%f\n",
	    //	   month,veg_con[iveg-1].Cv,veg_con[iveg-2].Cv,irrig_current,noirrig_current,irrig_new,noirrig_new);
	    dummy1=dummy2=0.;
	    for(layer = 0; layer < options.Nlayer; layer++) {
		soilmoist_noirr_current = prcp.cell[WET][iveg-2][0].layer[layer].moist; 
		soilmoist_irr_current = prcp.cell[WET][iveg-1][0].layer[layer].moist; 
		dummy1+=soilmoist_noirr_current*noirrig_current+soilmoist_irr_current*irrig_current;
		//if(rec==30 || rec==31) 
		// printf("vicNl_irrig rec%d day%d layer%d soilmoist_noirr:%f soilmoist_irr%f totalmoist:%f totalold(dummy1)%f\n",rec,dmy[rec].day,layer,soilmoist_noirr_current,soilmoist_irr_current,soilmoist_noirr_current*noirrig_current+soilmoist_irr_current*irrig_current,dummy1);
		
		if(noirrig_new>noirrig_current) { // irrigated fraction decreases
		    soilmoist_irr_new = soilmoist_irr_current;
		    soilmoist_noirr_new = 
			soilmoist_noirr_current + 
			(soilmoist_irr_current - soilmoist_noirr_current) * (irrig_current-irrig_new) / (noirrig_new);
//For some weird reason the next statement can't be commented out. If so, you get nan-values in some cells!
		    if(rec==30) {
			printf("vicNl decrease rec=%d layer=%d noirr_curr=%f irr_curr=%f noirr_new:%f irr_new:%f soilmoist_tot_current:%f soilmoist_new:%f\n",
			  rec,layer,soilmoist_noirr_current,soilmoist_irr_current,soilmoist_noirr_new,soilmoist_irr_new,soilmoist_irr_current*irrig_current+soilmoist_noirr_current*noirrig_current,soilmoist_irr_new*irrig_new+soilmoist_noirr_new*noirrig_new);
			//printf("vicNl_irrig decrease %d %f %f snow:%f snow_canopy:%f\n",rec,atmos[rec].air_temp[0],atmos[rec].prec[0],prcp.snow[iveg-1][0].swq,prcp.snow[iveg-1][0].snow_canopy);  
			  }
		}
		else { // irrigated fraction increases or stays constant
		    soilmoist_noirr_new = soilmoist_noirr_current;
		    soilmoist_irr_new = 
			soilmoist_irr_current + 
			(soilmoist_noirr_current - soilmoist_irr_current) * (irrig_new-irrig_current) / (irrig_new);


		    //     if(rec==30 || rec==31) {
		    //printf("vicNl_irrig increase/constant %d %d %f %f %f %f soilmoist_current:%f soilmoist_adjusted:%f\n",
		    //    rec,layer,soilmoist_noirr_current,soilmoist_irr_current,soilmoist_noirr_new,soilmoist_irr_new,soilmoist_irr_current*irrig_current+soilmoist_noirr_current*noirrig_current,soilmoist_irr_new*irrig_new+soilmoist_noirr_new*noirrig_new);
		//   printf("vicNl_irrig increase/constant %d %d %f %f %f %f soilmoist_old:%f soilmoist_new:%f\n",
		//      rec,layer,soilmoist_noirr_current,soilmoist_irr_current,soilmoist_noirr_new,soilmoist_irr_new,soilmoist_irr_current*irrig_current+soilmoist_noirr_current*noirrig_current,soilmoist_irr_new*irrig_new+soilmoist_noirr_new*noirrig_new);
		//	}
		}
	   
		if(soilmoist_noirr_new<0 || soilmoist_irr_new<0) printf("vicNl_irrig NEGATIVE!!! %d\n",rec);
		if((soilmoist_irr_current*irrig_current+soilmoist_noirr_current*noirrig_current+0.0001)<=(soilmoist_irr_new*irrig_new+soilmoist_noirr_new*noirrig_new) || (soilmoist_irr_current*irrig_current+soilmoist_noirr_current*noirrig_current-0.0001)>=(soilmoist_irr_new*irrig_new+soilmoist_noirr_new*noirrig_new))  printf("Hei %d %f %f \n",rec,soilmoist_irr_current*irrig_current+soilmoist_noirr_current*noirrig_current,soilmoist_irr_new*irrig_new+soilmoist_noirr_new*noirrig_new); 

	      prcp.cell[WET][iveg-2][0].layer[layer].moist=soilmoist_noirr_new; 
	      prcp.cell[WET][iveg-1][0].layer[layer].moist=soilmoist_irr_new; 
	      dummy2+=prcp.cell[WET][iveg-2][0].layer[layer].moist*noirrig_new+prcp.cell[WET][iveg-1][0].layer[layer].moist*irrig_new;
	      if(fabs(dummy1-dummy2)>0.001) 
	          printf("vicNl_irrig soilmoist difference layer%d totalold%f totalnew(dummy2)%f\n",layer,dummy1,dummy2);
	
	     } //er det mulig at du kan få negative verdier? 

            //adjust Wdew and snow. assuming wdew on irrigated veg is always higher than on nonirr veg
	    wdew_temp1=prcp.veg_var[WET][iveg-2][0].Wdew; //nonirrigated
	    wdew_temp2=prcp.veg_var[WET][iveg-1][0].Wdew; //irrigated
	    snow_temp1=prcp.snow[iveg-2][0].swq; //nonirrigated
	    snow_temp2=prcp.snow[iveg-1][0].swq; //irrigated

	    /*if(rec==30 || rec==31) {
	     printf("vicNl_irrig wdew end of month: rec:%d month:%d day:%d wdew_noirrig:%f wdew_irrig:%f total wdew storage:%f current nonirrfraction: %f current irrfraction:%f\n",
	    	 rec,dmy[rec].month,dmy[rec].day,wdew_temp1,wdew_temp2,wdew_temp1*noirrig_current+wdew_temp2*irrig_current,noirrig_current,irrig_current);
	     printf("vicNl_irrig snow end of month: rec:%d month:%d day:%d snow_noirrig:%f snow_irrig:%f snowcanopy:%f %f\n",
	     rec,dmy[rec].month,dmy[rec].day,snow_temp1,snow_temp2,prcp.snow[iveg-1][0].snow_canopy,prcp.snow[iveg-2][0].snow_canopy);
	     }*/

	     if(noirrig_new>noirrig_current) { // irrigated fraction decreases
		 wdew_temp2 = prcp.veg_var[WET][iveg-1][0].Wdew; //irrig
		 wdew_temp1 = wdew_temp1 + 
		     (wdew_temp2 - wdew_temp1) * (noirrig_new - noirrig_current) / (noirrig_new);
		 snow_temp2 = prcp.snow[iveg-1][0].swq; //irrig
		 snow_temp1 = snow_temp1 + 
		     (snow_temp2 - snow_temp1) * (noirrig_new - noirrig_current) / (noirrig_new);
	     }
	     else { //irrigated fraction increases
		 wdew_temp1 = prcp.veg_var[WET][iveg-2][0].Wdew; //noirrig
		 wdew_temp2 =  wdew_temp2 * irrig_current / irrig_new + 
		     wdew_temp1*(irrig_new-irrig_current)/irrig_new;
		 snow_temp1 = prcp.snow[iveg-2][0].swq; //noirrig
		 snow_temp2 =  snow_temp2 * irrig_current / irrig_new + 
		     snow_temp1*(irrig_new-irrig_current)/irrig_new;
	     }
	     prcp.veg_var[WET][iveg-2][0].Wdew=wdew_temp1;
	     prcp.veg_var[WET][iveg-1][0].Wdew=wdew_temp2;
	     prcp.snow[iveg-2][0].swq=snow_temp1;
	     prcp.snow[iveg-1][0].swq=snow_temp2;

	     //if(rec==30 || rec==31) 
	     //printf("vicNl_irrig wdew start of next month: rec:%d month:%d day:%d wdew_nonirrig:%f wdew_irrig:%f total wdew storage:%f next month's nonirrfraction: %f next month's irrfraction:%f\n",
	     //	 rec,dmy[rec].month,dmy[rec].day,wdew_temp1,wdew_temp2,wdew_temp1*noirrig_new+wdew_temp2*irrig_new,noirrig_new,irrig_new);

	      //adjust fraction irr/noirr vegetation //
	      //veg_con[iveg-2].Cv=noirrig_new;
	      //veg_con[iveg-1].Cv=irrig_new;
	      //printf("vicNl_irrig rec:%d month:%d cv1:%f cv2:%f \n",
	      //	     rec,dmy[rec].month,veg_con[iveg-1].Cv,veg_con[iveg-2].Cv);

	      //if(rec==30) 
	      // printf("vicNl_irrig soilmoist next month: rec:%d month:%d day:%d %f \n",
	      //	 rec,dmy[rec].month,dmy[rec].day,(prcp.cell[WET][iveg-2][0].layer[0].moist+prcp.cell[WET][iveg-2][0].layer[1].moist+prcp.cell[WET][iveg-2][0].layer[2].moist)*noirrig_new+(prcp.cell[WET][iveg-1][0].layer[0].moist+prcp.cell[WET][iveg-1][0].layer[1].moist+prcp.cell[WET][iveg-2][0].layer[2].moist)*irrig_new);	      
	} // end adjust veg fractions if irrigation

      }	/* End Rec Loop */

#endif /* !OUTPUT_FORCE */

      close_files(&filep,out_data_files,&filenames); 

#if !OUTPUT_FORCE

#if QUICK_FS
      if(options.FROZEN_SOIL) {
	for(i=0;i<MAX_LAYERS;i++) {
	  for(j=0;j<6;j++) 
	    free((char *)soil_con.ufwc_table_layer[i][j]);
	  free((char *)soil_con.ufwc_table_layer[i]);
	}
	for(i=0;i<MAX_NODES;i++) {
	  for(j=0;j<6;j++) 
	    free((char *)soil_con.ufwc_table_node[i][j]);
	  free((char *)soil_con.ufwc_table_node[i]);
	}
      }
#endif /* QUICK_FS */
      free_dist_prcp(&prcp,veg_con[0].vegetat_type_num);
      free_vegcon(&veg_con);
      free((char *)soil_con.AreaFract);
      free((char *)soil_con.Tfactor);
      free((char *)soil_con.Pfactor);
      free((char *)soil_con.AboveTreeLine);
      for(index=0;index<=options.Nlayer;index++) 
	free((char*)soil_con.layer_node_fract[index]);
      free((char*)soil_con.layer_node_fract);
#endif /* !OUTPUT_FORCE */
    }	/* End Run Model Condition */
  } 	/* End Grid Loop */

  /** cleanup **/
  free_atmos(global_param.nrecs, &atmos);
  free_dmy(&dmy);
  free_out_data_files(&out_data_files);
  free_out_data(&out_data);
#if !OUTPUT_FORCE
  free_veglib(&veg_lib);
#endif /* !OUTPUT_FORCE */
  fclose(filep.soilparam);
#if !OUTPUT_FORCE
  fclose(filep.vegparam);
  fclose(filep.veglib);
  if (options.SNOW_BAND>1)
    fclose(filep.snowband);
#endif /* !OUTPUT_FORCE */

  return EXIT_SUCCESS;
}	/* End Main Program */
