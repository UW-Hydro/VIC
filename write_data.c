#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

static char vcid[] = "$Id$";

void write_data(out_data_struct *out_data,
		outfiles_struct *outfiles,
		dmy_struct      *dmy,
		int              dt)
/**********************************************************************
	write_data	Dag Lohmann		Janurary 1996

  This subroutine writes all energy and moisture balance parameters to
  output files.

  OUTPUT:
	evaporation and vapor fluxes in mm/time step
	layer moisture in mm/time step
	runoff in mm/time step
	baseflow in mm/time step
	freezing and thawing depths in cm
	snow depth in cm
	snow water equivlence in mm
	all energy fluxes are in W/m^2

  Modifications:
  5/20/96	Program was modified to account for a variable
		number of soil layers.  It was also modified to
		write out frozen soils data per time step.	KAC
  1/15/97	Program modified to output daily sums, or values
		independant of selected time step.  This aids in
		comparisons between model versions.		KAC
  3/98          Routine modified to output fluxes in PILPS2c 
                ASCII column format                             Dag
  4/30/98       Routine modified to add binary output options for
                improved file speed, and less disk usage for large
		model basins                                    KAC
  7/19/99       modified to output a single binary file containing
                the data selected for the LDAS project         KAC
  8/3/99        modified again to reduce the storage space needed
                for the LDAS output files.  
  1/4/2000      modified to allow both standard and LDAS formatted
                output using a compiler flag                    KAC
  3/20/2001     made hour a variable in all output data file formats
                even if the model is run at a daily time step.  Also
                modified all output files to account for new 
                variables introduced by the spatial frost and snow
                algorithms, the lake algorithm and the PILPS 2e 
                study.                                          KAC

**********************************************************************/
{
  extern option_struct options;
#if LINK_DEBUG
  extern debug_struct debug;
#endif

  char               *tmp_cptr;
  short int          *tmp_siptr;
  unsigned short int *tmp_usiptr;
  int                 band, j;
  int                 fidx;
  int                *tmp_iptr;
  float              *tmp_fptr;

#if OPTIMIZE

  /*****************************************************************
    Create optimization output files
      type: ASCII
      columns: 5
               year
	       month
	       day
	       runoff
	       baseflow
      comment: runoff and baseflow are output as daily sums for all
               defined model time steps.
  ******************************************************************/

  fprintf(outfiles->fluxes,"%04i\t%02i\t%02i\t%.4f\t%.4f\n",
	  dmy->year, dmy->month, dmy->day, out_data->runoff,
	  out_data->baseflow);
  
#elif LDAS_OUTPUT 

  /*************************************************************
    Write output files using LDAS binary format
    - single file, limited output variables, output compressed

  LDAS output file format:
  unsigned short int             year
  char                           month
  char                           day
 (char)                         (hour)
  unsigned short int             prec * 100
  short int                      evap * 100
  float                          runoff
  float                          baseflow
  unsigned short int * Nlayers   moist * 10
  unsigned short int             swq * 100
  short int                      net_short * 10
  short int                      in_long * 10
  short int                      r_net * 10
  short int                      latent * 10
  short int                      sensible * 10
  short int                      grnd_flux * 10
  unsigned short int             albedo * 10000
  short int                      surf_temp * 100
  unsigned short int             rel_humid * 100
  short int                      air_temp*100
  unsigned short int             wind*100
  short int                      coverage*100
  if FROZEN_SOIL activated
    unsigned short int         ice[Nlayers] * 10     mean
    for each NUM_FRONTS:
      unsigned short int         fdepth[] * 100     mean
      unsigned short int         tdepth[] * 100     mean

  *************************************************************/

  /* initialize pointers */
  tmp_cptr = (char *)calloc(1,sizeof(char));
  tmp_siptr = (short int *)calloc(1,sizeof(short int));
  tmp_usiptr = (unsigned short int *)calloc(1,sizeof(unsigned short int));
  tmp_iptr = (int *)calloc(1,sizeof(int));
  tmp_fptr = (float *)calloc(1,sizeof(float));

  /***** Write Date and Time *****/
  tmp_usiptr[0] = (unsigned short int)dmy->year;
  fwrite(tmp_usiptr,1,sizeof(unsigned short int),outfiles->fluxes);
  tmp_cptr[0] = (char)dmy->month;
  fwrite(tmp_cptr,1,sizeof(char),outfiles->fluxes);
  tmp_cptr[0] = (char)dmy->day;
  fwrite(tmp_cptr,1,sizeof(char),outfiles->fluxes);
  if(dt<24) {
    tmp_cptr[0] = (char)dmy->hour;
    fwrite(tmp_cptr,1,sizeof(char),outfiles->fluxes);
  }
  
  /** water balance components **/

  /* precipitation */
  tmp_usiptr[0] = (unsigned short int)(out_data->prec * 100.);
  fwrite(tmp_usiptr,1,sizeof(unsigned short int),outfiles->fluxes);

  /* evaporation */
  tmp_siptr[0] = (short int)(out_data->evap * 100.);
  fwrite(tmp_siptr,1,sizeof(short int),outfiles->fluxes);
  
  /* runoff */
  tmp_fptr[0] = (float)out_data->runoff;
  fwrite(tmp_fptr,1,sizeof(float),outfiles->fluxes);

  /* baseflow */
  tmp_fptr[0] = (float)out_data->baseflow;
  fwrite(tmp_fptr,1,sizeof(float),outfiles->fluxes);

  /* liquid soil moisture */
  for(j=0;j<options.Nlayer;j++) {    
    tmp_usiptr[0] = (unsigned short int)(out_data->moist[j]*10.);
    fwrite(tmp_usiptr,1,sizeof(unsigned short int),outfiles->fluxes);
  }

  /* snow water equivalence */
  tmp_usiptr[0] = (unsigned short int)(out_data->swq[0]*100.);
  fwrite(tmp_usiptr,1,sizeof(unsigned short int),outfiles->fluxes);
  
  /** energy balance components **/

  /* net shortwave radiation */
  tmp_siptr[0] = (short int)(out_data->net_short*10.);
  fwrite(tmp_siptr,1,sizeof(short int),outfiles->fluxes);

  /* incoming longwave radiation */
  tmp_siptr[0] = (short int)(out_data->in_long*10.);
  fwrite(tmp_siptr,1,sizeof(short int),outfiles->fluxes);

  /* net radiation */
  tmp_siptr[0] = (short int)(out_data->r_net*10.);
  fwrite(tmp_siptr,1,sizeof(short int),outfiles->fluxes);

  /* latent heat flux */
  tmp_siptr[0] = (short int)((out_data->latent 
			      + out_data->latent_sub[0])*10.);
  fwrite(tmp_siptr,1,sizeof(short int),outfiles->fluxes);

  /* sensible heat flux */
  tmp_siptr[0] = (short int)((out_data->sensible 
			      + out_data->advected_sensible[0])*10.);
  fwrite(tmp_siptr,1,sizeof(short int),outfiles->fluxes);

  /* ground heat flux */
  tmp_siptr[0] = (short int)((out_data->grnd_flux + out_data->deltaH + out_data->fusion)*10.);
  fwrite(tmp_siptr,1,sizeof(short int),outfiles->fluxes);
  
  /* albedo */
  tmp_usiptr[0] = (unsigned short int)(out_data->albedo*10000.);
  fwrite(tmp_usiptr,1,sizeof(unsigned short int),outfiles->fluxes);

  /* raditaive temperature of surface */
  tmp_siptr[0] = (short int)((out_data->rad_temp-KELVIN)*100.);
  fwrite(tmp_siptr,1,sizeof(short int),outfiles->fluxes);

  /* relative humidity */
  tmp_usiptr[0] = (unsigned short int)(out_data->rel_humid*100.);
  fwrite(tmp_usiptr,1,sizeof(unsigned short int),outfiles->fluxes);

  /* air temperature */
  tmp_siptr[0] = (short int)(out_data->air_temp * 100.);
  fwrite(tmp_siptr,1,sizeof(short int),outfiles->fluxes);

  /* wind speed */
  tmp_siptr[0] = (short int)(out_data->wind * 100.);
  fwrite(tmp_siptr,1,sizeof(short int),outfiles->fluxes);

  /* snow cover fraction */
  tmp_siptr[0] = (short int)(out_data->coverage[0] * 100.);
  fwrite(tmp_siptr,1,sizeof(short int),outfiles->fluxes);

  /* frozen soil depths */
  if(options.FROZEN_SOIL) {

    /* total soil moisture content */
    for ( j = 0; j < options.Nlayer; j++ ) {    
      tmp_usiptr[0] = (unsigned short int)(out_data->ice[j]*10.);
      fwrite(tmp_usiptr,1,sizeof(unsigned short int),outfiles->fluxes);
    }
    for ( fidx = 0; fidx < MAX_FRONTS; fidx++ ) {
      /* frost front depths */
      tmp_usiptr[0] = (unsigned short int)(out_data->fdepth[fidx]*100.);
      fwrite(tmp_usiptr,1,sizeof(unsigned short int),outfiles->fluxes);
      /* thaw front depths */
      tmp_usiptr[0] = (unsigned short int)(out_data->tdepth[fidx]*100.);
      fwrite(tmp_usiptr,1,sizeof(unsigned short int),outfiles->fluxes);
    }
  }

  /* free pointers */
  free((char *)tmp_cptr);
  free((char *)tmp_siptr);
  free((char *)tmp_usiptr);
  free((char *)tmp_iptr);
  free((char *)tmp_fptr);

#else // !OPTIMIZE && !LDAS_OUTPUT

  /***************************************************************
    Write output files using default VIC ASCII or BINARY formats
    - multiple files, all variables, no truncation

    see VIC web page for format details:
      www.hydro.washington.edu/Lettenmaier/Models/VIC/VIChome.html
  ***************************************************************/

  /* initialize pointers */
  tmp_iptr = (int *)calloc(1,sizeof(int));
  tmp_fptr = (float *)calloc(1,sizeof(float));

  /************************************
    Output Frozen Soil Variables
  ************************************/
  if(options.FROZEN_SOIL && options.BINARY_OUTPUT) {
    /***** Write Binary Frozen Soil Output File *****/

    // year
    tmp_iptr[0] = dmy->year;
    fwrite(tmp_iptr,1,sizeof(int),outfiles->fdepth);

    // month
    tmp_iptr[0] = dmy->month;
    fwrite(tmp_iptr,1,sizeof(int),outfiles->fdepth);

    // day
    tmp_iptr[0] = dmy->day;
    fwrite(tmp_iptr,1,sizeof(int),outfiles->fdepth);

    // hour
    tmp_iptr[0] = dmy->hour;
    fwrite(tmp_iptr,1,sizeof(int),outfiles->fdepth);

    for ( fidx = 0; fidx < MAX_FRONTS; fidx++ ) {

      // frost depth [fidx]
      tmp_fptr[0] = (float)out_data->fdepth[fidx];
      fwrite(tmp_fptr,1,sizeof(float),outfiles->fdepth);

      // thaw depth [fidx]
      tmp_fptr[0] = (float)out_data->tdepth[fidx];
      fwrite(tmp_fptr,1,sizeof(float),outfiles->fdepth);
    }

    for ( j = 0; j < options.Nlayer; j++ ) {

      // total soil moisture
      tmp_fptr[0] = (float)(out_data->ice[j] + out_data->moist[j]);
      fwrite(tmp_fptr,1,sizeof(float),outfiles->fdepth);

    }

    // frozen fraction of ground surface
    tmp_fptr[0] = (float)(out_data->surf_frost_fract);
    fwrite(tmp_fptr,1,sizeof(float),outfiles->fdepth);
  }
  else if(options.FROZEN_SOIL) {
    /***** Write ASCII Frozen Soil Output File *****/
    fprintf(outfiles->fdepth  ,"%04i\t%02i\t%02i\t%02i",
	    dmy->year, dmy->month, dmy->day, dmy->hour);
    for ( fidx = 0; fidx < MAX_FRONTS; fidx++ ) {
      fprintf(outfiles->fdepth  ,"\t%.4f\t%.4f",out_data->fdepth[fidx], 
            out_data->tdepth[fidx]);
    }
    for ( j = 0; j < options.Nlayer; j++ ) {
      fprintf(outfiles->fdepth,"\t%f", out_data->ice[j] + out_data->moist[j]);
    }
    fprintf(outfiles->fdepth,"\t%f\n", out_data->surf_frost_fract); 
  }

#if LAKE_MODEL
  /************************
    Output Lake Variables
  ************************/
  if(options.LAKES && options.BINARY_OUTPUT) {
    /***** Write Binary Lake Output File *****/

    // year
    tmp_iptr[0] = dmy->year;
    fwrite(tmp_iptr,1,sizeof(int),outfiles->lake);

    // month
    tmp_iptr[0] = dmy->month;
    fwrite(tmp_iptr,1,sizeof(int),outfiles->lake);

    // day
    tmp_iptr[0] = dmy->day;
    fwrite(tmp_iptr,1,sizeof(int),outfiles->lake);

    // hour
    tmp_iptr[0] = dmy->hour;
    fwrite(tmp_iptr,1,sizeof(int),outfiles->lake);

    // lake ice temperature
    tmp_fptr[0] = (float)out_data->lake_ice_temp;
    fwrite(tmp_fptr,1,sizeof(float),outfiles->lake);

    // lake ice height
    tmp_fptr[0] = (float)out_data->lake_ice_height;
    fwrite(tmp_fptr,1,sizeof(float),outfiles->lake);

    // fractional lake ice coverage
    tmp_fptr[0] = (float)out_data->lake_ice_fract;
    fwrite(tmp_fptr,1,sizeof(float),outfiles->lake);

    // lake dpeth
    tmp_fptr[0] = (float)out_data->lake_depth;
    fwrite(tmp_fptr,1,sizeof(float),outfiles->lake);

    // lake surface area
    tmp_fptr[0] = (float)out_data->lake_surf_area;
    fwrite(tmp_fptr,1,sizeof(float),outfiles->lake);

    // lake volume
    tmp_fptr[0] = (float)out_data->lake_volume;
    fwrite(tmp_fptr,1,sizeof(float),outfiles->lake);

    // lake surface temperature
    tmp_fptr[0] = (float)out_data->lake_surf_temp;
    fwrite(tmp_fptr,1,sizeof(float),outfiles->lake);

  }
  else if(options.LAKES) {
    /***** Write ASCII Lake Output File *****/
    fprintf(outfiles->lake  ,"%04i\t%02i\t%02i\t%02i\t",
	    dmy->year,  dmy->month,  dmy->day,  dmy->hour);
    fprintf(outfiles->lake  ,"%04f\t%04f\t%04f\t%04f\t%04f\t%04f\t%04f\n",
	    out_data->lake_ice_temp, out_data->lake_ice_height, 
	    out_data->lake_ice_fract, out_data->lake_depth, 
	    out_data->lake_surf_area, out_data->lake_volume, 
	    out_data->lake_surf_temp);
  }
#endif // LAKE_MODEL

  /**************************************
    Ouput Snow Variables
  **************************************/
  if(options.BINARY_OUTPUT) {
    /***** Write Binary Snow Output File *****/

    // year
    tmp_iptr[0] = dmy->year;
    fwrite(tmp_iptr,1,sizeof(int),outfiles->snow);

    // month
    tmp_iptr[0] = dmy->month;
    fwrite(tmp_iptr,1,sizeof(int),outfiles->snow);

    // day
    tmp_iptr[0] = dmy->day;
    fwrite(tmp_iptr,1,sizeof(int),outfiles->snow);

    // hour
    tmp_iptr[0] = dmy->hour;
    fwrite(tmp_iptr,1,sizeof(int),outfiles->snow);

    // snow water equivalence
    tmp_fptr[0] = (float)out_data->swq[0];
    fwrite(tmp_fptr,1,sizeof(float),outfiles->snow);

    // snow depth
    tmp_fptr[0] = (float)out_data->snow_depth[0];
    fwrite(tmp_fptr,1,sizeof(float),outfiles->snow);

    // intercepted snow
    tmp_fptr[0] = (float)out_data->snow_canopy[0];
    fwrite(tmp_fptr,1,sizeof(float),outfiles->snow);

    // fractional snow coverage
    tmp_fptr[0] = (float)out_data->coverage[0];
    fwrite(tmp_fptr,1,sizeof(float),outfiles->snow);

    // advection
    tmp_fptr[0] = (float)out_data->advection[0];
    fwrite(tmp_fptr,1,sizeof(float),outfiles->snow);

    // change in snowpack cold content
    tmp_fptr[0] = (float)out_data->deltaCC[0];
    fwrite(tmp_fptr,1,sizeof(float),outfiles->snow);

    // snow heat flux
    tmp_fptr[0] = (float)out_data->snow_flux[0];
    fwrite(tmp_fptr,1,sizeof(float),outfiles->snow);

    // snowpack refreeze energy
    tmp_fptr[0] = (float)out_data->refreeze_energy[0];
    fwrite(tmp_fptr,1,sizeof(float),outfiles->snow);

    // snowpack melt energy
    tmp_fptr[0] = (float)out_data->melt_energy[0];
    fwrite(tmp_fptr,1,sizeof(float),outfiles->snow);

    // advected sensible heat
    tmp_fptr[0] = (float)out_data->advected_sensible[0];
    fwrite(tmp_fptr,1,sizeof(float),outfiles->snow);

    // latent heat of sublimation
    tmp_fptr[0] = (float)out_data->latent_sub[0];
    fwrite(tmp_fptr,1,sizeof(float),outfiles->snow);

    // snow surface temperature
    tmp_fptr[0] = (float)out_data->snow_surf_temp[0];
    fwrite(tmp_fptr,1,sizeof(float),outfiles->snow);

    // snow pack layer temperature
    tmp_fptr[0] = (float)out_data->snow_pack_temp[0];
    fwrite(tmp_fptr,1,sizeof(float),outfiles->snow);

    // snow melt
    tmp_fptr[0] = (float)out_data->melt[0];
    fwrite(tmp_fptr,1,sizeof(float),outfiles->snow);

  }
  else {
    /***** Write ASCII full energy snow output file *****/
    fprintf(outfiles->snow, "%04i\t%02i\t%02i\t", dmy->year, dmy->month, 
	    dmy->day);
    fprintf(outfiles->snow, "%02i\t", dmy->hour);
    fprintf(outfiles->snow, "%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n",
	    out_data->swq[0], out_data->snow_depth[0], 
	    out_data->snow_canopy[0], out_data->coverage[0], 
	    out_data->advection[0], out_data->deltaCC[0], 
	    out_data->snow_flux[0], out_data->refreeze_energy[0], 
	    out_data->melt_energy[0], out_data->advected_sensible[0], 
	    out_data->latent_sub[0], out_data->snow_surf_temp[0], 
	    out_data->snow_pack_temp[0], out_data->melt[0] );
  }

  /**************************************
    Ouput Snow Band Variables
  **************************************/
  if(options.PRT_SNOW_BAND && options.BINARY_OUTPUT) {
    /***** Write Binary Snow Band Output File *****/

    // year
    tmp_iptr[0] = dmy->year;
    fwrite(tmp_iptr,1,sizeof(int),outfiles->snowband);

    // month
    tmp_iptr[0] = dmy->month;
    fwrite(tmp_iptr,1,sizeof(int),outfiles->snowband);

    // day
    tmp_iptr[0] = dmy->day;
    fwrite(tmp_iptr,1,sizeof(int),outfiles->snowband);

    // hour
    tmp_iptr[0] = dmy->hour;
    fwrite(tmp_iptr,1,sizeof(int),outfiles->snowband);

    for ( band = 0; band < options.SNOW_BAND; band++ ) {

      // snow water equivalence
      tmp_fptr[0] = (float)out_data->swq[band+1];
      fwrite(tmp_fptr,1,sizeof(float),outfiles->snowband);

      // snow depth
      tmp_fptr[0] = (float)out_data->snow_depth[band+1];
      fwrite(tmp_fptr,1,sizeof(float),outfiles->snowband);

      // intercepted snow
      tmp_fptr[0] = (float)out_data->snow_canopy[band+1];
      fwrite(tmp_fptr,1,sizeof(float),outfiles->snowband);

      // fractional snow coverage
      tmp_fptr[0] = (float)out_data->coverage[0];
      fwrite(tmp_fptr,1,sizeof(float),outfiles->snowband);
      
      // advection
      tmp_fptr[0] = (float)out_data->advection[band+1];
      fwrite(tmp_fptr,1,sizeof(float),outfiles->snowband);

      // change in snowpack cold content
      tmp_fptr[0] = (float)out_data->deltaCC[band+1];
      fwrite(tmp_fptr,1,sizeof(float),outfiles->snowband);

      // snow heat flux
      tmp_fptr[0] = (float)out_data->snow_flux[band+1];
      fwrite(tmp_fptr,1,sizeof(float),outfiles->snowband);

      // snowpack refreeze energy
      tmp_fptr[0] = (float)out_data->refreeze_energy[band+1];
      fwrite(tmp_fptr,1,sizeof(float),outfiles->snowband);

      // snowpack melt energy
      tmp_fptr[0] = (float)out_data->melt_energy[band+1];
      fwrite(tmp_fptr,1,sizeof(float),outfiles->snow);

      // advected sensible heat
      tmp_fptr[0] = (float)out_data->advected_sensible[band+1];
      fwrite(tmp_fptr,1,sizeof(float),outfiles->snow);

      // latent heat of sublimation
      tmp_fptr[0] = (float)out_data->latent_sub[band+1];
      fwrite(tmp_fptr,1,sizeof(float),outfiles->snow);

      // snowpack surface layer temperature
      tmp_fptr[0] = (float)out_data->snow_surf_temp[band+1];
      fwrite(tmp_fptr,1,sizeof(float),outfiles->snow);

      // snowpack pack layer temperature
      tmp_fptr[0] = (float)out_data->snow_pack_temp[band+1];
      fwrite(tmp_fptr,1,sizeof(float),outfiles->snow);

      // snowpack melt
      tmp_fptr[0] = (float)out_data->melt[band+1];
      fwrite(tmp_fptr,1,sizeof(float),outfiles->snow);

    }
  }
  else if(options.PRT_SNOW_BAND) {
    /***** Write ASCII snow band output file *****/
    fprintf(outfiles->snowband, "%04i\t%02i\t%02i\t", dmy->year, 
	    dmy->month, dmy->day);
    fprintf(outfiles->snowband, "%02i\t", dmy->hour);
    for ( band = 0; band < options.SNOW_BAND; band++) 
      fprintf(outfiles->snowband  ,"\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f",
	      out_data->swq[band+1], out_data->snow_depth[band+1], 
	      out_data->snow_canopy[band+1], out_data->coverage[band+1], 
	      out_data->advection[band+1], out_data->deltaCC[band+1], 
	      out_data->snow_flux[band+1], out_data->refreeze_energy[band+1], 
	      out_data->melt_energy[band+1], 
	      out_data->advected_sensible[band+1], 
	      out_data->latent_sub[band+1], out_data->snow_surf_temp[band+1], 
	      out_data->snow_pack_temp[band+1], out_data->melt[band+1] );
    fprintf(outfiles->snowband,"\n");
  }

  /************************************
    Output Standard Energy and Moisture Flux Variables
  ************************************/
  if(options.BINARY_OUTPUT) {
    /***** Write Binary Fluxes File *****/

    // year
    tmp_iptr[0] = dmy->year;
    fwrite(tmp_iptr,1,sizeof(int),outfiles->fluxes);

    // month
    tmp_iptr[0] = dmy->month;
    fwrite(tmp_iptr,1,sizeof(int),outfiles->fluxes);

    // day
    tmp_iptr[0] = dmy->day;
    fwrite(tmp_iptr,1,sizeof(int),outfiles->fluxes);

    // hour
    tmp_iptr[0] = dmy->hour;
    fwrite(tmp_iptr,1,sizeof(int),outfiles->fluxes);
    
    // precipitation
    tmp_fptr[0] = (float)out_data->prec;
    fwrite(tmp_fptr,1,sizeof(float),outfiles->fluxes);
    
    // evaporation
    tmp_fptr[0] = (float)out_data->evap;
    fwrite(tmp_fptr,1,sizeof(float),outfiles->fluxes);

    // runoff
    tmp_fptr[0] = (float)out_data->runoff;
    fwrite(tmp_fptr,1,sizeof(float),outfiles->fluxes);
    
    // baseflow
    tmp_fptr[0] = (float)out_data->baseflow;
    fwrite(tmp_fptr,1,sizeof(float),outfiles->fluxes);

    // Wdew
    tmp_fptr[0] = (float)out_data->Wdew;
    fwrite(tmp_fptr,1,sizeof(float),outfiles->fluxes);

    for ( j = 0; j < options.Nlayer; j++ ) {    

      // unfrozen soil moisture
      tmp_fptr[0] = (float)out_data->moist[j];
      fwrite(tmp_fptr,1,sizeof(float),outfiles->fluxes);

    }

    // radiative temperature
    tmp_fptr[0] = (float)out_data->rad_temp;
    fwrite(tmp_fptr,1,sizeof(float),outfiles->fluxes);

    // net shortwave radiation
    tmp_fptr[0] = (float)out_data->net_short;
    fwrite(tmp_fptr,1,sizeof(float),outfiles->fluxes);

    // net radiation
    tmp_fptr[0] = (float)out_data->r_net;
    fwrite(tmp_fptr,1,sizeof(float),outfiles->fluxes);

    // latent heat
    tmp_fptr[0] = (float)out_data->latent;
    fwrite(tmp_fptr,1,sizeof(float),outfiles->fluxes);

    // canopy evaporation
    tmp_fptr[0] = (float)out_data->evap_canop;
    fwrite(tmp_fptr,1,sizeof(float),outfiles->fluxes);

    // vegetation transpiration
    tmp_fptr[0] = (float)out_data->evap_veg;
    fwrite(tmp_fptr,1,sizeof(float),outfiles->fluxes);

    // bare soil evaporation
    tmp_fptr[0] = (float)out_data->evap_bare;
    fwrite(tmp_fptr,1,sizeof(float),outfiles->fluxes);

    // canopy sublimation
    tmp_fptr[0] = (float)out_data->sub_canop;
    fwrite(tmp_fptr,1,sizeof(float),outfiles->fluxes);

    // snowpack sublimation
    tmp_fptr[0] = (float)out_data->sub_snow;
    fwrite(tmp_fptr,1,sizeof(float),outfiles->fluxes);

    // sensible heat
    tmp_fptr[0] = (float)out_data->sensible;
    fwrite(tmp_fptr,1,sizeof(float),outfiles->fluxes);

    // ground heat flux
    tmp_fptr[0] = (float)out_data->grnd_flux;
    fwrite(tmp_fptr,1,sizeof(float),outfiles->fluxes);

    // heat storage
    tmp_fptr[0] = (float)out_data->deltaH;
    fwrite(tmp_fptr,1,sizeof(float),outfiles->fluxes);

    // fusion
    tmp_fptr[0] = (float)out_data->fusion;
    fwrite(tmp_fptr,1,sizeof(float),outfiles->fluxes);

    // aerodynamic resistance
    tmp_fptr[0] = (float)out_data->aero_resist;
    fwrite(tmp_fptr,1,sizeof(float),outfiles->fluxes);

    // surface temperature
    tmp_fptr[0] = (float)out_data->surf_temp;
    fwrite(tmp_fptr,1,sizeof(float),outfiles->fluxes);

    // albedo
    tmp_fptr[0] = (float)out_data->albedo;
    fwrite(tmp_fptr,1,sizeof(float),outfiles->fluxes);

    // relative humidity
    tmp_fptr[0] = (float)out_data->rel_humid;
    fwrite(tmp_fptr,1,sizeof(float),outfiles->fluxes);

    // incoming longwave
    tmp_fptr[0] = (float)out_data->in_long;
    fwrite(tmp_fptr,1,sizeof(float),outfiles->fluxes);

    // air temperature
    tmp_fptr[0] = (float)out_data->air_temp;
    fwrite(tmp_fptr,1,sizeof(float),outfiles->fluxes);

    // wind speed
    tmp_fptr[0] = (float)out_data->wind;
    fwrite(tmp_fptr,1,sizeof(float),outfiles->fluxes);

  }
  else {
    /***** Write ASCII fluxes file *****/
    fprintf(outfiles->fluxes,"%04i\t%02i\t%02i\t", dmy->year, dmy->month, 
	    dmy->day);
    fprintf(outfiles->fluxes,"%02i\t", dmy->hour); 
    fprintf(outfiles->fluxes,"%.4f\t%.4f\t%.4f\t%.4f\t%.4f",
	    out_data->prec, out_data->evap, out_data->runoff,
	    out_data->baseflow, out_data->Wdew);
    for ( j = 0; j < options.Nlayer; j++ ) 
      fprintf(outfiles->fluxes,"\t%.4f", out_data->moist[j]);
    fprintf(outfiles->fluxes, "\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f",
	    out_data->rad_temp, out_data->net_short, out_data->r_net, 
	    out_data->latent, out_data->evap_canop, out_data->evap_veg );
    fprintf(outfiles->fluxes, "\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f",
	    out_data->evap_bare, out_data->sub_canop, out_data->sub_snow, 
	    out_data->sensible, out_data->grnd_flux, out_data->deltaH );
    fprintf(outfiles->fluxes, "\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f",
	    out_data->fusion, out_data->aero_resist, out_data->surf_temp, 
	    out_data->albedo, out_data->rel_humid, out_data->in_long );
    fprintf(outfiles->fluxes, "\t%.4f\t%.4f\n",
	    out_data->air_temp, out_data->wind );
  }

  free((char *)tmp_iptr);
  free((char *)tmp_fptr);

#endif

}

void calc_water_balance_error(int    rec,
			      double inflow,
			      double outflow,
			      double storage) {
/***************************************************************
  calc_water_balance_error  Keith Cherkauer        April 1998

  This subroutine computes the overall model water balance, and 
  warns the model user if large errors are found.
***************************************************************/

  static double last_storage;
  static double cum_error;
  static double max_error;
  static int    error_cnt;
  static int    Nrecs;

  double error;

  if(rec<0) {
    last_storage = storage;
    cum_error    = 0.;
    max_error    = 0.;
    error_cnt    = 0;
    Nrecs        = -rec;
  }
  else {
    error = inflow - outflow - (storage - last_storage);
    cum_error += error;
    if(fabs(error)>fabs(max_error) && fabs(error)>1e-5) {
      max_error = error;
      fprintf(stderr,"Maximum Moist Error:\t%i\t%.5f\t%.5f\n",
	      rec,error,cum_error);
    }
    if(rec==Nrecs-1) {
      fprintf(stderr,"Total Cumulative Water Error for Grid Cell = %.4f\n",
	      cum_error);
    }
    last_storage = storage;
  }

}

void calc_energy_balance_error(int    rec,
			       double net_rad,
			       double latent,
			       double sensible,
			       double grnd_flux,
			       double snow_fluxes) {
/***************************************************************
  calc_energy_balance_error   Keith Cherkauer     April 1998

  This subroutine computes the overall model energy balance, and
  reports the maximum time step error above a thresehold to the
  user.  The total cumulative error for the grid cell is also 
  computed and reported at the end of the model run.
***************************************************************/

  static double cum_error;
  static double max_error;
  static int    Nrecs;

  double error;

  if(rec<0) {
    cum_error = 0;
    Nrecs     = -rec;
    max_error = 0;
  }
  else {
    error = net_rad - latent - sensible - grnd_flux + snow_fluxes;
    cum_error += error;
    if(fabs(error)>fabs(max_error) && fabs(error)>0.001) {
      max_error = error;
      if ( rec > 0 ) 
	fprintf(stderr,"Maximum Energy Error:\t%i\t%.4f\t%.4f\n",
		rec,error,cum_error/(double)rec);
      else 
	fprintf(stderr,"Maximum Energy Error:\t%i\t%.4f\t%.4f\n",
		rec,error,cum_error);
    }
    if(rec==Nrecs-1) {
      fprintf(stderr,"Total Cumulative Energy Error for Grid Cell = %.4f\n",
	      cum_error/(double)rec);
    }
  }
}

