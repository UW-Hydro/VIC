#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
#include <string.h>
 
atmos_data_struct *read_forcing_data(infiles_struct      inf,
				     int                 starthour,
                                     int                *nrecs,
                                     int                 dt,
				     int                *force_dt)
/**********************************************************************
  read_forcing_data    Keith Cherkauer      April 6, 1998

  This subroutine was added to give users a place to add the routines
  to read the atmospheric forcing data used by the versions of their 
  models.

**********************************************************************/
{
  extern option_struct options;

  atmos_data_struct *temp;
  char              errorstr[MAXSTRING];

  temp = (atmos_data_struct*)calloc(*nrecs,sizeof(atmos_data_struct));

  /***** Read Forcing Files in PILPS2c Format *****/
  if(strcasecmp(options.FORCE_TYPE,"PILPS")==0) {
    read_PILPS2c(temp,inf.forcing[0],nrecs,dt,force_dt[0]);
  }
  /***** Read Forcing Files in SAWD Format *****/
  else if(strcasecmp(options.FORCE_TYPE,"SAWD")==0) {
    read_sawd(temp,inf.forcing[0],nrecs,dt,force_dt[0]);
    if(!options.HP) read_atmosdata(temp,inf.forcing[1],nrecs,dt,force_dt[1]);
  }
  /***** Read Forcing Files in SAWD Binary Format *****/
  else if(strcasecmp(options.FORCE_TYPE,"SAWD_BIN")==0) {
    read_sawd_binary(temp,inf.forcing[0],nrecs,dt,force_dt[0]);
    read_atmosdata(temp,inf.forcing[1],nrecs,dt,force_dt[1]);
  }
  /***** Read Daily Forcing Files *****/
  else if(strcasecmp(options.FORCE_TYPE,"DAILY")==0) {
    read_atmosdata(temp,inf.forcing[0],nrecs,dt,force_dt[0]);
    if(!options.SNOW_MODEL && inf.forcing[1] != NULL)
      read_snowmodel(temp,inf.forcing[1],nrecs[0],dt,force_dt[1]);
    else if(!options.SNOW_MODEL) {
      options.SNOW_MODEL = TRUE;
      fprintf(stderr,"WARNING: Check global parameter file, SNOW_MODEL ");
      fprintf(stderr,"not turned on, and no NWS snow model file defined\n");
      fprintf(stderr,"\tSetting SNOW_MODEL to TRUE to continue model run.\n");
    }
  }
  /***** Read SHAW Hourly Forcing Files *****/
  else if(strcasecmp(options.FORCE_TYPE,"SHAW")==0) {
    read_rosemount(temp,inf.forcing[0],nrecs,starthour,dt,force_dt[0]);
  }
  /***** Unrecognized Forcing File Type *****/
  else {
    sprintf(errorstr,"ERROR: Unrecognized forcing file type %s\n",
	    options.FORCE_TYPE);
    vicerror(errorstr);
  }

  return(temp);

}
