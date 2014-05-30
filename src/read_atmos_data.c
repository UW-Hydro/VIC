#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vicNl.h>

static char vcid[] = "$Id$";

void read_atmos_data(FILE                 *infile,
		     global_param_struct   global_param,
		     int                   file_num,
		     int                   forceskip,
		     double              **forcing_data,
		     double             ***veg_hist_data)
/**********************************************************************
  read_atmos_data
  
  This routine reads in atmospheric data values from a binary/ascii file.

  BINARY
  Binary data are always specified as unsigned or signed ints, and a
  multiplier is used to convert to float.

  atmos variable: type:            model units:
  
  precipitation   unsigned short   mm per file_dt
  temperature       signed short   C
  wind              signed short   m/s

  swap bytes code from Kernighan, Brian W. and Rob Pike, "The practice of
  programming", Addison-Wesley, Reading, Massachusetts, 1999, 267 pp,
  page 206.   		
  
  ASCII
  ASCII data should have the same units as given in the table above.

  
  Supported Input Field Combinations, options in parenthesis optional:
  
  Daily met data and daily model timestep - prec, tmax, tmin, (wind)
  Daily met data and subdaily model timestep - prec, tmax, tmin, (wind)
  
  If the code is modified check;
  * for BINARY, the number of fields is correct
  * get_global flags are implemented
  
  Modifications:
  01/10/00 Modified to read a generic Binary, or ASCII column 
           data file and read the contents into the provided  
           data arrays.                                    KAC
  ??-???-?? Replaced NF with global_param.dt in condition checking
	    whether forcing file contains enough records to cover
	    the time range of the simulation.	(documented by TJB)
  2006-08-23 Changed order of fread/fwrite statements from ...1, sizeof...
             to ...sizeof, 1,... 					GCT
  2006-Sep-23 Fixed RCS ID string.					TJB
  2007-Jan-15 Added PRT_HEADER option; now binary forcing files
	      might have headers, and these need to be skipped.		TJB
  2011-Nov-04 Fixed warning message dealing with insufficient
	      records.							TJB
  2014-Apr-25 Added non-climatological veg parameters (as forcing
	      variables).						TJB
  2014-Apr-25 Added partial vegcover fraction.				TJB

  **********************************************************************/
{
  
  extern option_struct options;
  extern param_set_struct param_set;
  
  int             rec;
  int             skip_recs;
  int             i,j;
  int             endian;
  int             fields;
  int             Nfields;
  int             day=0;
  int            *field_index;
  unsigned short  ustmp;
  signed short    stmp;
  char            str[MAXSTRING+1];
  char            ErrStr[MAXSTRING+1];
  unsigned short  Identifier[4];
  int             Nbytes;

  Nfields     = param_set.N_TYPES[file_num];
  field_index = param_set.FORCE_INDEX[file_num];

  /** locate starting record **/
  /* if ascii then the following refers to the number of lines to skip,
     if binary the following needs multiplying by the number of input fields */
  skip_recs = (int)((float)(global_param.dt * forceskip)) 
    / (float)param_set.FORCE_DT[file_num];
  if((((global_param.dt < 24 && (param_set.FORCE_DT[file_num] * forceskip) 
	% global_param.dt) > 0)) 
     || (global_param.dt == 24 && (global_param.dt 
				   % param_set.FORCE_DT[file_num] > 0)))
    nrerror("Currently unable to handle a model starting date that does not correspond to a line in the forcing file.");

  /** Error checking - Model can be run at any time step using daily forcing
      data, but if sub-daily data is used, the model must be run at the
      same time step as the data.  That way aggregation and disaggragation 
      techniques are left to the user. **/
  if(param_set.FORCE_DT[file_num] < 24 
     && global_param.dt != param_set.FORCE_DT[file_num]) {
    sprintf(ErrStr,"When forcing the model with sub-daily data, the model must be run at the same time step as the forcing data.  Currently the model time step is %i hours, while forcing file %i has a time step of %i hours.",global_param.dt,file_num,param_set.FORCE_DT[file_num]);
    nrerror(ErrStr);
  }

  if(infile==NULL)fprintf(stderr,"NULL file\n");
  
  /***************************
    Read BINARY Forcing Data
  ***************************/

  if(param_set.FORCE_FORMAT[file_num] == BINARY){
	  
    /** test whether the machine is little-endian or big-endian **/
    i = 1;
    if(*(char *)&i == 1)
      endian = LITTLE;
    else    
      endian = BIG;
	  
    // Check for presence of a header, & skip over it if appropriate.
    // A VIC header will start with 4 instances of the identifier,
    // followed by number of bytes in the header (Nbytes).
    // Nbytes is assumed to be the byte offset at which the data records start.
    fseek(infile,0,SEEK_SET);
    if (feof(infile))
      nrerror("No data in the forcing file.  Model stopping...");
    for (i=0; i<4; i++) {
      fread(&ustmp,sizeof(unsigned short),1,infile);
      if (endian != param_set.FORCE_ENDIAN[file_num]) {
        ustmp = ((ustmp & 0xFF) << 8) | ((ustmp >> 8) & 0xFF);
      }
      Identifier[i] = ustmp;
    }
    if (Identifier[0] != 0xFFFF || Identifier[1] != 0xFFFF || Identifier[2] != 0xFFFF || Identifier[3] != 0xFFFF) {
      Nbytes = 0;
    }
    else {
      fread(&ustmp,sizeof(unsigned short),1,infile);
      if (endian != param_set.FORCE_ENDIAN[file_num]) {
        ustmp = ((ustmp & 0xFF) << 8) | ((ustmp >> 8) & 0xFF);
      }
      Nbytes = (int)ustmp;
    }
    fseek(infile,Nbytes,SEEK_SET);


    /** if forcing file starts before the model simulation, 
	skip over its starting records **/
    fseek(infile,skip_recs*Nfields*sizeof(short),SEEK_CUR);
    if (feof(infile))
      nrerror("No data for the specified time period in the forcing file.  Model stopping...");
	  
    /** Read BINARY forcing data **/
    rec = 0;
	  
    while ( !feof(infile) && (rec * param_set.FORCE_DT[file_num] 
			      < global_param.nrecs * global_param.dt) ) {

      for(i=0;i<Nfields;i++) {
        if (field_index[i] != ALBEDO && field_index[i] != LAI_IN && field_index[i] != VEGCOVER) {
	  if(param_set.TYPE[field_index[i]].SIGNED) {
	    fread(&stmp,sizeof(short int),1,infile);
	    if (endian != param_set.FORCE_ENDIAN[file_num]) {
	      stmp = ((stmp & 0xFF) << 8) | ((stmp >> 8) & 0xFF);
	    }
	    forcing_data[field_index[i]][rec] 
	      = (double)stmp / param_set.TYPE[field_index[i]].multiplier;
	  }
	  else {
	    fread(&ustmp,sizeof(unsigned short int),1,infile);
	    if (endian != param_set.FORCE_ENDIAN[file_num]) {
	      ustmp = ((ustmp & 0xFF) << 8) | ((ustmp >> 8) & 0xFF);
	    }
	    forcing_data[field_index[i]][rec] 
	      = (double)ustmp / param_set.TYPE[field_index[i]].multiplier;
	  }
	}
	else {
          for(j=0;j<param_set.TYPE[field_index[i]].N_ELEM;j++) {
	    if(param_set.TYPE[field_index[i]].SIGNED) {
	      fread(&stmp,sizeof(short int),1,infile);
	      if (endian != param_set.FORCE_ENDIAN[file_num]) {
	        stmp = ((stmp & 0xFF) << 8) | ((stmp >> 8) & 0xFF);
	      }
	      veg_hist_data[field_index[i]][j][rec] 
	        = (double)stmp / param_set.TYPE[field_index[i]].multiplier;
	    }
	    else {
	      fread(&ustmp,sizeof(unsigned short int),1,infile);
	      if (endian != param_set.FORCE_ENDIAN[file_num]) {
	        ustmp = ((ustmp & 0xFF) << 8) | ((ustmp >> 8) & 0xFF);
	      }
	      veg_hist_data[field_index[i]][j][rec] 
	        = (double)ustmp / param_set.TYPE[field_index[i]].multiplier;
	    }
          }
	}
      }
			
      rec++;
			
    }
  }

  /**************************
    Read ASCII Forcing Data 
  **************************/

  else{
	  
    // No need to skip over a header here, since ascii file headers are skipped
    // in open_file().  However, if we wanted to read information from the header,
    // we'd want to do it here, after rewinding to the beginning of the file (or
    // moving the code that deals with headers from open_file() to this function
    // and to any other functions that read the files, so that those functions could
    // also read the headers if necessary).

    /* skip to the beginning of the required met data */
    for(i=0;i<skip_recs;i++){
      if( fgets(str, MAXSTRING, infile) == NULL )
	nrerror("No data for the specified time period in the forcing file.  Model stopping...");
    }
	  
    /* read forcing data */
    rec=0;

    while( !feof(infile) && (rec * param_set.FORCE_DT[file_num] 
			      < global_param.nrecs * global_param.dt ) ) {
      for(i=0;i<Nfields;i++) {
        if (field_index[i] != ALBEDO && field_index[i] != LAI_IN && field_index[i] != VEGCOVER) {
	  fscanf(infile,"%lf", &forcing_data[field_index[i]][rec]);
        }
        else {
          for(j=0;j<param_set.TYPE[field_index[i]].N_ELEM;j++) {
	    fscanf(infile,"%lf", &veg_hist_data[field_index[i]][j][rec]);
          }
        }
      }
      fgets(str, MAXSTRING, infile);
      rec++;
    }
  }
  
  if(rec * param_set.FORCE_DT[file_num] 
     < global_param.nrecs * global_param.dt ) {
    sprintf(ErrStr,"Not enough records in forcing file %i (%i * %i = %i) to run the number of records defined in the global file (%i * %i = %i).  Check forcing file time step, and global file", file_num+1, rec, param_set.FORCE_DT[file_num],
	    rec*param_set.FORCE_DT[file_num], global_param.nrecs, 
	    global_param.dt, global_param.nrecs*global_param.dt);
    nrerror(ErrStr);
  }
  
}
