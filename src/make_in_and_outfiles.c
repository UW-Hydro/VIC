#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vicNl.h>

static char vcid[] = "$Id$";

void make_in_and_outfiles(filep_struct         *filep, 
			  filenames_struct     *filenames,
			  soil_con_struct      *soil,
			  out_data_file_struct *out_data_files)
/**********************************************************************
	make_in_and_outfile	Dag Lohman	January 1996

  This program builds the files names for input and output of grided
  data files.

  Modifications:
  5/20/96	The routine was modified to accept a variable
		number of layers, as well as to work with 
		frozen soils					KAC
  11-18-02 Modified to print notification that the output fluxes file
           will be in a binary format.                          LCB
  29-Oct-03 Distinguishing between input lakeparam file and output
	    lake file.							TJB
  2005-Mar-24 Modified to handle ALMA output files.			TJB
  2006-Sep-23 Implemented flexible output configuration; uses the new
              out_data and out_data_files structures; removed the
              OPTIMIZE and LDAS_OUTPUT options.				TJB
  2006-Oct-16 Merged infiles and outfiles structs into filep_struct;
	      Merged builtnames into filenames->			TJB
  2007-Oct-31 Append "/" to result_dir so that this need not be done
	      in global parameter file.					TJB
  2011-May-25 Expanded latchar, lngchar, and junk allocations to handle
	      GRID_DECIMAL > 4.						TJB

**********************************************************************/
{
  extern option_struct    options;
  extern param_set_struct param_set;
  extern FILE *open_file(char string[], char type[]);

  char   latchar[20], lngchar[20], junk[6];
  int filenum;

  sprintf(junk, "%%.%if", options.GRID_DECIMAL);
  sprintf(latchar, junk, soil->lat);
  sprintf(lngchar, junk, soil->lng);
 
  /********************************
  Input Forcing Files
  ********************************/

  strcpy(filenames->forcing[0], filenames->f_path_pfx[0]);
  strcat(filenames->forcing[0], latchar);
  strcat(filenames->forcing[0], "_");
  strcat(filenames->forcing[0], lngchar);
  if(param_set.FORCE_FORMAT[0] == BINARY)
    filep->forcing[0] = open_file(filenames->forcing[0], "rb");
  else
    filep->forcing[0] = open_file(filenames->forcing[0], "r");

  filep->forcing[1] = NULL;
  if(strcasecmp(filenames->f_path_pfx[1],"MISSING")!=0) {
    strcpy(filenames->forcing[1], filenames->f_path_pfx[1]);
    strcat(filenames->forcing[1], latchar);
    strcat(filenames->forcing[1], "_");
    strcat(filenames->forcing[1], lngchar);
    if(param_set.FORCE_FORMAT[0] == BINARY) 
      filep->forcing[1] = open_file(filenames->forcing[1], "rb");
    else 
      filep->forcing[1] = open_file(filenames->forcing[1], "r");
  }
  filep->forcing[2] = NULL;
  if(strcasecmp(filenames->f_path_pfx[2],"MISSING")!=0) {
    strcpy(filenames->forcing[2], filenames->f_path_pfx[2]);
    strcat(filenames->forcing[2], latchar);
    strcat(filenames->forcing[2], "_");
    strcat(filenames->forcing[2], lngchar);
    if(param_set.FORCE_FORMAT[0] == BINARY) 
      filep->forcing[2] = open_file(filenames->forcing[2], "rb");
    else 
      filep->forcing[2] = open_file(filenames->forcing[2], "r");
  }

  /********************************
  Output Files
  ********************************/

  for (filenum=0; filenum<options.Noutfiles; filenum++) {
    strcpy(out_data_files[filenum].filename, filenames->result_dir);
    strcat(out_data_files[filenum].filename, "/");
    strcat(out_data_files[filenum].filename, out_data_files[filenum].prefix);
    strcat(out_data_files[filenum].filename, "_");
    strcat(out_data_files[filenum].filename, latchar);
    strcat(out_data_files[filenum].filename, "_");
    strcat(out_data_files[filenum].filename, lngchar);
    if(options.BINARY_OUTPUT)
      out_data_files[filenum].fh = open_file(out_data_files[filenum].filename, "wb");
    else out_data_files[filenum].fh = open_file(out_data_files[filenum].filename, "w");
  }

} 
