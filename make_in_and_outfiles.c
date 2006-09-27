#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vicNl.h>

static char vcid[] = "$Id$";

filenames_struct make_in_and_outfiles(infiles_struct       *infp, 
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
	    lake file.						TJB
  2005-Mar-24 Modified to handle ALMA output files.		TJB
  2006-Sep-23 Implemented flexible output configuration; uses the new
              out_data and out_data_files structures; removed the
              OPTIMIZE and LDAS_OUTPUT options. TJB

**********************************************************************/
{
  extern option_struct    options;
  extern param_set_struct param_set;
  extern FILE *open_file(char string[], char type[]);

  char             latchar[10], lngchar[10], junk[5];
  filenames_struct fnames;
  int filenum;

  fnames = *filenames;

  sprintf(junk, "%%.%if", options.GRID_DECIMAL);
  sprintf(latchar, junk, soil->lat);
  sprintf(lngchar, junk, soil->lng);
 
  /********************************
  Input Forcing Files
  ********************************/

  strcat(fnames.forcing[0], latchar);
  strcat(fnames.forcing[0], "_");
  strcat(fnames.forcing[0], lngchar);
  if(param_set.FORCE_FORMAT[0] == BINARY)
    infp->forcing[0] = open_file(fnames.forcing[0], "rb");
  else
    infp->forcing[0] = open_file(fnames.forcing[0], "r");

  infp->forcing[1] = NULL;
  if(strcasecmp(fnames.forcing[1],"FALSE")!=0) {
    strcat(fnames.forcing[1], latchar);
    strcat(fnames.forcing[1], "_");
    strcat(fnames.forcing[1], lngchar);
    if(param_set.FORCE_FORMAT[0] == BINARY) 
      infp->forcing[1] = open_file(fnames.forcing[1], "rb");
    else 
      infp->forcing[1] = open_file(fnames.forcing[1], "r");
  }

  /********************************
  Output Files
  ********************************/

  for (filenum=0; filenum<options.Noutfiles; filenum++) {
    strcpy(out_data_files[filenum].filename, fnames.result_dir);
    strcat(out_data_files[filenum].filename, out_data_files[filenum].prefix);
    strcat(out_data_files[filenum].filename, "_");
    strcat(out_data_files[filenum].filename, latchar);
    strcat(out_data_files[filenum].filename, "_");
    strcat(out_data_files[filenum].filename, lngchar);
    if(options.BINARY_OUTPUT)
      out_data_files[filenum].fh = open_file(out_data_files[filenum].filename, "wb");
    else out_data_files[filenum].fh = open_file(out_data_files[filenum].filename, "w");
  }

  return (fnames);

} 
