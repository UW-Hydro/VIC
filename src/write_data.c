#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

static char vcid[] = "$Id$";

void write_data(out_data_file_struct *out_data_files,
		out_data_struct *out_data,
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
  3-12-03       added energy fluxes to snow band output files   KAC
  04-22-03      Updated output of model for lakes and wetlands algorithm.
                Added output of blowing snow sublimation to LDAS and
                standard snow output files.  ** No Lake Variables are
                included in the LDAS output format. **         KAC
  04-23-2003    modified LDAS SWQ output, so that it is multiplied by
                10 instead of 100 before being converted to a short
                integer.  This reduces stored value precision to 0.1,
                but increases the maximum storable SWQ, which was
                exceeded in previous LDAS simulations.          KAC
  2003-Jul-07 Corrected output of sub_snow variable to item [0]
              rather than a point - will need to decide what
              parts of this array are important to output.		KAC
  2003-Oct-30 Replaced output of sub_snow[0] in fluxes file with
              sub_total.						TJB
  2005-Mar-24 Added support for ALMA variables.				TJB
  2005-Nov-08 Corrected outfiles from fluxes to snow for blowing snow
              sublimation. Corrected outfiles from snow to snowband
              for options.PRT_SNOW_BAND. Removed the following from
              snowband output: net sw radiation, net lw, albedo,
              latent heat flux, sensible heat flux, ground heat flux.	GCT
  2006-Aug-23 Changed order of fread/fwrite statements from ...1, sizeof...
              to ...sizeof, 1,...					GCT
  2006-Sep-23 Implemented flexible output configuration; uses the new
              out_data and out_data_files structures; moved the functions
              calc_energy_balance_error and calc_water_balance_error to
              the file calc_water_energy_balance_errors.c; implemented
	      aggregation of output variables.				TJB
  2012-Jan-16 Removed LINK_DEBUG code					BN
  2013-Dec-27 Moved OUTPUT_FORCE to options_struct.			TJB
**********************************************************************/
{
  extern option_struct options;
  int                 file_idx;
  int                 var_idx;
  int                 elem_idx;
  int                 ptr_idx;
  char               *tmp_cptr;
  short int          *tmp_siptr;
  unsigned short int *tmp_usiptr;
  int                *tmp_iptr;
  float              *tmp_fptr;
  double             *tmp_dptr;

  /***************************************************************
    Write output files using default VIC ASCII or BINARY formats
    - multiple files, all variables, no truncation

    see VIC web page for format details:
      www.hydro.washington.edu/Lettenmaier/Models/VIC/VIChome.html
  ***************************************************************/

  if(options.BINARY_OUTPUT) {  // BINARY

    // Initialize pointers
    tmp_cptr = (char *)calloc(N_OUTVAR_TYPES*options.Nlayer*options.SNOW_BAND,sizeof(char));
    tmp_siptr = (short int *)calloc(N_OUTVAR_TYPES*options.Nlayer*options.SNOW_BAND,sizeof(short int));
    tmp_usiptr = (unsigned short int *)calloc(N_OUTVAR_TYPES*options.Nlayer*options.SNOW_BAND,sizeof(unsigned short int));
    tmp_iptr = (int *)calloc(N_OUTVAR_TYPES*options.Nlayer*options.SNOW_BAND,sizeof(int));
    tmp_fptr = (float *)calloc(N_OUTVAR_TYPES*options.Nlayer*options.SNOW_BAND,sizeof(float));
    tmp_dptr = (double *)calloc(N_OUTVAR_TYPES*options.Nlayer*options.SNOW_BAND,sizeof(double));

    // Time
    tmp_iptr[0] = dmy->year;
    tmp_iptr[1] = dmy->month;
    tmp_iptr[2] = dmy->day;
    tmp_iptr[3] = dmy->hour;

    // Loop over output files
    for (file_idx = 0; file_idx < options.Noutfiles; file_idx++) {

      if (!options.OUTPUT_FORCE) {

        // Write the date
        if (dt < 24) {
          // Write year, month, day, and hour
          fwrite(tmp_iptr, sizeof(int), 4, out_data_files[file_idx].fh);
        }
        else {
          // Only write year, month, and day
          fwrite(tmp_iptr, sizeof(int), 3, out_data_files[file_idx].fh);
        }

      }

      // Loop over this output file's data variables
      for (var_idx = 0; var_idx < out_data_files[file_idx].nvars; var_idx++) {
        // Loop over this variable's elements
        ptr_idx = 0;
        if (out_data[out_data_files[file_idx].varid[var_idx]].type == OUT_TYPE_CHAR) {
          for (elem_idx = 0; elem_idx < out_data[out_data_files[file_idx].varid[var_idx]].nelem; elem_idx++) {
            tmp_cptr[ptr_idx++] = (char)out_data[out_data_files[file_idx].varid[var_idx]].aggdata[elem_idx];
          }
          fwrite(tmp_cptr, sizeof(char), ptr_idx, out_data_files[file_idx].fh);
        }
        else if (out_data[out_data_files[file_idx].varid[var_idx]].type == OUT_TYPE_SINT) {
          for (elem_idx = 0; elem_idx < out_data[out_data_files[file_idx].varid[var_idx]].nelem; elem_idx++) {
            tmp_siptr[ptr_idx++] = (short int)out_data[out_data_files[file_idx].varid[var_idx]].aggdata[elem_idx];
          }
          fwrite(tmp_siptr, sizeof(short int), ptr_idx, out_data_files[file_idx].fh);
        }
        else if (out_data[out_data_files[file_idx].varid[var_idx]].type == OUT_TYPE_USINT) {
          for (elem_idx = 0; elem_idx < out_data[out_data_files[file_idx].varid[var_idx]].nelem; elem_idx++) {
            tmp_usiptr[ptr_idx++] = (unsigned short int)out_data[out_data_files[file_idx].varid[var_idx]].aggdata[elem_idx];
          }
          fwrite(tmp_usiptr, sizeof(unsigned short int), ptr_idx, out_data_files[file_idx].fh);
        }
        else if (out_data[out_data_files[file_idx].varid[var_idx]].type == OUT_TYPE_INT) {
          for (elem_idx = 0; elem_idx < out_data[out_data_files[file_idx].varid[var_idx]].nelem; elem_idx++) {
            tmp_iptr[ptr_idx++] = (int)out_data[out_data_files[file_idx].varid[var_idx]].aggdata[elem_idx];
          }
          fwrite(tmp_iptr, sizeof(int), ptr_idx, out_data_files[file_idx].fh);
        }
        else if (out_data[out_data_files[file_idx].varid[var_idx]].type == OUT_TYPE_FLOAT) {
          for (elem_idx = 0; elem_idx < out_data[out_data_files[file_idx].varid[var_idx]].nelem; elem_idx++) {
            tmp_fptr[ptr_idx++] = (float)out_data[out_data_files[file_idx].varid[var_idx]].aggdata[elem_idx];
          }
          fwrite(tmp_fptr, sizeof(float), ptr_idx, out_data_files[file_idx].fh);
        }
        else if (out_data[out_data_files[file_idx].varid[var_idx]].type == OUT_TYPE_DOUBLE) {
          for (elem_idx = 0; elem_idx < out_data[out_data_files[file_idx].varid[var_idx]].nelem; elem_idx++) {
            tmp_dptr[ptr_idx++] = (double)out_data[out_data_files[file_idx].varid[var_idx]].aggdata[elem_idx];
          }
          fwrite(tmp_dptr, sizeof(double), ptr_idx, out_data_files[file_idx].fh);
        }
      }

    }

    // Free the arrays
    free((char *)tmp_cptr);
    free((char *)tmp_siptr);
    free((char *)tmp_usiptr);
    free((char *)tmp_iptr);
    free((char *)tmp_fptr);
    free((char *)tmp_dptr);

  }

  else {  // ASCII

    // Loop over output files
    for (file_idx = 0; file_idx < options.Noutfiles; file_idx++) {

      if (!options.OUTPUT_FORCE) {

        // Write the date
        if (dt < 24) {
          // Write year, month, day, and hour
          fprintf(out_data_files[file_idx].fh, "%04i\t%02i\t%02i\t%02i\t",
                  dmy->year, dmy->month, dmy->day, dmy->hour);
        }
        else {
          // Only write year, month, and day
          fprintf(out_data_files[file_idx].fh, "%04i\t%02i\t%02i\t",
                  dmy->year, dmy->month, dmy->day);
        }

      }

      // Loop over this output file's data variables
      for (var_idx = 0; var_idx < out_data_files[file_idx].nvars; var_idx++) {
        // Loop over this variable's elements
        for (elem_idx = 0; elem_idx < out_data[out_data_files[file_idx].varid[var_idx]].nelem; elem_idx++) {
          if (!(var_idx == 0 && elem_idx == 0)) {
            fprintf(out_data_files[file_idx].fh, "\t ");
          }
          fprintf(out_data_files[file_idx].fh, out_data[out_data_files[file_idx].varid[var_idx]].format, out_data[out_data_files[file_idx].varid[var_idx]].aggdata[elem_idx]);
        }
      }
      fprintf(out_data_files[file_idx].fh, "\n");

    }

  }

}

