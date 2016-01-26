/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine reads the VIC model global control file, getting information
 * for output variables list (if any).
 *
 * @section LICENSE
 *
 * The Variable Infiltration Capacity (VIC) macroscale hydrological model
 * Copyright (C) 2014 The Land Surface Hydrology Group, Department of Civil
 * and Environmental Engineering, University of Washington.
 *
 * The VIC model is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 *****************************************************************************/

#include <vic_def.h>
#include <vic_run.h>
#include <vic_driver_classic.h>

/******************************************************************************
 * @brief    Get output info from global parameter file.
 *****************************************************************************/
void
parse_output_info(FILE                  *gp,
                  out_data_file_struct **out_data_files,
                  out_data_struct       *out_data)
{
    extern option_struct options;

    char                 cmdstr[MAXSTRING];
    char                 optstr[MAXSTRING];
    short int            outfilenum;
    char                 varname[MAXSTRING];
    int                  outvarnum;
    char                 format[MAXSTRING];
    char                 typestr[MAXSTRING];
    int                  type;
    char                 multstr[MAXSTRING];
    double               mult;
    int                  tmp_noutfiles;

    strcpy(format, "*");

    /** Read through global control file to find output info **/

    fgets(cmdstr, MAXSTRING, gp);

    outfilenum = -1;
    outvarnum = 0;
    while (!feof(gp)) {
        if (cmdstr[0] != '#' && cmdstr[0] != '\n' && cmdstr[0] != '\0') {
            sscanf(cmdstr, "%s", optstr);

            if (strcasecmp("N_OUTFILES", optstr) == 0) {
                sscanf(cmdstr, "%*s %d", &tmp_noutfiles);
                free_out_data_files(out_data_files);
                options.Noutfiles = tmp_noutfiles;
                *out_data_files = calloc(options.Noutfiles,
                                         sizeof(*(*out_data_files)));
                outfilenum = -1;
                init_output_list(out_data, false, "%.4f", OUT_TYPE_FLOAT, 1);
                // PRT_SNOW_BAND is ignored if N_OUTFILES has been specified
                options.PRT_SNOW_BAND = false;
            }
            else if (strcasecmp("OUTFILE", optstr) == 0) {
                outfilenum++;
                if (!options.Noutfiles) {
                    log_err("in global param file: \"N_OUTFILES\" must be "
                            "specified before you can specify \"OUTFILE\".");
                }
                if (outfilenum >= (short int)options.Noutfiles) {
                    log_err("number of output files specified in N_OUTFILES "
                            "(%zu) is less than actual number of output files "
                            "defined in the global param file.",
                            options.Noutfiles);
                }
                sscanf(cmdstr, "%*s %s", (*out_data_files)[outfilenum].prefix);

                // determine how many variable will be in this file before
                // allocating (GH: 209)
                (*out_data_files)[outfilenum].nvars = count_outfile_nvars(gp);

                (*out_data_files)[outfilenum].varid =
                    calloc((*out_data_files)[outfilenum].nvars,
                           sizeof(*((*out_data_files)[outfilenum].varid)));
                outvarnum = 0;
            }
            else if (strcasecmp("OUTVAR", optstr) == 0) {
                if (outfilenum < 0) {
                    log_err("Error in global param file: \"OUTFILE\" must be "
                            "specified before you can specify \"OUTVAR\".");
                }
                strcpy(format, "");
                strcpy(typestr, "");
                strcpy(multstr, "");
                sscanf(cmdstr, "%*s %s %s %s %s", varname, format, typestr,
                       multstr);
                if (strcasecmp("", format) == 0) {
                    strcpy(format, "*");
                    type = OUT_TYPE_DEFAULT;
                    mult = 0; // 0 means default multiplier
                }
                else {
                    if (strcasecmp("OUT_TYPE_USINT", typestr) == 0) {
                        type = OUT_TYPE_USINT;
                    }
                    else if (strcasecmp("OUT_TYPE_SINT", typestr) == 0) {
                        type = OUT_TYPE_SINT;
                    }
                    else if (strcasecmp("OUT_TYPE_FLOAT", typestr) == 0) {
                        type = OUT_TYPE_FLOAT;
                    }
                    else if (strcasecmp("OUT_TYPE_DOUBLE", typestr) == 0) {
                        type = OUT_TYPE_DOUBLE;
                    }
                    else {
                        type = OUT_TYPE_DEFAULT;
                    }
                    if (strcmp("*", multstr) == 0) {
                        mult = 0; // 0 means use default multiplier
                    }
                    else {
                        mult = (double) atof(multstr);
                    }
                }
                if (set_output_var((*out_data_files), true, outfilenum,
                                   out_data, varname, outvarnum, format, type,
                                   mult) != 0) {
                    log_err("Invalid output variable specification.");
                }
                strcpy(format, "");
                outvarnum++;
            }
        }
        fgets(cmdstr, MAXSTRING, gp);
    }
    fclose(gp);
}

size_t
count_outfile_nvars(FILE *gp)
{
    size_t        nvars;
    unsigned long start_position;
    char          cmdstr[MAXSTRING];
    char          optstr[MAXSTRING];
    // Figure out where we are in the input file
    fflush(gp);
    start_position = ftell(gp);

    // read the first line
    fgets(cmdstr, MAXSTRING, gp);

    // initalize nvars
    nvars = 0;

    // Loop through the lines
    while (!feof(gp)) {
        if (cmdstr[0] != '#' && cmdstr[0] != '\n' && cmdstr[0] != '\0') {
            // line is not blank or a comment
            sscanf(cmdstr, "%s", optstr);

            // if the line starts with OUTFILE
            if (strcasecmp("OUTVAR", optstr) == 0) {
                nvars++;
            }
            // else we're done with this file so break out of loop
            else {
                break;
            }
        }
        fgets(cmdstr, MAXSTRING, gp);
    }

    // put the position in the file back to where we started
    fseek(gp, start_position, SEEK_SET);

    return nvars;
}
