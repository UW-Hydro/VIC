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

#include <vic_driver_classic.h>

/******************************************************************************
 * @brief    Get output info from global parameter file.
 *****************************************************************************/
void
parse_output_info(FILE               *gp,
                  stream_struct      **output_streams,
                  stream_file_struct **out_data_files)
{
    extern option_struct       options;
    extern global_param_struct global_param;

    char                 cmdstr[MAXSTRING];
    char                 optstr[MAXSTRING];
    short int            streamnum;
    char                 varname[MAXSTRING];
    int                  outvarnum;
    char                 format[MAXSTRING];
    char                 typestr[MAXSTRING];
    int                  type;
    char                 multstr[MAXSTRING];
    double               mult;
    size_t               nstreams;
    size_t               nvars;
    unsigned int         nextagg;

    strcpy(format, "*");

    /** Read through global control file to find output info **/

    fgets(cmdstr, MAXSTRING, gp);

    streamnum = -1;
    outvarnum = 0;

    // Count the number of output files listed in the global param file
    nstreams = count_n_outfiles(gp);

    // only parse the output info if there are output files to parse
    if (nstreams > 0) {
        options.Noutstreams = nstreams;

        output_streams = calloc(options.Noutstreams, sizeof(*(output_streams)));
        if (output_streams == NULL) {
            log_err("Memory allocation error in parse_output_info().");
        }
        out_data_files = calloc(options.Noutstreams, sizeof(*(out_data_files)));
        if (out_data_files == NULL) {
            log_err("Memory allocation error in parse_output_info().");
        }

        // PRT_SNOW_BAND is ignored if options.Noutstreams > 0
        options.PRT_SNOW_BAND = false;

        while (!feof(gp)) {
            if (cmdstr[0] != '#' && cmdstr[0] != '\n' && cmdstr[0] != '\0') {
                sscanf(cmdstr, "%s", optstr);

                if (strcasecmp("OUTFILE", optstr) == 0) {
                    streamnum++;
                    if (streamnum >= (short int) options.Noutstreams) {
                        log_err("Found too many output files, was expecting "
                                "%zu but found %hu", options.Noutstreams,
                                streamnum);
                    }
                    sscanf(cmdstr, "%*s %s", out_data_files[streamnum]->prefix);

                    // determine how many variable will be in this file before
                    // allocating (GH: 209)
                    nvars = count_outfile_nvars(gp);

                    // TODO: parse stream specific agg period info from global param
                    nextagg = (unsigned int) ((int) global_param.out_dt / (int) global_param.dt);

                    setup_stream(output_streams[streamnum],
                                 out_data_files[streamnum],
                                 nvars, nextagg);

                    outvarnum = 0;
                }
                else if (strcasecmp("OUTVAR", optstr) == 0) {
                    if (streamnum < 0) {
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
                    set_output_var(output_streams[streamnum],
                                   out_data_files[streamnum],
                                   varname, outvarnum, format, type, mult);
                    strcpy(format, "");
                    outvarnum++;
                }
            }
            fgets(cmdstr, MAXSTRING, gp);
        }
    }
    // Otherwise, set output files and their contents to default configuration
    else {
        set_output_defaults(output_streams, out_data_files);
    }
    fclose(gp);
}

/******************************************************************************
 * @brief    This routine determines the counts the number of output variables
             in each output file specified in the global parameter file.
 *****************************************************************************/
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

/******************************************************************************
 * @brief    This routine determines the counts the number of output files
             specified in the global parameter file.
 *****************************************************************************/
size_t
count_n_outfiles(FILE *gp)
{
    size_t        n_outfiles;
    unsigned long start_position;
    char          cmdstr[MAXSTRING];
    char          optstr[MAXSTRING];
    // Figure out where we are in the input file
    fflush(gp);
    start_position = ftell(gp);

    // read the first line
    fgets(cmdstr, MAXSTRING, gp);

    // initalize n_outfiles
    n_outfiles = 0;

    // Loop through the lines
    while (!feof(gp)) {
        if (cmdstr[0] != '#' && cmdstr[0] != '\n' && cmdstr[0] != '\0') {
            // line is not blank or a comment
            sscanf(cmdstr, "%s", optstr);

            // if the line starts with OUTFILE
            if (strcasecmp("OUTFILE", optstr) == 0) {
                n_outfiles++;
            }
        }
        fgets(cmdstr, MAXSTRING, gp);
    }

    // put the position in the file back to where we started
    fseek(gp, start_position, SEEK_SET);

    return n_outfiles;
}
