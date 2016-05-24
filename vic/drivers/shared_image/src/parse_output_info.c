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

#include <vic_driver_shared_image.h>

/******************************************************************************
 * @brief    Get output info from global parameter file.
 *****************************************************************************/
void
parse_output_info(FILE           *gp,
                  stream_struct **streams)
{
    extern option_struct       options;
    extern global_param_struct global_param;

    char                       cmdstr[MAXSTRING];
    char                       optstr[MAXSTRING];
    char                       flgstr[MAXSTRING];
    short int                  streamnum;
    char                       varname[MAXSTRING];
    int                        outvarnum;
    char                       format[MAXSTRING];
    char                       typestr[MAXSTRING];
    int                        type;
    char                       multstr[MAXSTRING];
    char                       aggstr[MAXSTRING];
    double                     mult;
    size_t                     nstreams;
    size_t                     nvars;
    unsigned int               nextagg;
    unsigned short int         agg_type;

    /** Read through global control file to find output info **/

    fgets(cmdstr, MAXSTRING, gp);

    streamnum = -1;
    outvarnum = 0;

    // Count the number of output files listed in the global param file
    nstreams = count_n_outfiles(gp);

    // only parse the output info if there are output files to parse
    if (nstreams > 0) {
        options.Noutstreams = nstreams;

        *streams = calloc(options.Noutstreams, sizeof(*(*streams)));
        if (*streams == NULL) {
            log_err("Memory allocation error in parse_output_info().");
        }

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
                    sscanf(cmdstr, "%*s %s",
                           (*streams)[streamnum].prefix);

                    // determine how many variable will be in this file before
                    // allocating (GH: 209)
                    nvars = count_outfile_nvars(gp);

                    // TODO: parse stream specific agg period info from global param
                    // nextagg = (unsigned int) ((int) global_param.out_dt / (int) global_param.dt);
                    nextagg = 1;

                    setup_stream(streams[streamnum], nvars);

                    outvarnum = 0;
                }
                else if (strcasecmp("OUTPUT_STEPS_PER_DAY", optstr) == 0) {
                    sscanf(cmdstr, "%*s %zu",
                           (&(*streams)[streamnum].output_steps_per_day));
                }
                else if (strcasecmp("COMPRESS", optstr) == 0) {
                    sscanf(cmdstr, "%*s %s", flgstr);
                    if (strcasecmp("TRUE", flgstr) == 0) {
                        (*streams)[streamnum].compress =
                            COMPRESSION_LVL_DEFAULT;
                    }
                    else if (strcasecmp("FALSE", flgstr) == 0) {
                        (*streams)[streamnum].compress = 0;
                    }
                    else {
                        (*streams)[streamnum].compress = atoi(flgstr);
                    }
                }
                else if (strcasecmp("OUT_FORMAT", optstr) == 0) {
                    sscanf(cmdstr, "%*s %s", flgstr);
                    if (strcasecmp("ASCII", flgstr) == 0) {
                        (*streams)[streamnum].file_format =
                            COMPRESSION_LVL_DEFAULT;
                    }
                    else if (strcasecmp("BINARY", flgstr) == 0) {
                        (*streams)[streamnum].file_format = 0;
                    }
                    else {
                        log_err(
                            "File format must be ASCII or BINARY [stream=%hu]",
                            streamnum);
                    }
                }
                else if (strcasecmp("OUTVAR", optstr) == 0) {
                    if (streamnum < 0) {
                        log_err("Error in global param file: \"OUTFILE\" must be "
                                "specified before you can specify \"OUTVAR\".");
                    }
                    strcpy(format, "");
                    strcpy(typestr, "");
                    strcpy(multstr, "");
                    sscanf(cmdstr, "%*s %s %s %s %s %s", varname, format, typestr,
                           multstr, aggstr);
                    if (strcasecmp("", format) == 0) {
                        strcpy(format, "*");
                    }

                    agg_type = str_to_agg_type(aggstr);
                    type = str_to_out_type(typestr);
                    mult = str_to_out_mult(multstr);

                    set_output_var(streams[streamnum],
                                   varname, outvarnum, format, type, mult,
                                   agg_type);
                    strcpy(format, "");
                    outvarnum++;
                }
            }
            fgets(cmdstr, MAXSTRING, gp);
        }
    }
    // Otherwise, set output files and their contents to default configuration
    else {
        set_output_defaults(streams);
    }
    fclose(gp);

    for (streamnum = 0; streamnum < (short int) options.Noutstreams; streamnum++) {
        // Validate the streams
        validate_stream_settings(&((*streams)[streamnum]));

        // Allocate memory for the stream aggdata arrays
        alloc_aggdata(&(*streams)[streamnum]);
    }
}
