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
parse_output_info(FILE           *gp,
                  stream_struct **streams,
                  dmy_struct     *dmy_current)
{
    extern option_struct       options;
    extern global_param_struct global_param;

    char                       cmdstr[MAXSTRING];
    char                       optstr[MAXSTRING];
    char                       flgstr[MAXSTRING];
    short int                  streamnum;
    char                       varname[MAXSTRING];
    int                        outvarnum;
    char                       freq_type_str[MAXSTRING];
    char                       freq_value_str[MAXSTRING];
    char                       format[MAXSTRING];
    char                       typestr[MAXSTRING];
    int                        type;
    char                       multstr[MAXSTRING];
    char                       aggstr[MAXSTRING];
    double                     mult;
    size_t                     nstreams;
    size_t                     nvars;
    unsigned short int         freq;
    int                        freq_n;
    dmy_struct                 freq_dmy;
    unsigned short int         agg_type;
    int                        found;

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
                    // determine how many variable will be in this file before
                    // allocating (GH: 209)
                    nvars = count_outfile_nvars(gp);

                    // Setup stream
                    setup_stream(streams[streamnum], nvars, 1);

                    if (sscanf(cmdstr, "%*s %s",
                               (*streams)[streamnum].prefix) != 1) {
                        log_err("Invalid specification for OUTFILE");
                    }

                    outvarnum = 0;
                }
                else if (strcasecmp("AGGFREQ", optstr) == 0) {
                    if (streamnum < 0) {
                        log_err("Error in global param file: \"OUTFILE\" must be "
                                "specified before you can specify \"AGGFREQ\".");
                    }
                    found = sscanf(cmdstr, "%*s %s %s", freq_type_str,
                                   freq_value_str);

                    if (!found) {
                        log_err("No arguments found after OUTFREQ");
                    }
                    // parse the frequency string to an enum value
                    freq = str_to_freq_flag(freq_type_str);

                    if (freq == FREQ_DATE) {
                        // Make sure we have a datestring
                        if (found != 2) {
                            log_err(
                                "AGGFREQ was set to DATE but no date string was found");
                        }
                        // parse date from freq_value_str
                        strpdmy(freq_value_str, "%Y-%m-%d", &freq_dmy);
                        // set the alarm
                        set_alarm(dmy_current, freq, &freq_dmy,
                                  (&(*streams)[streamnum].agg_alarm));
                    }
                    else {
                        if (found != 2) {
                            // Default frequency is 1
                            freq_n = 1;
                        }
                        else {
                            // get the frequency value as an integer
                            freq_n = atoi(freq_value_str);
                        }
                        // set the alarm
                        set_alarm(dmy_current, freq, &freq_n,
                                  (&(*streams)[streamnum].agg_alarm));
                    }
                }
                else if (strcasecmp("COMPRESS", optstr) == 0) {
                    if (streamnum < 0) {
                        log_err("Error in global param file: \"OUTFILE\" must be "
                                "specified before you can specify \"COMPRESS\".");
                    }
                    sscanf(cmdstr, "%*s %s", flgstr);
                    if (strcasecmp("TRUE", flgstr) == 0) {
                        (*streams)[streamnum].compress = COMPRESSION_LVL_UNSET;
                    }
                    else if (strcasecmp("FALSE", flgstr) == 0) {
                        (*streams)[streamnum].compress = 0;
                    }
                    else {
                        (*streams)[streamnum].compress = atoi(flgstr);
                    }
                }
                else if (strcasecmp("OUT_FORMAT", optstr) == 0) {
                    if (streamnum < 0) {
                        log_err("Error in global param file: \"OUTFILE\" must be "
                                "specified before you can specify \"OUT_FORMAT\".");
                    }
                    sscanf(cmdstr, "%*s %s", flgstr);
                    if (strcasecmp("ASCII", flgstr) == 0) {
                        (*streams)[streamnum].file_format = ASCII;
                    }
                    else if (strcasecmp("BINARY", flgstr) == 0) {
                        (*streams)[streamnum].file_format = BINARY;
                    }
                    else {
                        log_err("Classic driver file format must be ASCII or "
                                "BINARY [stream=%hu]", streamnum);
                    }
                }
                else if (strcasecmp("OUTVAR", optstr) == 0) {
                    if (streamnum < 0) {
                        log_err("Error in global param file: \"OUTFILE\" must be "
                                "specified before you can specify \"OUTVAR\".");
                    }
                    // parse outvar options
                    strcpy(format, "");
                    strcpy(typestr, "");
                    strcpy(multstr, "");
                    found = sscanf(cmdstr, "%*s %s %s %s %s %s", varname,
                                   format, typestr, multstr, aggstr);
                    if (!found) {
                        log_err("OUTVAR specified but no variable was listed");
                    }
                    // interpret string options, set defaults if necessary
                    str_to_ascii_format(format);
                    agg_type = str_to_agg_type(aggstr);
                    type = str_to_out_type(typestr);
                    mult = str_to_out_mult(multstr);

                    // Add OUTVAR to stream
                    set_output_var(streams[streamnum], varname, outvarnum,
                                   format, type, mult, agg_type);
                    outvarnum++;
                }
            }
            fgets(cmdstr, MAXSTRING, gp);
        }
    }
    // Otherwise, set output files and their contents to default configuration
    else {
        set_output_defaults(streams, 1, dmy_current);
    }
    fclose(gp);

    for (streamnum = 0;
         streamnum < (short int) options.Noutstreams;
         streamnum++) {
        // Allocate memory for the stream aggdata arrays
        alloc_aggdata(&(*streams)[streamnum]);
    }
}
