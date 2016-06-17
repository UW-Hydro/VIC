/******************************************************************************
 * @section DESCRIPTION
 *
 * This subroutine writes a header for all output files.
 *
 * @section LICENSE
 *
 * The Variable Infiltration Capacity (VIC) macroscale hydrological model
 * Copyright (C) 2016 The Computational Hydrology Group, Department of Civil
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
 * @brief    Write a header for all output files.
 *****************************************************************************/
void
write_header(stream_struct **streams,
             dmy_struct     *dmy)
{
    extern option_struct       options;
    extern global_param_struct global_param;
    extern metadata_struct     out_metadata[N_OUTVAR_TYPES];

    size_t                     stream_idx;
    size_t                     var_idx;
    unsigned                   elem_idx;
    size_t                     i;
    unsigned int               varid;
    unsigned short int         Identifier;
    unsigned short int         Nbytes;
    unsigned short int         Nbytes1;
    unsigned short int         Nbytes2;
    size_t                     nvars;
    char                       tmp_len;
    char                      *tmp_str;
    char                       tmp_type;
    float                      tmp_mult;


    // Loop over output files
    for (stream_idx = 0; stream_idx < options.Noutstreams; stream_idx++) {
        if ((*streams)[stream_idx].file_format == BINARY) {
            tmp_str = calloc(BINHEADERSIZE, sizeof(*tmp_str));

            // Binary header format:
            //
            // Data        Stored As           Comment
            //
            // Identifier  (unsigned short)*4  0xFFFF, repeated 4 times
            // Nbytes      (unsigned short)*1  Number of bytes in the header,
            // INCLUDING THE IDENTIFIER
            //
            // Part 1: Global Attributes
            // Nbytes1     (unsigned short)*1  Number of bytes in part 1
            // nrecs       (int)*1             Number of records in the file
            // startyear   (int)*1             Year of first record
            // startmonth  (int)*1             Month of first record
            // startday    (int)*1             Day of first record
            // startsec    (int)*1             Second of first record
            // nvars       (char)*1            Number of variables in the file, including date fields
            //
            // Part 2: Variables
            // Nbytes2     (unsigned short)*1  Number of bytes in part 2
            // For each variable, the following fields: { len varname type mult }
            // len       (char)*1            Number of characters in varname
            // varname   (char)*len          Variable name
            // type      (char)*1            Code identifying variable type
            // mult      (float)*1           Multiplier for variable

            // Identifier
            Identifier = 0xFFFF;

            // ***** Compute the number of bytes in part 1 *****

            // 1 instance of Nbytes1
            Nbytes1 = sizeof(unsigned short int);

            // nrecs
            Nbytes1 += sizeof(size_t);

            // start date (year, month, day, sec)
            Nbytes1 += sizeof(int) + 2 * sizeof(unsigned short int) +
                       sizeof(unsigned int);

            // nvars
            Nbytes1 += sizeof(size_t);

            // ***** Compute the number of bytes in part 2 *****

            // 1 instance of Nbytes2
            Nbytes2 = sizeof(unsigned short);

            // Date fields
            Nbytes2 += sizeof(char) + 4 * sizeof(char) + sizeof(char) +
                       sizeof(float);                                        // year
            Nbytes2 += sizeof(char) + 5 * sizeof(char) + sizeof(char) +
                       sizeof(float);                                        // month
            Nbytes2 += sizeof(char) + 3 * sizeof(char) + sizeof(char) +
                       sizeof(float);                                        // day
            if ((*streams)[stream_idx].agg_alarm.is_subdaily) {
                Nbytes2 += sizeof(char) + 4 * sizeof(char) + sizeof(char) +
                           sizeof(float);                                      // sec
            }

            // Loop over this output file's data variables
            for (var_idx = 0; var_idx < (*streams)[stream_idx].nvars;
                 var_idx++) {
                varid = (*streams)[stream_idx].varid[var_idx];
                // Loop over this variable's elements
                for (elem_idx = 0;
                     elem_idx < out_metadata[varid].nelem;
                     elem_idx++) {
                    if (out_metadata[varid].nelem > 1) {
                        sprintf(tmp_str, "%s_%d", out_metadata[varid].varname,
                                elem_idx);
                    }
                    else {
                        strcpy(tmp_str,
                               out_metadata[varid].varname);
                    }
                    Nbytes2 += sizeof(char) + strlen(tmp_str) * sizeof(char) +
                               sizeof(char) + sizeof(float);
                }
            }

            // ***** Compute the total number of bytes in the header *****

            // 4 instances of Identifier, plus 1 instance of Nbytes, plus number of bytes in parts 1 and 2
            Nbytes = 4 * sizeof(unsigned short) + sizeof(unsigned short) +
                     Nbytes1 + Nbytes2;

            // ***** Write the header *****

            // 4 instances of Identifier
            for (i = 0; i < 4; i++) {
                fwrite(&Identifier, sizeof(unsigned short), 1,
                       (*streams)[stream_idx].fh);
            }

            // Nbytes
            fwrite(&Nbytes, sizeof(unsigned short), 1,
                   (*streams)[stream_idx].fh);

            // Nbytes1
            fwrite(&Nbytes1, sizeof(unsigned short), 1,
                   (*streams)[stream_idx].fh);

            // nrecs
            fwrite(&(global_param.nrecs), sizeof(size_t), 1,
                   (*streams)[stream_idx].fh);

            // start date (year, month, day, sec)
            fwrite(&(dmy->year), sizeof(int), 1,
                   (*streams)[stream_idx].fh);
            fwrite(&(dmy->month), sizeof(unsigned short int), 1,
                   (*streams)[stream_idx].fh);
            fwrite(&(dmy->day), sizeof(unsigned short int), 1,
                   (*streams)[stream_idx].fh);
            fwrite(&(dmy->dayseconds), sizeof(unsigned int), 1,
                   (*streams)[stream_idx].fh);

            // nvars
            nvars = (*streams)[stream_idx].nvars;
            if ((*streams)[stream_idx].agg_alarm.is_subdaily) {
                nvars += 4;
            }
            else {
                nvars += 3;
            }
            fwrite(&nvars, sizeof(size_t), 1, (*streams)[stream_idx].fh);

            // Nbytes2
            fwrite(&Nbytes2, sizeof(unsigned short), 1,
                   (*streams)[stream_idx].fh);

            // Date fields
            tmp_type = OUT_TYPE_INT;
            tmp_mult = 1.;

            // year
            strcpy(tmp_str, "YEAR");
            tmp_len = strlen(tmp_str);
            fwrite(&tmp_len, sizeof(char), 1, (*streams)[stream_idx].fh);
            fwrite(tmp_str, sizeof(char), tmp_len,
                   (*streams)[stream_idx].fh);
            fwrite(&tmp_type, sizeof(char), 1,
                   (*streams)[stream_idx].fh);
            fwrite(&tmp_mult, sizeof(float), 1,
                   (*streams)[stream_idx].fh);

            // month
            strcpy(tmp_str, "MONTH");
            tmp_len = strlen(tmp_str);
            fwrite(&tmp_len, sizeof(char), 1, (*streams)[stream_idx].fh);
            fwrite(tmp_str, sizeof(char), tmp_len,
                   (*streams)[stream_idx].fh);
            fwrite(&tmp_type, sizeof(char), 1,
                   (*streams)[stream_idx].fh);
            fwrite(&tmp_mult, sizeof(float), 1,
                   (*streams)[stream_idx].fh);

            // day
            strcpy(tmp_str, "DAY");
            tmp_len = strlen(tmp_str);
            fwrite(&tmp_len, sizeof(char), 1, (*streams)[stream_idx].fh);
            fwrite(tmp_str, sizeof(char), tmp_len,
                   (*streams)[stream_idx].fh);
            fwrite(&tmp_type, sizeof(char), 1,
                   (*streams)[stream_idx].fh);
            fwrite(&tmp_mult, sizeof(float), 1,
                   (*streams)[stream_idx].fh);

            if ((*streams)[stream_idx].agg_alarm.is_subdaily) {
                // sec
                strcpy(tmp_str, "SEC");
                tmp_len = strlen(tmp_str);
                fwrite(&tmp_len, sizeof(char), 1,
                       (*streams)[stream_idx].fh);
                fwrite(tmp_str, sizeof(char), tmp_len,
                       (*streams)[stream_idx].fh);
                fwrite(&tmp_type, sizeof(char), 1,
                       (*streams)[stream_idx].fh);
                fwrite(&tmp_mult, sizeof(float), 1,
                       (*streams)[stream_idx].fh);
            }

            // Loop over this output file's data variables
            for (var_idx = 0; var_idx < (*streams)[stream_idx].nvars;
                 var_idx++) {
                varid = (*streams)[stream_idx].varid[var_idx];
                // Loop over this variable's elements
                for (elem_idx = 0;
                     elem_idx < out_metadata[varid].nelem;
                     elem_idx++) {
                    if (out_metadata[varid].nelem > 1) {
                        sprintf(tmp_str, "%s_%d", out_metadata[varid].varname,
                                elem_idx);
                    }
                    else {
                        strcpy(tmp_str, out_metadata[varid].varname);
                    }
                    tmp_len = strlen(tmp_str);
                    fwrite(&tmp_len, sizeof(char), 1,
                           (*streams)[stream_idx].fh);
                    fwrite(tmp_str, sizeof(char), tmp_len,
                           (*streams)[stream_idx].fh);
                    tmp_type = (*streams)[stream_idx].type[var_idx];
                    fwrite(&tmp_type, sizeof(char), 1,
                           (*streams)[stream_idx].fh);
                    tmp_mult = (*streams)[stream_idx].mult[var_idx];
                    fwrite(&tmp_mult, sizeof(float), 1,
                           (*streams)[stream_idx].fh);
                }
            }
        }
        else if ((*streams)[stream_idx].file_format == ASCII) {
            // ASCII header format:
            //
            // # NRECS: (nrecs)
            // # STARTDATE: yyyy-mm-dd hh:00:00
            // # NVARS: (nvars)
            // # VARNAME    VARNAME   VARNAME   ...
            //
            // where
            // nrecs       = Number of records in the file
            // start date  = Date and time of first record of file
            // nvars       = Number of variables in the file, including date fields

            // Header part 1: Global attributes
            nvars = (*streams)[stream_idx].nvars;
            if ((*streams)[stream_idx].agg_alarm.is_subdaily) {
                nvars += 4;
            }
            else {
                nvars += 3;
            }
            fprintf((*streams)[stream_idx].fh, "# SIMULATION: %s\n",
                    (*streams)[stream_idx].prefix);
            fprintf((*streams)[stream_idx].fh, "# MODEL_VERSION: %s\n",
                    SHORT_VERSION);

            // Header part 2: Variables
            // Write the date
            if ((*streams)[stream_idx].agg_alarm.is_subdaily) {
                // Write year, month, day, and sec
                fprintf((*streams)[stream_idx].fh,
                        "YEAR\tMONTH\tDAY\tSEC\t");
            }
            else {
                // Only write year, month, and day
                fprintf((*streams)[stream_idx].fh, "YEAR\tMONTH\tDAY\t");
            }

            // Loop over this output file's data variables
            for (var_idx = 0; var_idx < (*streams)[stream_idx].nvars;
                 var_idx++) {
                varid = (*streams)[stream_idx].varid[var_idx];
                // Loop over this variable's elements
                for (elem_idx = 0;
                     elem_idx < out_metadata[varid].nelem;
                     elem_idx++) {
                    if (!(var_idx == 0 && elem_idx == 0)) {
                        fprintf((*streams)[stream_idx].fh, "\t ");
                    }
                    fprintf((*streams)[stream_idx].fh, "%s",
                            out_metadata[varid].varname);
                    if (out_metadata[varid].nelem > 1) {
                        fprintf((*streams)[stream_idx].fh, "_%d",
                                elem_idx);
                    }
                }
            }
            fprintf((*streams)[stream_idx].fh, "\n");
        }
        else {
            log_err("Unrecognized OUT_FORMAT option");
        }
    }
}
