/******************************************************************************
 * @section DESCRIPTION
 *
 * This subroutine writes a header for all output files.
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
 * @brief    Write a header for all output files.
 *****************************************************************************/
void
write_header(out_data_file_struct *out_data_files,
             out_data_struct      *out_data,
             dmy_struct           *dmy,
             global_param_struct   global)
{
    extern option_struct options;
    size_t               file_idx;
    size_t               var_idx;
    unsigned             elem_idx;
    size_t               i;
    unsigned short       Identifier;
    unsigned short       Nbytes;
    unsigned short       Nbytes1;
    unsigned short       Nbytes2;
    char                 tmp_ALMA_OUTPUT;
    char                 Nvars;
    char                 tmp_len;
    char                *tmp_str;
    char                 tmp_type;
    float                tmp_mult;

    if (options.ALMA_OUTPUT) {
        tmp_ALMA_OUTPUT = 1;
    }
    else {
        tmp_ALMA_OUTPUT = 0;
    }

    if (options.BINARY_OUTPUT) { // BINARY
        tmp_str = (char *)calloc(256, sizeof(char));

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
        // dt          (int)*1             Output time step length in hours
        // startyear   (int)*1             Year of first record
        // startmonth  (int)*1             Month of first record
        // startday    (int)*1             Day of first record
        // starthour   (int)*1             Hour of first record
        // ALMA_OUTPUT (char)*1            0 = standard VIC units; 1 = ALMA units
        // Nvars       (char)*1            Number of variables in the file, including date fields
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

        // Loop over output files
        for (file_idx = 0; file_idx < options.Noutfiles; file_idx++) {
            // ***** Compute the number of bytes in part 1 *****

            // 1 instance of Nbytes1
            Nbytes1 = sizeof(unsigned short);

            // nrecs
            Nbytes1 += sizeof(int);

            // dt
            Nbytes1 += sizeof(int);

            // start date (year, month, day, hour)
            Nbytes1 += 4 * sizeof(int);

            // ALMA_OUTPUT
            Nbytes1 += sizeof(char);

            // Nvars
            Nbytes1 += sizeof(char);

            // ***** Compute the number of bytes in part 2 *****

            // 1 instance of Nbytes2
            Nbytes2 = sizeof(unsigned short);

            if (!options.OUTPUT_FORCE) {
                // Date fields
                Nbytes2 += sizeof(char) + 4 * sizeof(char) + sizeof(char) +
                           sizeof(float);                                        // year
                Nbytes2 += sizeof(char) + 5 * sizeof(char) + sizeof(char) +
                           sizeof(float);                                        // month
                Nbytes2 += sizeof(char) + 3 * sizeof(char) + sizeof(char) +
                           sizeof(float);                                        // day
                if (global.out_dt < HOURS_PER_DAY) {
                    Nbytes2 += sizeof(char) + 4 * sizeof(char) + sizeof(char) +
                               sizeof(float);                                      // hour
                }
            }

            // Loop over this output file's data variables
            for (var_idx = 0;
                 var_idx < out_data_files[file_idx].nvars;
                 var_idx++) {
                // Loop over this variable's elements
                for (elem_idx = 0;
                     elem_idx <
                     out_data[out_data_files[file_idx].varid[var_idx]].nelem;
                     elem_idx++) {
                    if (out_data[out_data_files[file_idx].varid[var_idx]].nelem
                        >
                        1) {
                        sprintf(tmp_str, "%s_%d",
                                out_data[out_data_files[file_idx].varid[var_idx]].varname,
                                elem_idx);
                    }
                    else {
                        strcpy(tmp_str,
                               out_data[out_data_files[file_idx].varid[var_idx]].varname);
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
                       out_data_files[file_idx].fh);
            }

            // Nbytes
            fwrite(&Nbytes, sizeof(unsigned short), 1,
                   out_data_files[file_idx].fh);

            // Nbytes1
            fwrite(&Nbytes1, sizeof(unsigned short), 1,
                   out_data_files[file_idx].fh);

            // nrecs
            fwrite(&(global.nrecs), sizeof(int), 1,
                   out_data_files[file_idx].fh);

            // dt
            fwrite(&(global.out_dt), sizeof(int), 1,
                   out_data_files[file_idx].fh);

            // start date (year, month, day, hour)
            fwrite(&(dmy->year), sizeof(int), 1, out_data_files[file_idx].fh);
            fwrite(&(dmy->month), sizeof(int), 1, out_data_files[file_idx].fh);
            fwrite(&(dmy->day), sizeof(int), 1, out_data_files[file_idx].fh);
            fwrite(&(dmy->hour), sizeof(int), 1, out_data_files[file_idx].fh);

            // ALMA_OUTPUT
            fwrite(&tmp_ALMA_OUTPUT, sizeof(char), 1,
                   out_data_files[file_idx].fh);

            // Nvars
            Nvars = out_data_files[file_idx].nvars;
            if (!options.OUTPUT_FORCE) {
                if (global.out_dt < HOURS_PER_DAY) {
                    Nvars += 4;
                }
                else {
                    Nvars += 3;
                }
            }
            fwrite(&Nvars, sizeof(char), 1, out_data_files[file_idx].fh);

            // Nbytes2
            fwrite(&Nbytes2, sizeof(unsigned short), 1,
                   out_data_files[file_idx].fh);

            if (!options.OUTPUT_FORCE) {
                // Date fields
                tmp_type = OUT_TYPE_INT;
                tmp_mult = 1.;

                // year
                strcpy(tmp_str, "YEAR");
                tmp_len = strlen(tmp_str);
                fwrite(&tmp_len, sizeof(char), 1, out_data_files[file_idx].fh);
                fwrite(tmp_str, sizeof(char), tmp_len,
                       out_data_files[file_idx].fh);
                fwrite(&tmp_type, sizeof(char), 1, out_data_files[file_idx].fh);
                fwrite(&tmp_mult, sizeof(float), 1,
                       out_data_files[file_idx].fh);

                // month
                strcpy(tmp_str, "MONTH");
                tmp_len = strlen(tmp_str);
                fwrite(&tmp_len, sizeof(char), 1, out_data_files[file_idx].fh);
                fwrite(tmp_str, sizeof(char), tmp_len,
                       out_data_files[file_idx].fh);
                fwrite(&tmp_type, sizeof(char), 1, out_data_files[file_idx].fh);
                fwrite(&tmp_mult, sizeof(float), 1,
                       out_data_files[file_idx].fh);

                // day
                strcpy(tmp_str, "DAY");
                tmp_len = strlen(tmp_str);
                fwrite(&tmp_len, sizeof(char), 1, out_data_files[file_idx].fh);
                fwrite(tmp_str, sizeof(char), tmp_len,
                       out_data_files[file_idx].fh);
                fwrite(&tmp_type, sizeof(char), 1, out_data_files[file_idx].fh);
                fwrite(&tmp_mult, sizeof(float), 1,
                       out_data_files[file_idx].fh);

                if (global.out_dt < HOURS_PER_DAY) {
                    // hour
                    strcpy(tmp_str, "HOUR");
                    tmp_len = strlen(tmp_str);
                    fwrite(&tmp_len, sizeof(char), 1,
                           out_data_files[file_idx].fh);
                    fwrite(tmp_str, sizeof(char), tmp_len,
                           out_data_files[file_idx].fh);
                    fwrite(&tmp_type, sizeof(char), 1,
                           out_data_files[file_idx].fh);
                    fwrite(&tmp_mult, sizeof(float), 1,
                           out_data_files[file_idx].fh);
                }
            }

            // Loop over this output file's data variables
            for (var_idx = 0;
                 var_idx < out_data_files[file_idx].nvars;
                 var_idx++) {
                // Loop over this variable's elements
                for (elem_idx = 0;
                     elem_idx <
                     out_data[out_data_files[file_idx].varid[var_idx]].nelem;
                     elem_idx++) {
                    if (out_data[out_data_files[file_idx].varid[var_idx]].nelem
                        >
                        1) {
                        sprintf(tmp_str, "%s_%d",
                                out_data[out_data_files[file_idx].varid[var_idx]].varname,
                                elem_idx);
                    }
                    else {
                        strcpy(tmp_str,
                               out_data[out_data_files[file_idx].varid[var_idx]].varname);
                    }
                    tmp_len = strlen(tmp_str);
                    fwrite(&tmp_len, sizeof(char), 1,
                           out_data_files[file_idx].fh);
                    fwrite(tmp_str, sizeof(char), tmp_len,
                           out_data_files[file_idx].fh);
                    tmp_type =
                        out_data[out_data_files[file_idx].varid[var_idx]].type;
                    fwrite(&tmp_type, sizeof(char), 1,
                           out_data_files[file_idx].fh);
                    tmp_mult =
                        out_data[out_data_files[file_idx].varid[var_idx]].mult;
                    fwrite(&tmp_mult, sizeof(float), 1,
                           out_data_files[file_idx].fh);
                }
            }
        }
    }
    else { // ASCII
           // ASCII header format:
           //
           // # NRECS: (nrecs)
           // # DT: (dt)
           // # STARTDATE: yyyy-mm-dd hh:00:00
           // # ALMA_OUTPUT: (0 or 1)
           // # NVARS: (Nvars)
           // # VARNAME    VARNAME   VARNAME   ...
           //
           // where
           // nrecs       = Number of records in the file
           // dt          = Output time step length in hours
           // start date  = Date and time of first record of file
           // ALMA_OUTPUT = Indicates units of the variables; 0 = standard VIC units; 1 = ALMA units
           // Nvars       = Number of variables in the file, including date fields

        // Loop over output files
        for (file_idx = 0; file_idx < options.Noutfiles; file_idx++) {
            // Header part 1: Global attributes
            Nvars = out_data_files[file_idx].nvars;
            if (!options.OUTPUT_FORCE) {
                if (global.out_dt < HOURS_PER_DAY) {
                    Nvars += 4;
                }
                else {
                    Nvars += 3;
                }
            }
            fprintf(out_data_files[file_idx].fh, "# NRECS: %d\n", global.nrecs);
            fprintf(out_data_files[file_idx].fh, "# DT: %d\n", global.out_dt);
            fprintf(out_data_files[file_idx].fh,
                    "# STARTDATE: %04d-%02d-%02d %02d:00:00\n",
                    dmy->year, dmy->month, dmy->day, dmy->hour);
            fprintf(out_data_files[file_idx].fh, "# ALMA_OUTPUT: %d\n",
                    tmp_ALMA_OUTPUT);
            fprintf(out_data_files[file_idx].fh, "# NVARS: %d\n", Nvars);

            // Header part 2: Variables
            fprintf(out_data_files[file_idx].fh, "# ");

            if (!options.OUTPUT_FORCE) {
                // Write the date
                if (global.out_dt < HOURS_PER_DAY) {
                    // Write year, month, day, and hour
                    fprintf(out_data_files[file_idx].fh,
                            "YEAR\tMONTH\tDAY\tHOUR\t");
                }
                else {
                    // Only write year, month, and day
                    fprintf(out_data_files[file_idx].fh, "YEAR\tMONTH\tDAY\t");
                }
            }

            // Loop over this output file's data variables
            for (var_idx = 0;
                 var_idx < out_data_files[file_idx].nvars;
                 var_idx++) {
                // Loop over this variable's elements
                for (elem_idx = 0;
                     elem_idx <
                     out_data[out_data_files[file_idx].varid[var_idx]].nelem;
                     elem_idx++) {
                    if (!(var_idx == 0 && elem_idx == 0)) {
                        fprintf(out_data_files[file_idx].fh, "\t ");
                    }
                    fprintf(out_data_files[file_idx].fh, "%s",
                            out_data[out_data_files[file_idx].varid[var_idx]].varname);
                    if (out_data[out_data_files[file_idx].varid[var_idx]].nelem
                        >
                        1) {
                        fprintf(out_data_files[file_idx].fh, "_%d", elem_idx);
                    }
                }
            }
            fprintf(out_data_files[file_idx].fh, "\n");
        }
    }
}
