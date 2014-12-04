/******************************************************************************
 * @section DESCRIPTION
 *
 * This subroutine writes all output variables to output files.
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
 * @brief    write all variables to output files.
 *****************************************************************************/
void
write_data(out_data_file_struct *out_data_files,
           out_data_struct      *out_data,
           dmy_struct           *dmy,
           int                   dt)
{
    extern option_struct options;
    size_t               file_idx;
    size_t               var_idx;
    size_t               elem_idx;
    size_t               ptr_idx;
    char                *tmp_cptr;
    short int           *tmp_siptr;
    unsigned short int  *tmp_usiptr;
    int                 *tmp_iptr;
    float               *tmp_fptr;
    double              *tmp_dptr;

    if (options.BINARY_OUTPUT) { // BINARY
        // Initialize pointers
        tmp_cptr = (char *)calloc(
            N_OUTVAR_TYPES * options.Nlayer * options.SNOW_BAND, sizeof(char));
        tmp_siptr = (short int *)calloc(
            N_OUTVAR_TYPES * options.Nlayer * options.SNOW_BAND,
            sizeof(short int));
        tmp_usiptr = (unsigned short int *)calloc(
            N_OUTVAR_TYPES * options.Nlayer * options.SNOW_BAND,
            sizeof(unsigned short int));
        tmp_iptr = (int *)calloc(
            N_OUTVAR_TYPES * options.Nlayer * options.SNOW_BAND, sizeof(int));
        tmp_fptr = (float *)calloc(
            N_OUTVAR_TYPES * options.Nlayer * options.SNOW_BAND, sizeof(float));
        tmp_dptr = (double *)calloc(
            N_OUTVAR_TYPES * options.Nlayer * options.SNOW_BAND,
            sizeof(double));

        // Time
        tmp_iptr[0] = dmy->year;
        tmp_iptr[1] = dmy->month;
        tmp_iptr[2] = dmy->day;
        tmp_iptr[3] = dmy->hour;

        // Loop over output files
        for (file_idx = 0; file_idx < options.Noutfiles; file_idx++) {
            if (!options.OUTPUT_FORCE) {
                // Write the date
                if (dt < HOURS_PER_DAY) {
                    // Write year, month, day, and hour
                    fwrite(tmp_iptr, sizeof(int), 4,
                           out_data_files[file_idx].fh);
                }
                else {
                    // Only write year, month, and day
                    fwrite(tmp_iptr, sizeof(int), 3,
                           out_data_files[file_idx].fh);
                }
            }

            // Loop over this output file's data variables
            for (var_idx = 0;
                 var_idx < out_data_files[file_idx].nvars;
                 var_idx++) {
                // Loop over this variable's elements
                ptr_idx = 0;
                if (out_data[out_data_files[file_idx].varid[var_idx]].type ==
                    OUT_TYPE_CHAR) {
                    for (elem_idx = 0;
                         elem_idx <
                         out_data[out_data_files[file_idx].varid[var_idx]].nelem;
                         elem_idx++) {
                        tmp_cptr[ptr_idx++] =
                            (char)out_data[out_data_files[file_idx].varid[
                                               var_idx]].
                            aggdata[elem_idx];
                    }
                    fwrite(tmp_cptr, sizeof(char), ptr_idx,
                           out_data_files[file_idx].fh);
                }
                else if (out_data[out_data_files[file_idx].varid[var_idx]].type
                         ==
                         OUT_TYPE_SINT) {
                    for (elem_idx = 0;
                         elem_idx <
                         out_data[out_data_files[file_idx].varid[var_idx]].nelem;
                         elem_idx++) {
                        tmp_siptr[ptr_idx++] =
                            (short int)out_data[out_data_files[file_idx].varid[
                                                    var_idx]].aggdata[elem_idx];
                    }
                    fwrite(tmp_siptr, sizeof(short int), ptr_idx,
                           out_data_files[file_idx].fh);
                }
                else if (out_data[out_data_files[file_idx].varid[var_idx]].type
                         ==
                         OUT_TYPE_USINT) {
                    for (elem_idx = 0;
                         elem_idx <
                         out_data[out_data_files[file_idx].varid[var_idx]].nelem;
                         elem_idx++) {
                        tmp_usiptr[ptr_idx++] =
                            (unsigned short int)out_data[out_data_files[file_idx
                                                         ].
                                                         varid[var_idx]].aggdata
                            [elem_idx];
                    }
                    fwrite(tmp_usiptr, sizeof(unsigned short int), ptr_idx,
                           out_data_files[file_idx].fh);
                }
                else if (out_data[out_data_files[file_idx].varid[var_idx]].type
                         ==
                         OUT_TYPE_INT) {
                    for (elem_idx = 0;
                         elem_idx <
                         out_data[out_data_files[file_idx].varid[var_idx]].nelem;
                         elem_idx++) {
                        tmp_iptr[ptr_idx++] =
                            (int)out_data[out_data_files[file_idx].varid[var_idx
                                          ]].
                            aggdata[elem_idx];
                    }
                    fwrite(tmp_iptr, sizeof(int), ptr_idx,
                           out_data_files[file_idx].fh);
                }
                else if (out_data[out_data_files[file_idx].varid[var_idx]].type
                         ==
                         OUT_TYPE_FLOAT) {
                    for (elem_idx = 0;
                         elem_idx <
                         out_data[out_data_files[file_idx].varid[var_idx]].nelem;
                         elem_idx++) {
                        tmp_fptr[ptr_idx++] =
                            (float)out_data[out_data_files[file_idx].varid[
                                                var_idx]]
                            .aggdata[elem_idx];
                    }
                    fwrite(tmp_fptr, sizeof(float), ptr_idx,
                           out_data_files[file_idx].fh);
                }
                else if (out_data[out_data_files[file_idx].varid[var_idx]].type
                         ==
                         OUT_TYPE_DOUBLE) {
                    for (elem_idx = 0;
                         elem_idx <
                         out_data[out_data_files[file_idx].varid[var_idx]].nelem;
                         elem_idx++) {
                        tmp_dptr[ptr_idx++] =
                            (double)out_data[out_data_files[file_idx].varid[
                                                 var_idx]
                            ].aggdata[elem_idx];
                    }
                    fwrite(tmp_dptr, sizeof(double), ptr_idx,
                           out_data_files[file_idx].fh);
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
    else { // ASCII
           // Loop over output files
        for (file_idx = 0; file_idx < options.Noutfiles; file_idx++) {
            if (!options.OUTPUT_FORCE) {
                // Write the date
                if (dt < HOURS_PER_DAY) {
                    // Write year, month, day, and hour
                    fprintf(out_data_files[file_idx].fh,
                            "%04i\t%02i\t%02i\t%02i\t",
                            dmy->year, dmy->month, dmy->day, dmy->hour);
                }
                else {
                    // Only write year, month, and day
                    fprintf(out_data_files[file_idx].fh, "%04i\t%02i\t%02i\t",
                            dmy->year, dmy->month, dmy->day);
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
                    fprintf(out_data_files[file_idx].fh,
                            out_data[out_data_files[file_idx].varid[var_idx]].format,
                            out_data[out_data_files[file_idx].varid[var_idx]].aggdata[
                                elem_idx]);
                }
            }
            fprintf(out_data_files[file_idx].fh, "\n");
        }
    }
}
