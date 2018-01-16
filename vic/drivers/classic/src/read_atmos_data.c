/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine reads in atmospheric data values from a binary/ascii file.
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
 * @brief    Read in atmospheric data values from a binary/ascii file.
 *****************************************************************************/
void
read_atmos_data(FILE               *infile,
                global_param_struct global_param,
                int                 file_num,
                int                 forceskip,
                double            **forcing_data,
                double           ***veg_hist_data)
{
    extern param_set_struct param_set;

    unsigned int            rec;
    unsigned int            skip_recs;
    unsigned int            i, j;
    int                     endian;
    unsigned int            Nfields;
    int                    *field_index;
    unsigned short int      ustmp;
    signed short            stmp;
    char                    str[MAXSTRING + 1];
    unsigned short int      Identifier[4];
    int                     Nbytes;

    Nfields = param_set.N_TYPES[file_num];
    field_index = param_set.FORCE_INDEX[file_num];

    /** locate starting record **/

    /* if ascii then the following refers to the number of lines to skip,
       if binary the following needs multiplying by the number of input fields */
    skip_recs = (unsigned int) ((global_param.dt * forceskip)) /
                param_set.FORCE_DT[file_num];
    if ((((global_param.dt < SEC_PER_DAY &&
           (unsigned int) (param_set.FORCE_DT[file_num] * forceskip) %
           (unsigned int) global_param.dt) > 0)) ||
        (global_param.dt == SEC_PER_DAY &&
         ((unsigned int) global_param.dt %
          (unsigned int) param_set.FORCE_DT[file_num] >
          0))) {
        log_err("Currently unable to handle a model starting date that does "
                "not correspond to a line in the forcing file.");
    }

    /** Error checking - Model can be run at any time step using daily forcing
        data, but if sub-daily data is used, the model must be run at the
        same time step as the data.  That way aggregation and disaggragation
        techniques are left to the user. **/
    if (param_set.FORCE_DT[file_num] < SEC_PER_DAY &&
        global_param.dt != param_set.FORCE_DT[file_num]) {
        log_err("When forcing the model with sub-daily data, the model must be "
                "run at the same time step as the forcing data.  Currently the "
                "model time step is %f seconds, while forcing file %i has a "
                "time step of %f seconds.", global_param.dt, file_num,
                param_set.FORCE_DT[file_num]);
    }

    if (infile == NULL) {
        log_info("NULL file");
    }

    /***************************
       Read BINARY Forcing Data
    ***************************/

    if (param_set.FORCE_FORMAT[file_num] == BINARY) {
        /** test whether the machine is little-endian or big-endian **/
        i = 1;
        if (*(char *)&i == 1) {
            endian = LITTLE;
        }
        else {
            endian = BIG;
        }

        // Check for presence of a header, & skip over it if appropriate.
        // A VIC header will start with 4 instances of the identifier,
        // followed by number of bytes in the header (Nbytes).
        // Nbytes is assumed to be the byte offset at which the data records start.
        fseek(infile, 0, SEEK_SET);
        if (feof(infile)) {
            log_err("No data in the forcing file.");
        }
        for (i = 0; i < 4; i++) {
            fread(&ustmp, sizeof(unsigned short int), 1, infile);
            if (endian != param_set.FORCE_ENDIAN[file_num]) {
                ustmp = ((ustmp & 0xFF) << 8) | ((ustmp >> 8) & 0xFF);
            }
            Identifier[i] = ustmp;
        }
        if (Identifier[0] != 0xFFFF || Identifier[1] != 0xFFFF ||
            Identifier[2] != 0xFFFF || Identifier[3] != 0xFFFF) {
            Nbytes = 0;
        }
        else {
            fread(&ustmp, sizeof(unsigned short int), 1, infile);
            if (endian != param_set.FORCE_ENDIAN[file_num]) {
                ustmp = ((ustmp & 0xFF) << 8) | ((ustmp >> 8) & 0xFF);
            }
            Nbytes = (int) ustmp;
        }
        fseek(infile, Nbytes, SEEK_SET);


        /** if forcing file starts before the model simulation,
            skip over its starting records **/
        fseek(infile, skip_recs * Nfields * sizeof(short int), SEEK_CUR);
        if (feof(infile)) {
            log_err("No data for the specified time period in the forcing "
                    "file.");
        }

        /** Read BINARY forcing data **/
        rec = 0;

        while (!feof(infile) && (rec * param_set.FORCE_DT[file_num] <
                                 global_param.nrecs * global_param.dt)) {
            for (i = 0; i < Nfields; i++) {
                if (field_index[i] != ALBEDO && field_index[i] != LAI &&
                    field_index[i] != FCANOPY) {
                    if (param_set.TYPE[field_index[i]].SIGNED) {
                        fread(&stmp, sizeof(short int), 1, infile);
                        if (endian != param_set.FORCE_ENDIAN[file_num]) {
                            stmp = ((stmp & 0xFF) << 8) | ((stmp >> 8) & 0xFF);
                        }
                        forcing_data[field_index[i]][rec] =
                            (double) stmp /
                            param_set.TYPE[field_index[i]].multiplier;
                    }
                    else {
                        fread(&ustmp, sizeof(unsigned short int), 1, infile);
                        if (endian != param_set.FORCE_ENDIAN[file_num]) {
                            ustmp =
                                ((ustmp & 0xFF) << 8) | ((ustmp >> 8) & 0xFF);
                        }
                        forcing_data[field_index[i]][rec] =
                            (double) ustmp /
                            param_set.TYPE[field_index[i]].multiplier;
                    }
                }
                else {
                    for (j = 0; j < param_set.TYPE[field_index[i]].N_ELEM;
                         j++) {
                        if (param_set.TYPE[field_index[i]].SIGNED) {
                            fread(&stmp, sizeof(short int), 1, infile);
                            if (endian != param_set.FORCE_ENDIAN[file_num]) {
                                stmp =
                                    ((stmp & 0xFF) << 8) | ((stmp >> 8) & 0xFF);
                            }
                            veg_hist_data[field_index[i]][j][rec] =
                                (double) stmp /
                                param_set.TYPE[field_index[i]].multiplier;
                        }
                        else {
                            fread(&ustmp, sizeof(unsigned short int), 1,
                                  infile);
                            if (endian != param_set.FORCE_ENDIAN[file_num]) {
                                ustmp =
                                    ((ustmp &
                                      0xFF) << 8) | ((ustmp >> 8) & 0xFF);
                            }
                            veg_hist_data[field_index[i]][j][rec] =
                                (double) ustmp /
                                param_set.TYPE[field_index[i]].multiplier;
                        }
                    }
                }
            }

            rec++;
        }
    }

    /**************************
       Read ASCII Forcing Data
    **************************/

    else {
        // No need to skip over a header here, since ascii file headers are skipped
        // in open_file().  However, if we wanted to read information from the header,
        // we'd want to do it here, after rewinding to the beginning of the file (or
        // moving the code that deals with headers from open_file() to this function
        // and to any other functions that read the files, so that those functions could
        // also read the headers if necessary).

        /* skip to the beginning of the required met data */
        for (i = 0; i < skip_recs; i++) {
            if (fgets(str, MAXSTRING, infile) == NULL) {
                log_err("No data for the specified time period in the forcing "
                        "file.");
            }
        }

        /* read forcing data */
        rec = 0;

        while (!feof(infile) && (rec * param_set.FORCE_DT[file_num] <
                                 global_param.nrecs * global_param.dt)) {
            for (i = 0; i < Nfields; i++) {
                if (field_index[i] != ALBEDO && field_index[i] != LAI &&
                    field_index[i] != FCANOPY) {
                    fscanf(infile, "%lf", &forcing_data[field_index[i]][rec]);
                }
                else {
                    for (j = 0; j < param_set.TYPE[field_index[i]].N_ELEM;
                         j++) {
                        fscanf(infile, "%lf",
                               &veg_hist_data[field_index[i]][j][rec]);
                    }
                }
            }
            fgets(str, MAXSTRING, infile);
            rec++;
        }
    }

    if (rec * param_set.FORCE_DT[file_num] <
        global_param.nrecs * global_param.dt) {
        log_err("Not enough records in forcing file %i (%u * %f = %f) to run "
                "the number of records defined in the global file "
                "(%zu * %f = %f).  Check forcing file time step, and global "
                "file", file_num + 1, rec, param_set.FORCE_DT[file_num],
                rec * param_set.FORCE_DT[file_num], global_param.nrecs,
                global_param.dt,
                global_param.nrecs * global_param.dt);
    }
}
