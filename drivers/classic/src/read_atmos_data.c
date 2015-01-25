/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine reads in atmospheric data values from a binary/ascii file.
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
 * @brief    Read in atmospheric data values from a binary/ascii file.
 *****************************************************************************/
void
read_atmos_data(FILE               *infile,
                global_param_struct global_param,
                int                 file_num,
                double            **forcing_data,
                double           ***veg_hist_data)
{
    extern param_set_struct param_set;

    size_t                  rec;
    size_t                  force_nrecs;
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

    // Calculate the number of forcing records that need to be read
    force_nrecs = global_param.nrecs * param_set.force_steps_per_day[file_num] /
                  global_param.model_steps_per_day;

    /** Error checking - Model can be run at any time step using daily forcing
        data, but if sub-daily data is used, the model must be run at the
        same time step as the data.  That way aggregation and disaggragation
        techniques are left to the user. **/
    if (param_set.force_steps_per_day[file_num] > 1 &&
        global_param.model_steps_per_day !=
        param_set.force_steps_per_day[file_num]) {
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
            log_err("No data in the forcing file.  Model stopping...");
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

        /** Read BINARY forcing data **/
        rec = 0;

        while (!feof(infile) && (rec < force_nrecs)) {
            for (i = 0; i < Nfields; i++) {
                if (field_index[i] != ALBEDO && field_index[i] != LAI_IN &&
                    field_index[i] != VEGCOVER) {
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
        /* read forcing data */
        rec = 0;

        while (!feof(infile) && (rec < force_nrecs)) {
            for (i = 0; i < Nfields; i++) {
                if (field_index[i] != ALBEDO && field_index[i] != LAI_IN &&
                    field_index[i] != VEGCOVER) {
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

    if (rec < force_nrecs) {
        log_err("Not enough records in forcing file %i (%zu * %f = %f) to run "
                "the number of records defined in the global file "
                "(%zu * %f = %f).  Check forcing file time step, and global "
                "file", file_num + 1, rec, param_set.FORCE_DT[file_num],
                rec * param_set.FORCE_DT[file_num], global_param.nrecs,
                global_param.dt,
                global_param.nrecs * global_param.dt);
    }
}
