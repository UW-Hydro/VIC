/******************************************************************************
 * @section DESCRIPTION
 *
 * This subroutine controls the order and number of forcing variables read from
 * the forcing data files.  Two forcing files are allowed, variables, time step
 * and file format must be defined in the global control file.
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
 * @brief    Control the order and number of forcing variables read from the
 *           forcing data files.
 *****************************************************************************/
double **
read_forcing_data(FILE              **infile,
                  global_param_struct global_param,
                  double          ****veg_hist_data)
{
    extern param_set_struct param_set;
    extern size_t           NF;

    size_t                  i, j;
    double                **forcing_data;

    /** Allocate data arrays for input forcing data **/
    forcing_data = (double **)calloc(N_FORCING_TYPES, sizeof(double*));
    (*veg_hist_data) = (double ***)calloc(N_FORCING_TYPES, sizeof(double**));
    for (i = 0; i < N_FORCING_TYPES; i++) {
        if (param_set.TYPE[i].SUPPLIED) {
            if (i != ALBEDO && i != LAI_IN && i != VEGCOVER) {
                forcing_data[i] = (double *)calloc((global_param.nrecs * NF),
                                                   sizeof(double));
            }
            else {
                (*veg_hist_data)[i] = (double **)calloc(
                    param_set.TYPE[i].N_ELEM, sizeof(double*));
                for (j = 0; j < param_set.TYPE[i].N_ELEM; j++) {
                    (*veg_hist_data)[i][j] =
                        (double *)calloc((global_param.nrecs * NF),
                                         sizeof(double));
                }
            }
        }
    }

    /** Read First Forcing Data File **/
    if (param_set.FORCE_DT[0] > 0) {
        read_atmos_data(infile[0], global_param, 0, global_param.forceskip[0],
                        forcing_data, (*veg_hist_data));
    }
    else {
        log_err("File time step must be defined for at least the first "
                "forcing file (FILE_DT).\n");
    }

    /** Read Second Forcing Data File **/
    if (param_set.FORCE_DT[1] > 0) {
        read_atmos_data(infile[1], global_param, 1, global_param.forceskip[1],
                        forcing_data, (*veg_hist_data));
    }

    return(forcing_data);
}
