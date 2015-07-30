/******************************************************************************
 * @section DESCRIPTION
 *
 * Initialize routing model parameters
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
#include <vic_driver_image.h>

#include <rout.h>

/******************************************************************************
 * @brief    Initialize routing model parameters
 *****************************************************************************/
void
rout_init(void)
{
    extern rout_struct         rout;
    int                       *ivar = NULL;
    double                    *dvar = NULL;

    size_t                     i, j;
    size_t                     i1start;
    size_t                     d3count[3];
    size_t                     d3start[3];

    i1start = 0;
    
    d3start[0] = 0;
    d3start[1] = 0;
    d3start[2] = 0;
    d3count[0] = rout.rout_param.iSubsetLength;
    d3count[1] = rout.rout_param.iSources;
    d3count[2] = 1; // tracers dimension


    // allocate memory for variables to be read
    ivar = (int *) malloc(sizeof(int) * rout.rout_param.iSources);
    if (ivar == NULL) {
        log_err("Memory allocation error in vic_init().");
    }
    
    // allocate memory for variables to be read
    dvar = (double *) malloc(rout.rout_param.iSubsetLength * rout.rout_param.iSources *
                             sizeof(double));
    if (dvar == NULL) {
        log_err("Memory allocation error in vic_init().");
    }
    
    // The Ring
    for (j=0;j<rout.rout_param.iOutlets;j++) {
      for (i=0;i<rout.rout_param.iSubsetLength;i++) {
        rout.ring[j*rout.rout_param.iSubsetLength + i] = 0.0;
      }
    }

    // source2outlet_ind: source to outlet index mapping
    get_nc_field_int(rout.param_filename,
                                "source2outlet_ind",
                                &i1start, &rout.rout_param.iSources, ivar);
    for (i = 0; i < rout.rout_param.iSources; i++) {
        rout.rout_param.source2outlet_ind[i] = (int) ivar[i];
    }

    // source_time_offset: Number of leading timesteps ommited
    get_nc_field_int(rout.param_filename,
                                "source_time_offset",
                                &i1start, &rout.rout_param.iSources, ivar);
    for (i = 0; i < rout.rout_param.iSources; i++) {
        rout.rout_param.source_time_offset[i] = (int) ivar[i];
    }

    // source_x_ind: x grid coordinate of source grid cell
    get_nc_field_int(rout.param_filename,
                                "source_x_ind",
                                &i1start, &rout.rout_param.iSources, ivar);
    for (i = 0; i < rout.rout_param.iSources; i++) {
        rout.rout_param.source_x_ind[i] = (int) ivar[i];
    }

    // source_y_ind: y grid coordinate of source grid cell
    get_nc_field_int(rout.param_filename,
                                "source_y_ind",
                                &i1start, &rout.rout_param.iSources, ivar);
    for (i = 0; i < rout.rout_param.iSources; i++) {
        rout.rout_param.source_y_ind[i] = (int) ivar[i];
    }

    // Unit Hydrograph:
    get_nc_field_double(rout.param_filename,
                                "unit_hydrograph",
                                d3start, d3count, dvar);
    for (i = 0; i < (rout.rout_param.iSubsetLength * rout.rout_param.iSources); i++) {
        rout.rout_param.unit_hydrograph[i] = (double) dvar[i];
    }

    // cleanup
    free(ivar);
    free(dvar);
}
