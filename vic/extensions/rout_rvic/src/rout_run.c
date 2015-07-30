/******************************************************************************
 * @section DESCRIPTION
 *
 * Run routing over the domain.
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
#include <rout.h>
#include <assert.h>
#include <vic_driver_image.h>

void rout_run(void)
{
    log_info("In Routing Lohmann Model");

    extern rout_struct         rout;
    extern domain_struct       global_domain;  
    extern out_data_struct   **out_data;
    size_t i;

    // even worse stuff here 
    for (i = 0; i < global_domain.ncells; i++) {
       rout.rout_param.aggrunin[i]= out_data[i][OUT_RUNOFF].data[0] + out_data[i][OUT_BASEFLOW].data[0];
    }

    //TODO: Check if nOutlets in ROUTING is the same (same lat lons and same length) as VIC
    // completely useless double stuff here    
 
    log_info("Tolkien here");
    size_t  j;
    size_t  ni = rout.rout_param.iSubsetLength;
    size_t  nj = rout.rout_param.iOutlets;
    
//    printf("\nThe Ring...\n");
//    for (j=0;j<nj;j++) {
//      for (i=0;i<ni;i++) {
//        rout.ring[j*ni + i] = j*ni + i;
//      }
//    }
//    for (j=0;j<nj;j++) {
//      for (i=0;i<ni;i++) {
//        printf("%7.0f",rout.ring[j*ni + i]);
//      }
//      printf("\n");
//    }

    // Zero out current ring
    // in python: (from variables.py) self.ring[tracer][0, :] = 0.
    printf("\nThe Ring (ZERO)...\n");
    for (j=0;j<nj;j++) {
        rout.ring[j*ni] = 0.0;
    }
    for (j=0;j<nj;j++) {
      for (i=0;i<ni;i++) {
        printf("%7.5f",rout.ring[j*ni + i]);
      }
      printf("\n");
    }

    // Equivalent to Fortran 90 cshift function
    // in python: (from variables.py) self.ring[tracer] = np.roll(self.ring[tracer], -1, axis=0)
    for (j=0;j<nj;j++) {
        cshift(rout.ring, rout.rout_param.iSubsetLength, j);
    }
    
    printf("\nThe Ring, c-shifted...\n");
    for (j=0;j<nj;j++) {
      for (i=0;i<ni;i++) {
        printf("%7.5f",rout.ring[j*ni + i]);
      }
      printf("\n");
    }

    convolve(rout.rout_param.iSources,               /*scalar - number of sources*/
                rout.rout_param.iOutlets,            /*scalar - length of subset*/
                rout.rout_param.iSubsetLength,       /*scalar - length of subset*/
                global_domain.n_nx,
                rout.rout_param.source2outlet_ind,   /*1d array - source to outlet mapping*/
                rout.rout_param.source_y_ind,        /*1d array - source y location*/
                rout.rout_param.source_x_ind,        /*1d array - source x location*/
                rout.rout_param.source_time_offset,  /*1d array - source time offset*/
                rout.rout_param.unit_hydrograph,     /*2d array[times][sources] - unit hydrographs*/
                rout.rout_param.aggrunin,            /*2d array[ysize][xsize] - vic runoff flux*/
                rout.ring);                          /*2d array[times][outlets] - convolution ring*/
}

