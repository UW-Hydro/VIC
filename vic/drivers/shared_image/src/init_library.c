/******************************************************************************
 * @section DESCRIPTION
 *
 * Initilization library.
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

#include <vic_driver_shared_image.h>

/******************************************************************************
 * @brief    Initialize soil con sructure.
 *****************************************************************************/
void
initialize_soil_con(soil_con_struct *soil_con)
{
    extern option_struct options;
    size_t               i;
    size_t               j;


    soil_con->FS_ACTIVE = 0;
    soil_con->gridcel = MISSING_USI;

    soil_con->AlbedoPar = 0.;
    soil_con->elevation = 0.;
    soil_con->lat = 0.;
    soil_con->lng = 0.;
    soil_con->time_zone_lng = 0.;

    soil_con->annual_prec = 0.;
    soil_con->aspect = 0.;
    soil_con->avg_temp = 0.;
    soil_con->avgJulyAirTemp = 0.;
    soil_con->b_infilt = 0.;
    soil_con->c = 0.;
    soil_con->cell_area = 0.;
    soil_con->dp = 0.;
    soil_con->Ds = 0.;
    soil_con->Dsmax = 0.;
    soil_con->ehoriz = 0.;
    soil_con->frost_slope = 0.;
    soil_con->max_infil = 0.;
    soil_con->max_snow_distrib_slope = 0.;
    soil_con->rough = 0.;
    soil_con->slope = 0.;
    soil_con->snow_rough = 0.;
    soil_con->whoriz = 0.;
    soil_con->Ws = 0.;

    for (i = 0; i < MAX_LAYERS; i++) {
        soil_con->bubble[i] = 0.;
        soil_con->bulk_density[i] = 0.;
        soil_con->bulk_dens_min[i] = 0.;
        soil_con->bulk_dens_org[i] = 0.;
        soil_con->depth[i] = 0.;
        soil_con->expt[i] = 0.;
        soil_con->init_moist[i] = 0.;
        soil_con->Ksat[i] = 0.;
        soil_con->max_moist[i] = 0.;
        soil_con->phi_s[i] = 0.;
        soil_con->porosity[i] = 0.;
        soil_con->quartz[i] = 0.;
        soil_con->organic[i] = 0.;
        soil_con->resid_moist[i] = 0.;
        soil_con->soil_density[i] = 0.;
        soil_con->soil_dens_min[i] = 0.;
        soil_con->soil_dens_org[i] = 0.;
        soil_con->Wcr[i] = 0.;
        soil_con->Wpwp[i] = 0.;
    }

    for (i = 0; i < MAX_NODES; i++) {
        soil_con->alpha[i] = 0.;
        soil_con->beta[i] = 0.;
        soil_con->bubble_node[i] = 0.;
        soil_con->dz_node[i] = 0.;
        soil_con->Zsum_node[i] = 0.;
        soil_con->expt_node[i] = 0.;
        soil_con->gamma[i] = 0.;
        soil_con->max_moist_node[i] = 0.;
    }

    for (i = 0; i < MAX_FROST_AREAS; i++) {
        soil_con->frost_fract[i] = 0.;
    }

    for (i = 0; i < options.SNOW_BAND; i++) {
        soil_con->AboveTreeLine[i] = false;
        soil_con->BandElev[i] = 0.;
        soil_con->AreaFract[i] = 1.;
        soil_con->Pfactor[i] = 0.;
        soil_con->Tfactor[i] = 0.;
    }

    for (i = 0; i < MAX_LAYERS + 2; i++) {
        for (j = 0; j < MAX_ZWTVMOIST; j++) {
            soil_con->zwtvmoist_zwt[i][j] = 0.;
            soil_con->zwtvmoist_moist[i][j] = 0.;
        }
    }
}

/******************************************************************************
 * @brief    Initialize veg con sructure.
 *****************************************************************************/
void
initialize_veg_con(veg_con_struct *veg_con)
{
    extern option_struct options;
    size_t               i;

    veg_con->Cv = 0.;
    veg_con->veg_class = NODATA_VEG; // -1 to force a crash if inappropriate
    veg_con->vegetat_type_num = 0.;
    veg_con->sigma_slope = 0.;
    veg_con->lag_one = 0.;
    veg_con->fetch = 0.;
    veg_con->LAKE = 0;
    for (i = 0; i < MAX_LAYERS; i++) {
        veg_con->root[i] = 0.;
    }
    for (i = 0; i < options.ROOT_ZONES; i++) {
        veg_con->zone_depth[i] = 0.;
        veg_con->zone_fract[i] = 0.;
    }
    if (options.CARBON) {
        for (i = 0; i < options.Ncanopy; i++) {
            veg_con->CanopLayerBnd[i] = 0.;
        }
    }
}

/******************************************************************************
 * @brief    Initialize domain info stucture
 *****************************************************************************/
void
initialize_domain_info(domain_info_struct *info)
{
    strcpy(info->lat_var, "MISSING");
    strcpy(info->lon_var, "MISSING");
    strcpy(info->mask_var, "MISSING");
    strcpy(info->area_var, "MISSING");
    strcpy(info->frac_var, "MISSING");
    strcpy(info->y_dim, "MISSING");
    strcpy(info->x_dim, "MISSING");
    info->n_coord_dims = 0;
}

/******************************************************************************
 * @brief    Initialize global structures
 *****************************************************************************/
void
initialize_global_structures(void)
{
    extern domain_struct global_domain;
    extern domain_struct local_domain;
    extern int           mpi_rank;

    initialize_domain_info(&local_domain.info);
    if (mpi_rank == VIC_MPI_ROOT) {
        initialize_options();
        initialize_global();
        initialize_parameters();
        initialize_filenames();
        initialize_domain_info(&global_domain.info);
        initialize_domain(&global_domain);
    }
}
