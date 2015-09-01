/******************************************************************************
 * @section DESCRIPTION
 *
 * Initilization library.
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
#include <vic_driver_cesm.h>

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
    soil_con->gridcel = -1;

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
        soil_con->AboveTreeLine[i] = 0;
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
    veg_con->Cv_sum = 0.;
    veg_con->veg_class = -1; // -1 to force a crash if inappropriate
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
 * @brief    Initialize x2l_data_struct.
 *****************************************************************************/
void
initialize_x2l_data()
{
    extern x2l_data_struct *x2l_vic;
    extern domain_struct    local_domain;

    size_t                  i;

    log_info("Setting all x2l fields to %f", MISSING);

    for (i = 0; i < local_domain.ncells; i++) {
        x2l_vic[i].x2l_Sa_z = MISSING;
        x2l_vic[i].x2l_Sa_u = MISSING;
        x2l_vic[i].x2l_Sa_v = MISSING;
        x2l_vic[i].x2l_Sa_ptem = MISSING;
        x2l_vic[i].x2l_Sa_shum = MISSING;
        x2l_vic[i].x2l_Sa_pbot = MISSING;
        x2l_vic[i].x2l_Sa_tbot = MISSING;
        x2l_vic[i].x2l_Faxa_lwdn = MISSING;
        x2l_vic[i].x2l_Faxa_rainc = MISSING;
        x2l_vic[i].x2l_Faxa_rainl = MISSING;
        x2l_vic[i].x2l_Faxa_snowc = MISSING;
        x2l_vic[i].x2l_Faxa_snowl = MISSING;
        x2l_vic[i].x2l_Faxa_swndr = MISSING;
        x2l_vic[i].x2l_Faxa_swvdr = MISSING;
        x2l_vic[i].x2l_Faxa_swndf = MISSING;
        x2l_vic[i].x2l_Faxa_swvdf = MISSING;
        x2l_vic[i].x2l_Sa_co2prog = MISSING;
        x2l_vic[i].x2l_Sa_co2diag = MISSING;
        x2l_vic[i].x2l_Faxa_bcphidry = MISSING;
        x2l_vic[i].x2l_Faxa_bcphodry = MISSING;
        x2l_vic[i].x2l_Faxa_bcphiwet = MISSING;
        x2l_vic[i].x2l_Faxa_ocphidry = MISSING;
        x2l_vic[i].x2l_Faxa_ocphodry = MISSING;
        x2l_vic[i].x2l_Faxa_ocphiwet = MISSING;
        x2l_vic[i].x2l_Faxa_dstwet1 = MISSING;
        x2l_vic[i].x2l_Faxa_dstwet2 = MISSING;
        x2l_vic[i].x2l_Faxa_dstwet3 = MISSING;
        x2l_vic[i].x2l_Faxa_dstwet4 = MISSING;
        x2l_vic[i].x2l_Faxa_dstdry1 = MISSING;
        x2l_vic[i].x2l_Faxa_dstdry2 = MISSING;
        x2l_vic[i].x2l_Faxa_dstdry3 = MISSING;
        x2l_vic[i].x2l_Faxa_dstdry4 = MISSING;
        x2l_vic[i].x2l_Flrr_flood = MISSING;
        x2l_vic[i].x2l_vars_set = false;
    }
}

/******************************************************************************
 * @brief    Initialize l2x_data_struct.
 *****************************************************************************/
void
initialize_l2x_data()
{
    extern l2x_data_struct *l2x_vic;
    extern domain_struct    local_domain;

    size_t                  i;

    log_info("Setting all l2x fields to %f", MISSING);

    for (i = 0; i < local_domain.ncells; i++) {
        l2x_vic[i].l2x_Sl_t = MISSING;
        l2x_vic[i].l2x_Sl_tref = MISSING;
        l2x_vic[i].l2x_Sl_qref = MISSING;
        l2x_vic[i].l2x_Sl_avsdr = MISSING;
        l2x_vic[i].l2x_Sl_anidr = MISSING;
        l2x_vic[i].l2x_Sl_avsdf = MISSING;
        l2x_vic[i].l2x_Sl_anidf = MISSING;
        l2x_vic[i].l2x_Sl_snowh = MISSING;
        l2x_vic[i].l2x_Sl_u10 = MISSING;
        l2x_vic[i].l2x_Sl_ddvel = MISSING;
        l2x_vic[i].l2x_Sl_fv = MISSING;
        l2x_vic[i].l2x_Sl_ram1 = MISSING;
        l2x_vic[i].l2x_Sl_soilw = MISSING;
        l2x_vic[i].l2x_Sl_logz0 = MISSING;
        l2x_vic[i].l2x_Fall_taux = MISSING;
        l2x_vic[i].l2x_Fall_tauy = MISSING;
        l2x_vic[i].l2x_Fall_lat = MISSING;
        l2x_vic[i].l2x_Fall_sen = MISSING;
        l2x_vic[i].l2x_Fall_lwup = MISSING;
        l2x_vic[i].l2x_Fall_evap = MISSING;
        l2x_vic[i].l2x_Fall_swnet = MISSING;
        l2x_vic[i].l2x_Fall_fco2_lnd = MISSING;
        l2x_vic[i].l2x_Fall_flxdst1 = MISSING;
        l2x_vic[i].l2x_Fall_flxdst2 = MISSING;
        l2x_vic[i].l2x_Fall_flxdst3 = MISSING;
        l2x_vic[i].l2x_Fall_flxdst4 = MISSING;
        l2x_vic[i].l2x_Fall_flxvoc = MISSING;
        l2x_vic[i].l2x_Flrl_rofliq = MISSING;
        l2x_vic[i].l2x_Flrl_rofice = MISSING;
        l2x_vic[i].l2x_vars_set = false;
    }
}
