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

#include <vic_driver_cesm.h>

/******************************************************************************
 * @brief    Initialize x2l_data_struct.
 *****************************************************************************/
void
initialize_x2l_data(void)
{
    extern x2l_data_struct *x2l_vic;
    extern domain_struct    local_domain;

    size_t                  i;

    log_info("Setting all x2l fields to %f", SHR_CONST_SPVAL);

    for (i = 0; i < local_domain.ncells_active; i++) {
        x2l_vic[i].x2l_Sa_z = SHR_CONST_SPVAL;
        x2l_vic[i].x2l_Sa_u = SHR_CONST_SPVAL;
        x2l_vic[i].x2l_Sa_v = SHR_CONST_SPVAL;
        x2l_vic[i].x2l_Sa_ptem = SHR_CONST_SPVAL;
        x2l_vic[i].x2l_Sa_shum = SHR_CONST_SPVAL;
        x2l_vic[i].x2l_Sa_pbot = SHR_CONST_SPVAL;
        x2l_vic[i].x2l_Sa_tbot = SHR_CONST_SPVAL;
        x2l_vic[i].x2l_Faxa_lwdn = SHR_CONST_SPVAL;
        x2l_vic[i].x2l_Faxa_rainc = SHR_CONST_SPVAL;
        x2l_vic[i].x2l_Faxa_rainl = SHR_CONST_SPVAL;
        x2l_vic[i].x2l_Faxa_snowc = SHR_CONST_SPVAL;
        x2l_vic[i].x2l_Faxa_snowl = SHR_CONST_SPVAL;
        x2l_vic[i].x2l_Faxa_swndr = SHR_CONST_SPVAL;
        x2l_vic[i].x2l_Faxa_swvdr = SHR_CONST_SPVAL;
        x2l_vic[i].x2l_Faxa_swndf = SHR_CONST_SPVAL;
        x2l_vic[i].x2l_Faxa_swvdf = SHR_CONST_SPVAL;
        x2l_vic[i].x2l_Sa_co2prog = SHR_CONST_SPVAL;
        x2l_vic[i].x2l_Sa_co2diag = SHR_CONST_SPVAL;
        x2l_vic[i].x2l_Faxa_bcphidry = SHR_CONST_SPVAL;
        x2l_vic[i].x2l_Faxa_bcphodry = SHR_CONST_SPVAL;
        x2l_vic[i].x2l_Faxa_bcphiwet = SHR_CONST_SPVAL;
        x2l_vic[i].x2l_Faxa_ocphidry = SHR_CONST_SPVAL;
        x2l_vic[i].x2l_Faxa_ocphodry = SHR_CONST_SPVAL;
        x2l_vic[i].x2l_Faxa_ocphiwet = SHR_CONST_SPVAL;
        x2l_vic[i].x2l_Faxa_dstwet1 = SHR_CONST_SPVAL;
        x2l_vic[i].x2l_Faxa_dstwet2 = SHR_CONST_SPVAL;
        x2l_vic[i].x2l_Faxa_dstwet3 = SHR_CONST_SPVAL;
        x2l_vic[i].x2l_Faxa_dstwet4 = SHR_CONST_SPVAL;
        x2l_vic[i].x2l_Faxa_dstdry1 = SHR_CONST_SPVAL;
        x2l_vic[i].x2l_Faxa_dstdry2 = SHR_CONST_SPVAL;
        x2l_vic[i].x2l_Faxa_dstdry3 = SHR_CONST_SPVAL;
        x2l_vic[i].x2l_Faxa_dstdry4 = SHR_CONST_SPVAL;
        x2l_vic[i].x2l_Flrr_flood = SHR_CONST_SPVAL;
        x2l_vic[i].x2l_vars_set = false;
    }
}

/******************************************************************************
 * @brief    Initialize l2x_data_struct.
 *****************************************************************************/
void
initialize_l2x_data(void)
{
    extern l2x_data_struct *l2x_vic;
    extern domain_struct    local_domain;

    size_t                  i;

    log_info("Initializing l2x_data_struct");

    for (i = 0; i < local_domain.ncells_active; i++) {
        l2x_vic[i].l2x_Sl_t = SHR_CONST_SPVAL;
        l2x_vic[i].l2x_Sl_tref = SHR_CONST_SPVAL;
        l2x_vic[i].l2x_Sl_qref = SHR_CONST_SPVAL;
        l2x_vic[i].l2x_Sl_avsdr = SHR_CONST_SPVAL;
        l2x_vic[i].l2x_Sl_anidr = SHR_CONST_SPVAL;
        l2x_vic[i].l2x_Sl_avsdf = SHR_CONST_SPVAL;
        l2x_vic[i].l2x_Sl_anidf = SHR_CONST_SPVAL;
        l2x_vic[i].l2x_Sl_snowh = SHR_CONST_SPVAL;
        l2x_vic[i].l2x_Sl_u10 = SHR_CONST_SPVAL;
        l2x_vic[i].l2x_Sl_ddvel = SHR_CONST_SPVAL;
        l2x_vic[i].l2x_Sl_fv = SHR_CONST_SPVAL;
        l2x_vic[i].l2x_Sl_ram1 = SHR_CONST_SPVAL;
        l2x_vic[i].l2x_Sl_logz0 = SHR_CONST_SPVAL;
        l2x_vic[i].l2x_Fall_taux = SHR_CONST_SPVAL;
        l2x_vic[i].l2x_Fall_tauy = SHR_CONST_SPVAL;
        l2x_vic[i].l2x_Fall_lat = SHR_CONST_SPVAL;
        l2x_vic[i].l2x_Fall_sen = SHR_CONST_SPVAL;
        l2x_vic[i].l2x_Fall_lwup = SHR_CONST_SPVAL;
        l2x_vic[i].l2x_Fall_evap = SHR_CONST_SPVAL;
        l2x_vic[i].l2x_Fall_swnet = SHR_CONST_SPVAL;
        l2x_vic[i].l2x_Fall_fco2_lnd = SHR_CONST_SPVAL;
        l2x_vic[i].l2x_Fall_flxdst1 = SHR_CONST_SPVAL;
        l2x_vic[i].l2x_Fall_flxdst2 = SHR_CONST_SPVAL;
        l2x_vic[i].l2x_Fall_flxdst3 = SHR_CONST_SPVAL;
        l2x_vic[i].l2x_Fall_flxdst4 = SHR_CONST_SPVAL;
        l2x_vic[i].l2x_Fall_flxvoc = SHR_CONST_SPVAL;
        l2x_vic[i].l2x_Flrl_rofliq = SHR_CONST_SPVAL;
        l2x_vic[i].l2x_Flrl_rofice = SHR_CONST_SPVAL;
        l2x_vic[i].l2x_vars_set = true;
    }
}

/******************************************************************************
 * @brief      Initialize albedo values in l2x_data_struct.
 *****************************************************************************/
void
vic_initialize_albedo(void)
{
    extern l2x_data_struct *l2x_vic;
    extern domain_struct    local_domain;
    extern all_vars_struct *all_vars;

    size_t                  i;

    log_info("Initializing albedo values");

    for (i = 0; i < local_domain.ncells_active; i++) {
        l2x_vic[i].l2x_Sl_avsdr = all_vars[i].gridcell_avg.avg_albedo;
        l2x_vic[i].l2x_Sl_anidr = all_vars[i].gridcell_avg.avg_albedo;
        l2x_vic[i].l2x_Sl_avsdf = all_vars[i].gridcell_avg.avg_albedo;
        l2x_vic[i].l2x_Sl_anidf = all_vars[i].gridcell_avg.avg_albedo;
    }
}

/*****************************************************************************
 * @brief     Initialize temperature in l2x_data_struct.
 ****************************************************************************/
void
vic_initialize_temperature(void)
{
    extern l2x_data_struct *l2x_vic;
    extern domain_struct    local_domain;
    extern soil_con_struct *soil_con;

    size_t                  i;

    log_info("Initializing temperature");

    for (i = 0; i < local_domain.ncells_active; i++) {
        l2x_vic[i].l2x_Sl_t = soil_con[i].avg_temp + CONST_TKFRZ;
    }
}

/*****************************************************************************
 * @brief     Initialize upwelling longwave in l2x_data_struct.
 ****************************************************************************/
void
vic_initialize_lwup(void)
{
    extern l2x_data_struct  *l2x_vic;
    extern domain_struct     local_domain;
    extern parameters_struct param;

    size_t                   i;

    log_info("Initializing upwelling longwave");

    for (i = 0; i < local_domain.ncells_active; i++) {
        // adjust sign for CESM sign convention
        l2x_vic[i].l2x_Fall_lwup = -1 * param.EMISS_GRND * CONST_STEBOL * pow(
            l2x_vic[i].l2x_Sl_t, 4);
    }
}
