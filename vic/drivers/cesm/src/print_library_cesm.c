/******************************************************************************
 * @section DESCRIPTION
 *
 * Print library.
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
 * @brief    Print vic_clock structure.  See also cesm_print_library_f.F90
 *****************************************************************************/
void
print_vic_clock(vic_clock *vclock)
{
    fprintf(LOG_DEST, "vic_clock            :\n");
    fprintf(LOG_DEST, "\ttimestep           : %d\n", vclock->timestep);
    fprintf(LOG_DEST, "\tcurrent_year       : %d\n", vclock->current_year);
    fprintf(LOG_DEST, "\tcurrent_month      : %d\n", vclock->current_month);
    fprintf(LOG_DEST, "\tcurrent_day        : %d\n", vclock->current_day);
    fprintf(LOG_DEST, "\tcurrent_dayseconds : %d\n",
            vclock->current_dayseconds);
    fprintf(LOG_DEST, "\tstate_flag         : %d\n", vclock->state_flag);
    fprintf(LOG_DEST, "\tstop_flag          : %d\n", vclock->stop_flag);
    fprintf(LOG_DEST, "\tcalendar           : %s\n", trimstr(vclock->calendar));
}

/******************************************************************************
 * @brief    Print case_metadata structure.  See also cesm_print_library_f.F90
 *****************************************************************************/
void
print_case_metadata(case_metadata *cmeta)
{
    fprintf(LOG_DEST, "case_metadata   :\n");
    fprintf(LOG_DEST, "\tcaseid        : %s\n", trimstr(cmeta->caseid));
    fprintf(LOG_DEST, "\tcasedesc      : %s\n", trimstr(cmeta->casedesc));
    fprintf(LOG_DEST, "\tstarttype     : %s\n", trimstr(cmeta->starttype));
    fprintf(LOG_DEST, "\tmodel_version : %s\n", trimstr(cmeta->model_version));
    fprintf(LOG_DEST, "\thostname      : %s\n", trimstr(cmeta->hostname));
    fprintf(LOG_DEST, "\tusername      : %s\n", trimstr(cmeta->username));
}

/******************************************************************************
 * @brief    Print x2l_data_struct.  See also cesm_print_library_f.F90
 *****************************************************************************/
void
print_x2l_data(x2l_data_struct *x2l)
{
    fprintf(LOG_DEST, "x2l_data            :\n");
    fprintf(LOG_DEST, "\tx2l_Sa_z          : %f\n", x2l->x2l_Sa_z);
    fprintf(LOG_DEST, "\tx2l_Sa_u          : %f\n", x2l->x2l_Sa_u);
    fprintf(LOG_DEST, "\tx2l_Sa_v          : %f\n", x2l->x2l_Sa_v);
    fprintf(LOG_DEST, "\tx2l_Sa_ptem       : %f\n", x2l->x2l_Sa_ptem);
    fprintf(LOG_DEST, "\tx2l_Sa_shum       : %f\n", x2l->x2l_Sa_shum);
    fprintf(LOG_DEST, "\tx2l_Sa_pbot       : %f\n", x2l->x2l_Sa_pbot);
    fprintf(LOG_DEST, "\tx2l_Sa_tbot       : %f\n", x2l->x2l_Sa_tbot);
    fprintf(LOG_DEST, "\tx2l_Faxa_lwdn     : %f\n", x2l->x2l_Faxa_lwdn);
    fprintf(LOG_DEST, "\tx2l_Faxa_rainc    : %f\n", x2l->x2l_Faxa_rainc);
    fprintf(LOG_DEST, "\tx2l_Faxa_rainl    : %f\n", x2l->x2l_Faxa_rainl);
    fprintf(LOG_DEST, "\tx2l_Faxa_snowc    : %f\n", x2l->x2l_Faxa_snowc);
    fprintf(LOG_DEST, "\tx2l_Faxa_snowl    : %f\n", x2l->x2l_Faxa_snowl);
    fprintf(LOG_DEST, "\tx2l_Faxa_swndr    : %f\n", x2l->x2l_Faxa_swndr);
    fprintf(LOG_DEST, "\tx2l_Faxa_swvdr    : %f\n", x2l->x2l_Faxa_swvdr);
    fprintf(LOG_DEST, "\tx2l_Faxa_swndf    : %f\n", x2l->x2l_Faxa_swndf);
    fprintf(LOG_DEST, "\tx2l_Faxa_swvdf    : %f\n", x2l->x2l_Faxa_swvdf);
    fprintf(LOG_DEST, "\tx2l_Sa_co2prog    : %f\n", x2l->x2l_Sa_co2prog);
    fprintf(LOG_DEST, "\tx2l_Sa_co2diag    : %f\n", x2l->x2l_Sa_co2diag);
    fprintf(LOG_DEST, "\tx2l_Faxa_bcphidry : %f\n", x2l->x2l_Faxa_bcphidry);
    fprintf(LOG_DEST, "\tx2l_Faxa_bcphodry : %f\n", x2l->x2l_Faxa_bcphodry);
    fprintf(LOG_DEST, "\tx2l_Faxa_bcphiwet : %f\n", x2l->x2l_Faxa_bcphiwet);
    fprintf(LOG_DEST, "\tx2l_Faxa_ocphidry : %f\n", x2l->x2l_Faxa_ocphidry);
    fprintf(LOG_DEST, "\tx2l_Faxa_ocphodry : %f\n", x2l->x2l_Faxa_ocphodry);
    fprintf(LOG_DEST, "\tx2l_Faxa_ocphiwet : %f\n", x2l->x2l_Faxa_ocphiwet);
    fprintf(LOG_DEST, "\tx2l_Faxa_dstwet1  : %f\n", x2l->x2l_Faxa_dstwet1);
    fprintf(LOG_DEST, "\tx2l_Faxa_dstwet2  : %f\n", x2l->x2l_Faxa_dstwet2);
    fprintf(LOG_DEST, "\tx2l_Faxa_dstwet3  : %f\n", x2l->x2l_Faxa_dstwet3);
    fprintf(LOG_DEST, "\tx2l_Faxa_dstwet4  : %f\n", x2l->x2l_Faxa_dstwet4);
    fprintf(LOG_DEST, "\tx2l_Faxa_dstdry1  : %f\n", x2l->x2l_Faxa_dstdry1);
    fprintf(LOG_DEST, "\tx2l_Faxa_dstdry2  : %f\n", x2l->x2l_Faxa_dstdry2);
    fprintf(LOG_DEST, "\tx2l_Faxa_dstdry3  : %f\n", x2l->x2l_Faxa_dstdry3);
    fprintf(LOG_DEST, "\tx2l_Faxa_dstdry4  : %f\n", x2l->x2l_Faxa_dstdry4);
    fprintf(LOG_DEST, "\tx2l_Flrr_flood    : %f\n", x2l->x2l_Flrr_flood);
}

/******************************************************************************
 * @brief    Print l2x_data_struct.  See also cesm_print_library_f.F90
 *****************************************************************************/
void
print_l2x_data(l2x_data_struct *l2x)
{
    fprintf(LOG_DEST, "l2x_data            :\n");
    fprintf(LOG_DEST, "\tl2x_Sl_t          : %f\n", l2x->l2x_Sl_t);
    fprintf(LOG_DEST, "\tl2x_Sl_tref       : %f\n", l2x->l2x_Sl_tref);
    fprintf(LOG_DEST, "\tl2x_Sl_qref       : %f\n", l2x->l2x_Sl_qref);
    fprintf(LOG_DEST, "\tl2x_Sl_avsdr      : %f\n", l2x->l2x_Sl_avsdr);
    fprintf(LOG_DEST, "\tl2x_Sl_anidr      : %f\n", l2x->l2x_Sl_anidr);
    fprintf(LOG_DEST, "\tl2x_Sl_avsdf      : %f\n", l2x->l2x_Sl_avsdf);
    fprintf(LOG_DEST, "\tl2x_Sl_anidf      : %f\n", l2x->l2x_Sl_anidf);
    fprintf(LOG_DEST, "\tl2x_Sl_snowh      : %f\n", l2x->l2x_Sl_snowh);
    fprintf(LOG_DEST, "\tl2x_Sl_u10        : %f\n", l2x->l2x_Sl_u10);
    fprintf(LOG_DEST, "\tl2x_Sl_ddvel      : %f\n", l2x->l2x_Sl_ddvel);
    fprintf(LOG_DEST, "\tl2x_Sl_fv         : %f\n", l2x->l2x_Sl_fv);
    fprintf(LOG_DEST, "\tl2x_Sl_ram1       : %f\n", l2x->l2x_Sl_ram1);
    fprintf(LOG_DEST, "\tl2x_Sl_logz0      : %f\n", l2x->l2x_Sl_logz0);
    fprintf(LOG_DEST, "\tl2x_Fall_taux     : %f\n", l2x->l2x_Fall_taux);
    fprintf(LOG_DEST, "\tl2x_Fall_tauy     : %f\n", l2x->l2x_Fall_tauy);
    fprintf(LOG_DEST, "\tl2x_Fall_lat      : %f\n", l2x->l2x_Fall_lat);
    fprintf(LOG_DEST, "\tl2x_Fall_sen      : %f\n", l2x->l2x_Fall_sen);
    fprintf(LOG_DEST, "\tl2x_Fall_lwup     : %f\n", l2x->l2x_Fall_lwup);
    fprintf(LOG_DEST, "\tl2x_Fall_evap     : %f\n", l2x->l2x_Fall_evap);
    fprintf(LOG_DEST, "\tl2x_Fall_swnet    : %f\n", l2x->l2x_Fall_swnet);
    fprintf(LOG_DEST, "\tl2x_Fall_fco2_lnd : %f\n", l2x->l2x_Fall_fco2_lnd);
    fprintf(LOG_DEST, "\tl2x_Fall_flxdst1  : %f\n", l2x->l2x_Fall_flxdst1);
    fprintf(LOG_DEST, "\tl2x_Fall_flxdst2  : %f\n", l2x->l2x_Fall_flxdst2);
    fprintf(LOG_DEST, "\tl2x_Fall_flxdst3  : %f\n", l2x->l2x_Fall_flxdst3);
    fprintf(LOG_DEST, "\tl2x_Fall_flxdst4  : %f\n", l2x->l2x_Fall_flxdst4);
    fprintf(LOG_DEST, "\tl2x_Fall_flxvoc   : %f\n", l2x->l2x_Fall_flxvoc);
    fprintf(LOG_DEST, "\tl2x_Flrl_rofliq   : %f\n", l2x->l2x_Flrl_rofliq);
    fprintf(LOG_DEST, "\tl2x_Flrl_rofice   : %f\n", l2x->l2x_Flrl_rofice);
}
