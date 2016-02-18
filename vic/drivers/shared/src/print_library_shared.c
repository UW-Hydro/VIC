/******************************************************************************
 * @section DESCRIPTION
 *
 * Print library.
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
#include <vic_driver_shared.h>

/******************************************************************************
 * @brief    Print dell data structure.
 *****************************************************************************/
void
print_cell_data(cell_data_struct *cell,
                size_t            nlayers,
                size_t            nfrost)
{
    size_t i;

    fprintf(LOG_DEST, "cell_data:\n");
    fprintf(LOG_DEST, "\taero_resist :");
    for (i = 0; i < 2; i++) {
        fprintf(LOG_DEST, "\t%.4f", cell->aero_resist[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tasat        : %.4f\n", cell->asat);
    fprintf(LOG_DEST, "\tbaseflow    : %.4f\n", cell->baseflow);
    fprintf(LOG_DEST, "\tCLitter     : %.4f\n", cell->CLitter);
    fprintf(LOG_DEST, "\tCInter      : %.4f\n", cell->CInter);
    fprintf(LOG_DEST, "\tCSlow       : %.4f\n", cell->CSlow);
    fprintf(LOG_DEST, "\tinflow      : %.4f\n", cell->inflow);
    fprintf(LOG_DEST, "\tpot_evap    : %.4f\n", cell->pot_evap);
    fprintf(LOG_DEST, "\trunoff      : %.4f\n", cell->runoff);
    for (i = 0; i < nlayers; i++) {
        fprintf(LOG_DEST, "\tlayer %zd   :\n", i);
        print_layer_data(&(cell->layer[i]), nfrost);
    }
    fprintf(LOG_DEST, "\tRhLitter    : %.4f\n", cell->RhLitter);
    fprintf(LOG_DEST, "\tRhLitter2Atm: %.4f\n", cell->RhLitter2Atm);
    fprintf(LOG_DEST, "\tRhInter     : %.4f\n", cell->RhInter);
    fprintf(LOG_DEST, "\tRhSlow      : %.4f\n", cell->RhSlow);
    fprintf(LOG_DEST, "\tRhTot       : %.4f\n", cell->RhTot);
    fprintf(LOG_DEST, "\trootmoist   : %.4f\n", cell->rootmoist);
    fprintf(LOG_DEST, "\twetness     : %.4f\n", cell->wetness);
    fprintf(LOG_DEST, "\tzwt         : %.4f\n", cell->zwt);
    fprintf(LOG_DEST, "\tzwt_lumped  : %.4f\n", cell->zwt_lumped);
}

/******************************************************************************
 * @brief    Print day-month-year structure.
 *****************************************************************************/
void
print_dmy(dmy_struct *dmy)
{
    fprintf(LOG_DEST, "dmy:\n");
    fprintf(LOG_DEST, "\tday        : %hu\n", dmy->day);
    fprintf(LOG_DEST, "\tday_in_year: %hu\n", dmy->day_in_year);
    fprintf(LOG_DEST, "\tseconds    : %u\n", dmy->dayseconds);
    fprintf(LOG_DEST, "\tmonth      : %hu\n", dmy->month);
    fprintf(LOG_DEST, "\tyear       : %u\n", dmy->year);
}

/******************************************************************************
 * @brief    Print energy balance structure.
 *****************************************************************************/
void
print_energy_bal(energy_bal_struct *eb,
                 size_t             nnodes,
                 size_t             nfronts)
{
    size_t i;

    fprintf(LOG_DEST, "energy_bal:\n");
    fprintf(LOG_DEST, "\tAlbedoLake       : %.4f\n", eb->AlbedoLake);
    fprintf(LOG_DEST, "\tAlbedoOver       : %.4f\n", eb->AlbedoOver);
    fprintf(LOG_DEST, "\tAlbedoUnder      : %.4f\n", eb->AlbedoUnder);
    fprintf(LOG_DEST, "\tCs               :");
    for (i = 0; i < 2; i++) {
        fprintf(LOG_DEST, "\t%.4f", eb->Cs[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tCs_node          :");
    for (i = 0; i < nnodes; i++) {
        fprintf(LOG_DEST, "\t%.4f", eb->Cs_node[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tfdepth           :");
    for (i = 0; i < nfronts; i++) {
        fprintf(LOG_DEST, "\t%.4f", eb->fdepth[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tfrozen           : %d\n", eb->frozen);
    fprintf(LOG_DEST, "\tice              :");
    for (i = 0; i < nnodes; i++) {
        fprintf(LOG_DEST, "\t%.4f", eb->ice[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tkappa            :");
    for (i = 0; i < 2; i++) {
        fprintf(LOG_DEST, "\t%.4f", eb->kappa[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tkappa_node       :");
    for (i = 0; i < nnodes; i++) {
        fprintf(LOG_DEST, "\t%.4f", eb->kappa_node[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tmoist            :");
    for (i = 0; i < nnodes; i++) {
        fprintf(LOG_DEST, "\t%.4f", eb->moist[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tNfrost           : %zu\n", eb->Nfrost);
    fprintf(LOG_DEST, "\tNthaw            : %zu\n", eb->Nthaw);
    fprintf(LOG_DEST, "\tT                :");
    for (i = 0; i < nnodes; i++) {
        fprintf(LOG_DEST, "\t%.4f", eb->T[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tT_fbflag         :");
    for (i = 0; i < nnodes; i++) {
        fprintf(LOG_DEST, "\t%d", eb->T_fbflag[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tT_fbcount        :");
    for (i = 0; i < nnodes; i++) {
        fprintf(LOG_DEST, "\t%d", eb->T_fbcount[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tT1_index         : %d\n", eb->T1_index);
    fprintf(LOG_DEST, "\tTcanopy          : %.4f\n", eb->Tcanopy);
    fprintf(LOG_DEST, "\tTcanopy_fbflag   : %d\n", eb->Tcanopy_fbflag);
    fprintf(LOG_DEST, "\tTcanopy_fbcount  : %d\n", eb->Tcanopy_fbcount);
    fprintf(LOG_DEST, "\ttdepth           :");
    for (i = 0; i < nfronts; i++) {
        fprintf(LOG_DEST, "\t%.4f", eb->tdepth[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tTfoliage         : %.4f\n", eb->Tfoliage);
    fprintf(LOG_DEST, "\tTfoliage_fbflag  : %d\n", eb->Tfoliage_fbflag);
    fprintf(LOG_DEST, "\tTfoliage_fbcount : %d\n", eb->Tfoliage_fbcount);
    fprintf(LOG_DEST, "\tTsurf            : %.4f\n", eb->Tsurf);
    fprintf(LOG_DEST, "\tTsurf_fbflag     : %d\n", eb->Tsurf_fbflag);
    fprintf(LOG_DEST, "\tTsurf_fbcount    : %d\n", eb->Tsurf_fbcount);
    fprintf(LOG_DEST, "\tunfrozen         : %.4f\n", eb->unfrozen);
    fprintf(LOG_DEST, "\tadvected_sensible: %.4f\n", eb->advected_sensible);
    fprintf(LOG_DEST, "\tadvection        : %.4f\n", eb->advection);
    fprintf(LOG_DEST, "\tAtmosError       : %.4f\n", eb->AtmosError);
    fprintf(LOG_DEST, "\tAtmosLatent      : %.4f\n", eb->AtmosLatent);
    fprintf(LOG_DEST, "\tAtmosLatentSub   : %.4f\n", eb->AtmosLatentSub);
    fprintf(LOG_DEST, "\tAtmosSensible    : %.4f\n", eb->AtmosSensible);
    fprintf(LOG_DEST, "\tcanopy_advection : %.4f\n", eb->canopy_advection);
    fprintf(LOG_DEST, "\tcanopy_latent    : %.4f\n", eb->canopy_latent);
    fprintf(LOG_DEST, "\tcanopy_latent_sub: %.4f\n", eb->canopy_latent_sub);
    fprintf(LOG_DEST, "\tcanopy_refreeze  : %.4f\n", eb->canopy_refreeze);
    fprintf(LOG_DEST, "\tcanopy_sensible  : %.4f\n", eb->canopy_sensible);
    fprintf(LOG_DEST, "\tdeltaCC          : %.4f\n", eb->deltaCC);
    fprintf(LOG_DEST, "\tdeltaH           : %.4f\n", eb->deltaH);
    fprintf(LOG_DEST, "\terror            : %.4f\n", eb->error);
    fprintf(LOG_DEST, "\tfusion           : %.4f\n", eb->fusion);
    fprintf(LOG_DEST, "\tgrnd_flux        : %.4f\n", eb->grnd_flux);
    fprintf(LOG_DEST, "\tlatent           : %.4f\n", eb->latent);
    fprintf(LOG_DEST, "\tlatent_sub       : %.4f\n", eb->latent_sub);
    fprintf(LOG_DEST, "\tlongwave         : %.4f\n", eb->longwave);
    fprintf(LOG_DEST, "\tLongOverIn       : %.4f\n", eb->LongOverIn);
    fprintf(LOG_DEST, "\tLongUnderIn      : %.4f\n", eb->LongUnderIn);
    fprintf(LOG_DEST, "\tLongUnderOut     : %.4f\n", eb->LongUnderOut);
    fprintf(LOG_DEST, "\tmelt_energy      : %.4f\n", eb->melt_energy);
    fprintf(LOG_DEST, "\tNetLongAtmos     : %.4f\n", eb->NetLongAtmos);
    fprintf(LOG_DEST, "\tNetLongOver      : %.4f\n", eb->NetLongOver);
    fprintf(LOG_DEST, "\tNetLongUnder     : %.4f\n", eb->NetLongUnder);
    fprintf(LOG_DEST, "\tNetShortAtmos    : %.4f\n", eb->NetShortAtmos);
    fprintf(LOG_DEST, "\tNetShortGrnd     : %.4f\n", eb->NetShortGrnd);
    fprintf(LOG_DEST, "\tNetShortOver     : %.4f\n", eb->NetShortOver);
    fprintf(LOG_DEST, "\tNetShortUnder    : %.4f\n", eb->NetShortUnder);
    fprintf(LOG_DEST, "\tout_long_canopy  : %.4f\n", eb->out_long_canopy);
    fprintf(LOG_DEST, "\tout_long_surface : %.4f\n", eb->out_long_surface);
    fprintf(LOG_DEST, "\trefreeze_energy  : %.4f\n", eb->refreeze_energy);
    fprintf(LOG_DEST, "\tsensible         : %.4f\n", eb->sensible);
    fprintf(LOG_DEST, "\tshortwave        : %.4f\n", eb->shortwave);
    fprintf(LOG_DEST, "\tShortOverIn      : %.4f\n", eb->ShortOverIn);
    fprintf(LOG_DEST, "\tShortUnderIn     : %.4f\n", eb->ShortUnderIn);
    fprintf(LOG_DEST, "\tsnow_flux        : %.4f\n", eb->snow_flux);
}

/******************************************************************************
 * @brief    Print filenames structure.
 *****************************************************************************/
void
print_filenames(filenames_struct *fnames)
{
    fprintf(LOG_DEST, "filenames:\n");
    fprintf(LOG_DEST, "\tforcing[0]   : %s\n", fnames->forcing[0]);
    fprintf(LOG_DEST, "\tforcing[1]   : %s\n", fnames->forcing[1]);
    fprintf(LOG_DEST, "\tf_path_pfx[0]: %s\n", fnames->f_path_pfx[0]);
    fprintf(LOG_DEST, "\tf_path_pfx[1]: %s\n", fnames->f_path_pfx[1]);
    fprintf(LOG_DEST, "\tglobal       : %s\n", fnames->global);
    fprintf(LOG_DEST, "\tconstants    : %s\n", fnames->constants);
    fprintf(LOG_DEST, "\tdomain       : %s\n", fnames->domain);
    fprintf(LOG_DEST, "\tinit_state   : %s\n", fnames->init_state);
    fprintf(LOG_DEST, "\tlakeparam    : %s\n", fnames->lakeparam);
    fprintf(LOG_DEST, "\tresult_dir   : %s\n", fnames->result_dir);
    fprintf(LOG_DEST, "\tsnowband     : %s\n", fnames->snowband);
    fprintf(LOG_DEST, "\tsoil         : %s\n", fnames->soil);
    fprintf(LOG_DEST, "\tstatefile    : %s\n", fnames->statefile);
    fprintf(LOG_DEST, "\tveg          : %s\n", fnames->veg);
    fprintf(LOG_DEST, "\tveglib       : %s\n", fnames->veglib);
    fprintf(LOG_DEST, "\tlog_path     : %s\n", fnames->log_path);
}

/******************************************************************************
 * @brief    Print file path structure.
 *****************************************************************************/
void
print_filep(filep_struct *fp)
{
    fprintf(LOG_DEST, "filep:\n");
    fprintf(LOG_DEST, "\tforcing[0] : %p\n", fp->forcing[0]);
    fprintf(LOG_DEST, "\tforcing[1] : %p\n", fp->forcing[1]);
    fprintf(LOG_DEST, "\tglobalparam: %p\n", fp->globalparam);
    fprintf(LOG_DEST, "\tconstants  : %p\n", fp->constants);
    fprintf(LOG_DEST, "\tdomain     : %p\n", fp->domain);
    fprintf(LOG_DEST, "\tinit_state : %p\n", fp->init_state);
    fprintf(LOG_DEST, "\tlakeparam  : %p\n", fp->lakeparam);
    fprintf(LOG_DEST, "\tsnowband   : %p\n", fp->snowband);
    fprintf(LOG_DEST, "\tsoilparam  : %p\n", fp->soilparam);
    fprintf(LOG_DEST, "\tstatefile  : %p\n", fp->statefile);
    fprintf(LOG_DEST, "\tveglib     : %p\n", fp->veglib);
    fprintf(LOG_DEST, "\tvegparam   : %p\n", fp->vegparam);
    fprintf(LOG_DEST, "\tlogfile    : %p\n", fp->logfile);
}

/******************************************************************************
 * @brief    Print forcing type structure.
 *****************************************************************************/
void
print_force_type(force_type_struct *force_type)
{
    fprintf(LOG_DEST, "force_type:\n");
    fprintf(LOG_DEST, "\tSIGNED    : %d\n", force_type->SIGNED);
    fprintf(LOG_DEST, "\tSUPPLIED  : %d\n", force_type->SUPPLIED);
    fprintf(LOG_DEST, "\tmultiplier: %f\n", force_type->multiplier);
}

/******************************************************************************
 * @brief    Print global parameters structure.
 *****************************************************************************/
void
print_global_param(global_param_struct *gp)
{
    size_t i;

    fprintf(LOG_DEST, "global_param:\n");
    fprintf(LOG_DEST, "\twind_h              : %.4f\n", gp->wind_h);
    fprintf(LOG_DEST, "\tresolution          : %.4f\n", gp->resolution);
    fprintf(LOG_DEST, "\tdt                  : %.4f\n", gp->dt);
    fprintf(LOG_DEST, "\tsnow_dt             : %.4f\n", gp->snow_dt);
    fprintf(LOG_DEST, "\trunoff_dt           : %.4f\n", gp->runoff_dt);
    fprintf(LOG_DEST, "\tout_dt              : %.4f\n", gp->out_dt);
    fprintf(LOG_DEST, "\tmodel_steps_per_day : %zu\n", gp->model_steps_per_day);
    fprintf(LOG_DEST, "\tsnow_steps_per_day  : %zu\n", gp->snow_steps_per_day);
    fprintf(LOG_DEST, "\trunoff_steps_per_day: %zu\n",
            gp->runoff_steps_per_day);
    fprintf(LOG_DEST, "\tendday              : %hu\n", gp->endday);
    fprintf(LOG_DEST, "\tendmonth            : %hu\n", gp->endmonth);
    fprintf(LOG_DEST, "\tendyear             : %hu\n", gp->endyear);
    for (i = 0; i < 2; i++) {
        fprintf(LOG_DEST, "\tforceday[%zd]       : %hu\n", i, gp->forceday[i]);
        fprintf(LOG_DEST, "\tforcesec[%zd]       : %u\n", i, gp->forcesec[i]);
        fprintf(LOG_DEST, "\tforcemonth[%zd]     : %hu\n", i,
                gp->forcemonth[i]);
        fprintf(LOG_DEST, "\tforceoffset[%zd]    : %hu\n", i,
                gp->forceoffset[i]);
        fprintf(LOG_DEST, "\tforceskip[%zd]      : %u\n", i, gp->forceskip[i]);
        fprintf(LOG_DEST, "\tforceyear[%zd]      : %hu\n", i, gp->forceyear[i]);
    }
    fprintf(LOG_DEST, "\tnrecs               : %zu\n", gp->nrecs);
    fprintf(LOG_DEST, "\tskipyear            : %hu\n", gp->skipyear);
    fprintf(LOG_DEST, "\tstartday            : %hu\n", gp->startday);
    fprintf(LOG_DEST, "\tstartsec            : %u\n", gp->startsec);
    fprintf(LOG_DEST, "\tstartmonth          : %hu\n", gp->startmonth);
    fprintf(LOG_DEST, "\tstartyear           : %hu\n", gp->startyear);
    fprintf(LOG_DEST, "\tstateday            : %hu\n", gp->stateday);
    fprintf(LOG_DEST, "\tstatemonth          : %hu\n", gp->statemonth);
    fprintf(LOG_DEST, "\tstateyear           : %hu\n", gp->stateyear);
}

/******************************************************************************
 * @brief    Print lake_con_structure.
 *****************************************************************************/
void
print_lake_con(lake_con_struct *lcon,
               size_t           nlnodes)
{
    size_t i;

    fprintf(LOG_DEST, "lake_con:\n");
    fprintf(LOG_DEST, "\tnumnod   : %zu\n", lcon->numnod);
    fprintf(LOG_DEST, "\tz        :");
    for (i = 0; i < nlnodes; i++) {
        fprintf(LOG_DEST, "\t%.4f", lcon->z[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tbasin    :");
    for (i = 0; i < nlnodes; i++) {
        fprintf(LOG_DEST, "\t%.4f", lcon->basin[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tCl       :");
    for (i = 0; i < nlnodes; i++) {
        fprintf(LOG_DEST, "\t%.4f", lcon->Cl[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tb        : %.4f\n", lcon->b);
    fprintf(LOG_DEST, "\tmaxdepth : %.4f\n", lcon->maxdepth);
    fprintf(LOG_DEST, "\tmindepth : %.4f\n", lcon->mindepth);
    fprintf(LOG_DEST, "\tmaxvolume: %.4f\n", lcon->maxvolume);
    fprintf(LOG_DEST, "\tminvolume: %.4f\n", lcon->minvolume);
    fprintf(LOG_DEST, "\tbpercent : %.4f\n", lcon->bpercent);
    fprintf(LOG_DEST, "\trpercent : %.4f\n", lcon->rpercent);
    fprintf(LOG_DEST, "\twfrac    : %.4f\n", lcon->wfrac);
    fprintf(LOG_DEST, "\tdepth_in : %.4f\n", lcon->depth_in);
    fprintf(LOG_DEST, "\tlake_idx : %d\n", lcon->lake_idx);
}

/******************************************************************************
 * @brief    Print lake variables structure.
 *****************************************************************************/
void
print_lake_var(lake_var_struct *lvar,
               size_t           nlnodes,
               size_t           nfronts,
               size_t           nlayers,
               size_t           nnodes,
               size_t           nfrost)
{
    size_t i;

    fprintf(LOG_DEST, "lake_var:\n");
    fprintf(LOG_DEST, "\tactivenod      : %d\n", lvar->activenod);
    fprintf(LOG_DEST, "\tdz             : %.4f\n", lvar->dz);
    fprintf(LOG_DEST, "\tsurfdz         : %.4f\n", lvar->surfdz);
    fprintf(LOG_DEST, "\tldepth         : %.4f\n", lvar->ldepth);
    fprintf(LOG_DEST, "\tsurface        :");
    for (i = 0; i < nlnodes + 1; i++) {
        fprintf(LOG_DEST, "\t%.4f", lvar->surface[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tsarea          : %.4f\n", lvar->sarea);
    fprintf(LOG_DEST, "\tsarea_save     : %.4f\n", lvar->sarea_save);
    fprintf(LOG_DEST, "\tvolume         : %.4f\n", lvar->volume);
    fprintf(LOG_DEST, "\tvolume_save    : %.4f\n", lvar->volume_save);
    fprintf(LOG_DEST, "\ttemp           :");
    for (i = 0; i < nlnodes; i++) {
        fprintf(LOG_DEST, "\t%.4f", lvar->temp[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\ttempavg        : %.4f\n", lvar->tempavg);
    fprintf(LOG_DEST, "\tareai          : %.4f\n", lvar->areai);
    fprintf(LOG_DEST, "\tnew_ice_area   : %.4f\n", lvar->new_ice_area);
    fprintf(LOG_DEST, "\tice_water_eq   : %.4f\n", lvar->ice_water_eq);
    fprintf(LOG_DEST, "\thice           : %.4f\n", lvar->hice);
    fprintf(LOG_DEST, "\ttempi          : %.4f\n", lvar->tempi);
    fprintf(LOG_DEST, "\tswe            : %.4f\n", lvar->swe);
    fprintf(LOG_DEST, "\tswe_save       : %.4f\n", lvar->swe_save);
    fprintf(LOG_DEST, "\tsurf_temp      : %.4f\n", lvar->surf_temp);
    fprintf(LOG_DEST, "\tpack_temp      : %.4f\n", lvar->pack_temp);
    fprintf(LOG_DEST, "\tcoldcontent    : %.4f\n", lvar->coldcontent);
    fprintf(LOG_DEST, "\tsurf_water     : %.4f\n", lvar->surf_water);
    fprintf(LOG_DEST, "\tpack_water     : %.4f\n", lvar->pack_water);
    fprintf(LOG_DEST, "\tSAlbedo        : %.4f\n", lvar->SAlbedo);
    fprintf(LOG_DEST, "\tsdepth         : %.4f\n", lvar->sdepth);
    fprintf(LOG_DEST, "\taero_resist    : %.4f\n", lvar->aero_resist);
    fprintf(LOG_DEST, "\tdensity        :");
    for (i = 0; i < nlnodes; i++) {
        fprintf(LOG_DEST, "\t%.4f", lvar->density[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tbaseflow_in    : %.4f\n", lvar->baseflow_in);
    fprintf(LOG_DEST, "\tbaseflow_out   : %.4f\n", lvar->baseflow_out);
    fprintf(LOG_DEST, "\tchannel_in     : %.4f\n", lvar->channel_in);
    fprintf(LOG_DEST, "\tevapw          : %.4f\n", lvar->evapw);
    fprintf(LOG_DEST, "\tice_throughfall: %.4f\n", lvar->ice_throughfall);
    fprintf(LOG_DEST, "\tprec           : %.4f\n", lvar->prec);
    fprintf(LOG_DEST, "\trecharge       : %.4f\n", lvar->recharge);
    fprintf(LOG_DEST, "\trunoff_in      : %.4f\n", lvar->runoff_in);
    fprintf(LOG_DEST, "\trunoff_out     : %.4f\n", lvar->runoff_out);
    fprintf(LOG_DEST, "\tsnowmlt        : %.4f\n", lvar->snowmlt);
    fprintf(LOG_DEST, "\tvapor_flux     : %.4f\n", lvar->vapor_flux);
    print_snow_data(&(lvar->snow));
    print_energy_bal(&(lvar->energy), nnodes, nfronts);
    print_cell_data(&(lvar->soil), nlayers, nfrost);
}

/******************************************************************************
 * @brief    Print layer data structure.
 *****************************************************************************/
void
print_layer_data(layer_data_struct *ldata,
                 size_t             nfrost)
{
    size_t i;

    fprintf(LOG_DEST, "layer_data:\n");
    fprintf(LOG_DEST, "\tCs   : %.4f\n", ldata->Cs);
    fprintf(LOG_DEST, "\tT    : %.4f\n", ldata->T);
    fprintf(LOG_DEST, "\tevap : %.4f\n", ldata->evap);
    fprintf(LOG_DEST, "\tice  :");
    for (i = 0; i < nfrost; i++) {
        fprintf(LOG_DEST, "\t%.4f", ldata->ice[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tkappa: %.4f\n", ldata->kappa);
    fprintf(LOG_DEST, "\tmoist: %.4f\n", ldata->moist);
    fprintf(LOG_DEST, "\tphi  : %.4f\n", ldata->phi);
    fprintf(LOG_DEST, "\tzwt  : %.4f\n", ldata->zwt);
}

/******************************************************************************
 * @brief    Print options structure.
 *****************************************************************************/
void
print_option(option_struct *option)
{
    fprintf(LOG_DEST, "option:\n");
    fprintf(LOG_DEST, "\tAboveTreelineVeg     : %d\n",
            option->AboveTreelineVeg);
    fprintf(LOG_DEST, "\tAERO_RESIST_CANSNOW  : %d\n",
            option->AERO_RESIST_CANSNOW);
    fprintf(LOG_DEST, "\tBLOWING              : %d\n", option->BLOWING);
    fprintf(LOG_DEST, "\tBLOWING_VAR_THRESHOLD: %d\n",
            option->BLOWING_VAR_THRESHOLD);
    fprintf(LOG_DEST, "\tBLOWING_CALC_PROB    : %d\n",
            option->BLOWING_CALC_PROB);
    fprintf(LOG_DEST, "\tBLOWING_SIMPLE       : %d\n", option->BLOWING_SIMPLE);
    fprintf(LOG_DEST, "\tBLOWING_FETCH        : %d\n", option->BLOWING_FETCH);
    fprintf(LOG_DEST, "\tBLOWING_SPATIAL_WIND : %d\n",
            option->BLOWING_SPATIAL_WIND);
    fprintf(LOG_DEST, "\tCARBON               : %d\n", option->CARBON);
    fprintf(LOG_DEST, "\tCLOSE_ENERGY         : %d\n", option->CLOSE_ENERGY);
    fprintf(LOG_DEST, "\tCOMPUTE_TREELINE     : %d\n",
            option->COMPUTE_TREELINE);
    fprintf(LOG_DEST, "\tCONTINUEONERROR      : %d\n", option->CONTINUEONERROR);
    fprintf(LOG_DEST, "\tCORRPREC             : %d\n", option->CORRPREC);
    fprintf(LOG_DEST, "\tEQUAL_AREA           : %d\n", option->EQUAL_AREA);
    fprintf(LOG_DEST, "\tEXP_TRANS            : %d\n", option->EXP_TRANS);
    fprintf(LOG_DEST, "\tFROZEN_SOIL          : %d\n", option->FROZEN_SOIL);
    fprintf(LOG_DEST, "\tFULL_ENERGY          : %d\n", option->FULL_ENERGY);
    fprintf(LOG_DEST, "\tGRND_FLUX_TYPE       : %d\n", option->GRND_FLUX_TYPE);
    fprintf(LOG_DEST, "\tIMPLICIT             : %d\n", option->IMPLICIT);
    fprintf(LOG_DEST, "\tJULY_TAVG_SUPPLIED   : %d\n",
            option->JULY_TAVG_SUPPLIED);
    fprintf(LOG_DEST, "\tLAKES                : %d\n", option->LAKES);
    fprintf(LOG_DEST, "\tNcanopy              : %zu\n", option->Ncanopy);
    fprintf(LOG_DEST, "\tNfrost               : %zu\n", option->Nfrost);
    fprintf(LOG_DEST, "\tNlakenode            : %zu\n", option->Nlakenode);
    fprintf(LOG_DEST, "\tNlayer               : %zu\n", option->Nlayer);
    fprintf(LOG_DEST, "\tNnode                : %zu\n", option->Nnode);
    fprintf(LOG_DEST, "\tNOFLUX               : %d\n", option->NOFLUX);
    fprintf(LOG_DEST, "\tNVEGTYPES            : %zu\n", option->NVEGTYPES);
    fprintf(LOG_DEST, "\tRC_MODE              : %d\n", option->RC_MODE);
    fprintf(LOG_DEST, "\tROOT_ZONES           : %zu\n", option->ROOT_ZONES);
    fprintf(LOG_DEST, "\tQUICK_FLUX           : %d\n", option->QUICK_FLUX);
    fprintf(LOG_DEST, "\tQUICK_SOLVE          : %d\n", option->QUICK_SOLVE);
    fprintf(LOG_DEST, "\tSHARE_LAYER_MOIST    : %d\n",
            option->SHARE_LAYER_MOIST);
    fprintf(LOG_DEST, "\tSNOW_DENSITY         : %d\n", option->SNOW_DENSITY);
    fprintf(LOG_DEST, "\tSNOW_BAND            : %zu\n", option->SNOW_BAND);
    fprintf(LOG_DEST, "\tSPATIAL_FROST        : %d\n", option->SPATIAL_FROST);
    fprintf(LOG_DEST, "\tSPATIAL_SNOW         : %d\n", option->SPATIAL_SNOW);
    fprintf(LOG_DEST, "\tTFALLBACK            : %d\n", option->TFALLBACK);
    fprintf(LOG_DEST, "\tBASEFLOW             : %d\n", option->BASEFLOW);
    fprintf(LOG_DEST, "\tGRID_DECIMAL         : %d\n", option->GRID_DECIMAL);
    fprintf(LOG_DEST, "\tVEGLIB_PHOTO         : %d\n", option->VEGLIB_PHOTO);
    fprintf(LOG_DEST, "\tVEGLIB_FCAN      : %d\n", option->VEGLIB_FCAN);
    fprintf(LOG_DEST, "\tVEGPARAM_ALB         : %d\n", option->VEGPARAM_ALB);
    fprintf(LOG_DEST, "\tVEGPARAM_LAI         : %d\n", option->VEGPARAM_LAI);
    fprintf(LOG_DEST, "\tVEGPARAM_FCAN    : %d\n",
            option->VEGPARAM_FCAN);
    fprintf(LOG_DEST, "\tALB_SRC              : %d\n", option->ALB_SRC);
    fprintf(LOG_DEST, "\tLAI_SRC              : %d\n", option->LAI_SRC);
    fprintf(LOG_DEST, "\tFCAN_SRC         : %d\n", option->FCAN_SRC);
    fprintf(LOG_DEST, "\tLAKE_PROFILE         : %d\n", option->LAKE_PROFILE);
    fprintf(LOG_DEST, "\tORGANIC_FRACT        : %d\n", option->ORGANIC_FRACT);
    fprintf(LOG_DEST, "\tBINARY_STATE_FILE    : %d\n",
            option->BINARY_STATE_FILE);
    fprintf(LOG_DEST, "\tINIT_STATE           : %d\n", option->INIT_STATE);
    fprintf(LOG_DEST, "\tSAVE_STATE           : %d\n", option->SAVE_STATE);
    fprintf(LOG_DEST, "\tALMA_OUTPUT          : %d\n", option->ALMA_OUTPUT);
    fprintf(LOG_DEST, "\tBINARY_OUTPUT        : %d\n", option->BINARY_OUTPUT);
    fprintf(LOG_DEST, "\tCOMPRESS             : %d\n", option->COMPRESS);
    fprintf(LOG_DEST, "\tMOISTFRACT           : %d\n", option->MOISTFRACT);
    fprintf(LOG_DEST, "\tNoutfiles            : %zu\n", option->Noutfiles);
    fprintf(LOG_DEST, "\tPRT_HEADER           : %d\n", option->PRT_HEADER);
    fprintf(LOG_DEST, "\tPRT_SNOW_BAND        : %d\n", option->PRT_SNOW_BAND);
}

/******************************************************************************
 * @brief    Print out data structure.
 *****************************************************************************/
void
print_out_data(out_data_struct *out,
               size_t           nelem)
{
    size_t i;

    fprintf(LOG_DEST, "out_data:\n");
    fprintf(LOG_DEST, "\tvarname: %s\n", out->varname);
    fprintf(LOG_DEST, "\twrite: %d\n", out->write);
    fprintf(LOG_DEST, "\tformat: %s\n", out->format);
    fprintf(LOG_DEST, "\ttype: %d\n", out->type);
    fprintf(LOG_DEST, "\tmult: %.4f\n", out->mult);
    fprintf(LOG_DEST, "\tnelem: %d\n", out->nelem);
    fprintf(LOG_DEST, "\tdata:");
    for (i = 0; i < nelem; i++) {
        fprintf(LOG_DEST, "\t%.4f", out->data[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\taggdata:");
    for (i = 0; i < nelem; i++) {
        fprintf(LOG_DEST, "\t%.4f", out->aggdata[i]);
    }
    fprintf(LOG_DEST, "\n");
}

/******************************************************************************
 * @brief    Print out data file structure.
 *****************************************************************************/
void
print_out_data_file(out_data_file_struct *outf)
{
    fprintf(LOG_DEST, "\tprefix: %s\n", outf->prefix);
    fprintf(LOG_DEST, "\tfilename: %s\n", outf->filename);
    fprintf(LOG_DEST, "\tfh: %p\n", outf->fh);
    fprintf(LOG_DEST, "\tnvars: %zu\n", outf->nvars);
    fprintf(LOG_DEST, "\tvarid: %p\n", outf->varid);
}

/******************************************************************************
 * @brief    print param set structure.
 *****************************************************************************/
void
print_param_set(param_set_struct *param_set)
{
    size_t i;

    fprintf(LOG_DEST, "param_set:\n");
    for (i = 0; i < N_FORCING_TYPES; i++) {
        print_force_type(&(param_set->TYPE[i]));
    }
    fprintf(LOG_DEST, "\tFORCE_DT    : %.4f %.4f\n", param_set->FORCE_DT[0],
            param_set->FORCE_DT[1]);
    fprintf(LOG_DEST, "\tFORCE_ENDIAN: %d %d\n", param_set->FORCE_ENDIAN[0],
            param_set->FORCE_ENDIAN[1]);
    fprintf(LOG_DEST, "\tFORCE_FORMAT: %d %d\n", param_set->FORCE_FORMAT[0],
            param_set->FORCE_FORMAT[1]);
    fprintf(LOG_DEST, "\tFORCE_INDEX :\n");
    for (i = 0; i < N_FORCING_TYPES; i++) {
        fprintf(LOG_DEST, "\t\t%zd: %d %d\n", i, param_set->FORCE_INDEX[0][i],
                param_set->FORCE_INDEX[1][i]);
    }
    fprintf(LOG_DEST, "\tN_TYPES     : %zu %zu\n", param_set->N_TYPES[0],
            param_set->N_TYPES[1]);
}

/******************************************************************************
 * @brief    Print model parameters.
 *****************************************************************************/
void
print_parameters(parameters_struct *param)
{
    fprintf(LOG_DEST, "parameters:\n");
    fprintf(LOG_DEST, "\tLAPSE_RATE: %.4f\n", param->LAPSE_RATE);
    fprintf(LOG_DEST, "\tGAUGE_HEIGHT: %.4f\n", param->GAUGE_HEIGHT);
    fprintf(LOG_DEST, "\tHUGE_RESIST: %.4f\n", param->HUGE_RESIST);
    fprintf(LOG_DEST, "\tALBEDO_BARE_SOIL: %.4f\n", param->ALBEDO_BARE_SOIL);
    fprintf(LOG_DEST, "\tALBEDO_H20_SURF: %.4f\n", param->ALBEDO_H20_SURF);
    fprintf(LOG_DEST, "\tEMISS_GRND: %.4f\n", param->EMISS_GRND);
    fprintf(LOG_DEST, "\tEMISS_VEG: %.4f\n", param->EMISS_VEG);
    fprintf(LOG_DEST, "\tEMISS_ICE: %.4f\n", param->EMISS_ICE);
    fprintf(LOG_DEST, "\tEMISS_SNOW: %.4f\n", param->EMISS_SNOW);
    fprintf(LOG_DEST, "\tEMISS_H2O: %.4f\n", param->EMISS_H2O);
    fprintf(LOG_DEST, "\tSOIL_RESID_MOIST: %.4f\n", param->SOIL_RESID_MOIST);
    fprintf(LOG_DEST, "\tSOIL_SLAB_MOIST_FRACT: %.4f\n",
            param->SOIL_SLAB_MOIST_FRACT);
    fprintf(LOG_DEST, "\tVEG_LAI_SNOW_MULTIPLIER: %.4f\n",
            param->VEG_LAI_SNOW_MULTIPLIER);
    fprintf(LOG_DEST, "\tVEG_MIN_INTERCEPTION_STORAGE: %.4f\n",
            param->VEG_MIN_INTERCEPTION_STORAGE);
    fprintf(LOG_DEST, "\tVEG_LAI_WATER_FACTOR: %.4f\n",
            param->VEG_LAI_WATER_FACTOR);
    fprintf(LOG_DEST, "\tCANOPY_CLOSURE: %.4f\n", param->CANOPY_CLOSURE);
    fprintf(LOG_DEST, "\tCANOPY_RSMAX: %.4f\n", param->CANOPY_RSMAX);
    fprintf(LOG_DEST, "\tCANOPY_VPDMINFACTOR: %.4f\n",
            param->CANOPY_VPDMINFACTOR);
    fprintf(LOG_DEST, "\tLAKE_TMELT: %.4f\n", param->LAKE_TMELT);
    fprintf(LOG_DEST, "\tLAKE_MAX_SURFACE: %.4f\n", param->LAKE_MAX_SURFACE);
    fprintf(LOG_DEST, "\tLAKE_BETA: %.4f\n", param->LAKE_BETA);
    fprintf(LOG_DEST, "\tLAKE_FRACMIN: %.4f\n", param->LAKE_FRACMIN);
    fprintf(LOG_DEST, "\tLAKE_FRACLIM: %.4f\n", param->LAKE_FRACLIM);
    fprintf(LOG_DEST, "\tLAKE_DM: %.4f\n", param->LAKE_DM);
    fprintf(LOG_DEST, "\tLAKE_SNOWCRIT: %.4f\n", param->LAKE_SNOWCRIT);
    fprintf(LOG_DEST, "\tLAKE_ZWATER: %.4f\n", param->LAKE_ZWATER);
    fprintf(LOG_DEST, "\tLAKE_ZSNOW: %.4f\n", param->LAKE_ZSNOW);
    fprintf(LOG_DEST, "\tLAKE_RHOSNOW: %.4f\n", param->LAKE_RHOSNOW);
    fprintf(LOG_DEST, "\tLAKE_CONDI: %.4f\n", param->LAKE_CONDI);
    fprintf(LOG_DEST, "\tLAKE_CONDS: %.4f\n", param->LAKE_CONDS);
    fprintf(LOG_DEST, "\tLAKE_LAMISW: %.4f\n", param->LAKE_LAMISW);
    fprintf(LOG_DEST, "\tLAKE_LAMILW: %.4f\n", param->LAKE_LAMILW);
    fprintf(LOG_DEST, "\tLAKE_LAMSSW: %.4f\n", param->LAKE_LAMSSW);
    fprintf(LOG_DEST, "\tLAKE_LAMSLW: %.4f\n", param->LAKE_LAMSLW);
    fprintf(LOG_DEST, "\tLAKE_LAMWSW: %.4f\n", param->LAKE_LAMWSW);
    fprintf(LOG_DEST, "\tLAKE_LAMWLW: %.4f\n", param->LAKE_LAMWLW);
    fprintf(LOG_DEST, "\tLAKE_A1: %.4f\n", param->LAKE_A1);
    fprintf(LOG_DEST, "\tLAKE_A2: %.4f\n", param->LAKE_A2);
    fprintf(LOG_DEST, "\tLAKE_QWTAU: %.4f\n", param->LAKE_QWTAU);
    fprintf(LOG_DEST, "\tLAKE_MAX_ITER: %d\n", param->LAKE_MAX_ITER);
    fprintf(LOG_DEST, "\tSVP_A: %.4f\n", param->SVP_A);
    fprintf(LOG_DEST, "\tSVP_B: %.4f\n", param->SVP_B);
    fprintf(LOG_DEST, "\tSVP_C: %.4f\n", param->SVP_C);
    fprintf(LOG_DEST, "\tPHOTO_OMEGA: %.4f\n", param->PHOTO_OMEGA);
    fprintf(LOG_DEST, "\tPHOTO_LAIMAX: %.4f\n", param->PHOTO_LAIMAX);
    fprintf(LOG_DEST, "\tPHOTO_LAILIMIT: %.4f\n", param->PHOTO_LAILIMIT);
    fprintf(LOG_DEST, "\tPHOTO_LAIMIN: %.4f\n", param->PHOTO_LAIMIN);
    fprintf(LOG_DEST, "\tPHOTO_EPAR: %.4f\n", param->PHOTO_EPAR);
    fprintf(LOG_DEST, "\tPHOTO_FCMAX: %.4f\n", param->PHOTO_FCMAX);
    fprintf(LOG_DEST, "\tPHOTO_FCMIN: %.4f\n", param->PHOTO_FCMIN);
    fprintf(LOG_DEST, "\tPHOTO_ZENITHMIN: %.4f\n", param->PHOTO_ZENITHMIN);
    fprintf(LOG_DEST, "\tPHOTO_ZENITHMINPAR: %.4f\n",
            param->PHOTO_ZENITHMINPAR);
    fprintf(LOG_DEST, "\tPHOTO_ALBSOIPARMIN: %.4f\n",
            param->PHOTO_ALBSOIPARMIN);
    fprintf(LOG_DEST, "\tPHOTO_MINMAXETRANS: %.4f\n",
            param->PHOTO_MINMAXETRANS);
    fprintf(LOG_DEST, "\tPHOTO_MINSTOMCOND: %.4f\n", param->PHOTO_MINSTOMCOND);
    fprintf(LOG_DEST, "\tPHOTO_FCI1C3: %.4f\n", param->PHOTO_FCI1C3);
    fprintf(LOG_DEST, "\tPHOTO_FCI1C4: %.4f\n", param->PHOTO_FCI1C4);
    fprintf(LOG_DEST, "\tPHOTO_OX: %.4f\n", param->PHOTO_OX);
    fprintf(LOG_DEST, "\tPHOTO_KC: %.4f\n", param->PHOTO_KC);
    fprintf(LOG_DEST, "\tPHOTO_KO: %.4f\n", param->PHOTO_KO);
    fprintf(LOG_DEST, "\tPHOTO_EC: %.4f\n", param->PHOTO_EC);
    fprintf(LOG_DEST, "\tPHOTO_EO: %.4f\n", param->PHOTO_EO);
    fprintf(LOG_DEST, "\tPHOTO_EV: %.4f\n", param->PHOTO_EV);
    fprintf(LOG_DEST, "\tPHOTO_ER: %.4f\n", param->PHOTO_ER);
    fprintf(LOG_DEST, "\tPHOTO_ALC3: %.4f\n", param->PHOTO_ALC3);
    fprintf(LOG_DEST, "\tPHOTO_FRDC3: %.4f\n", param->PHOTO_FRDC3);
    fprintf(LOG_DEST, "\tPHOTO_EK: %.4f\n", param->PHOTO_EK);
    fprintf(LOG_DEST, "\tPHOTO_ALC4: %.4f\n", param->PHOTO_ALC4);
    fprintf(LOG_DEST, "\tPHOTO_FRDC4: %.4f\n", param->PHOTO_FRDC4);
    fprintf(LOG_DEST, "\tPHOTO_THETA: %.4f\n", param->PHOTO_THETA);
    fprintf(LOG_DEST, "\tPHOTO_FRLEAF: %.4f\n", param->PHOTO_FRLEAF);
    fprintf(LOG_DEST, "\tPHOTO_FRGROWTH: %.4f\n", param->PHOTO_FRGROWTH);
    fprintf(LOG_DEST, "\tSRESP_E0_LT: %.4f\n", param->SRESP_E0_LT);
    fprintf(LOG_DEST, "\tSRESP_T0_LT: %.4f\n", param->SRESP_T0_LT);
    fprintf(LOG_DEST, "\tSRESP_WMINFM: %.4f\n", param->SRESP_WMINFM);
    fprintf(LOG_DEST, "\tSRESP_WMAXFM: %.4f\n", param->SRESP_WMAXFM);
    fprintf(LOG_DEST, "\tSRESP_WOPTFM: %.4f\n", param->SRESP_WOPTFM);
    fprintf(LOG_DEST, "\tSRESP_RHSAT: %.4f\n", param->SRESP_RHSAT);
    fprintf(LOG_DEST, "\tSRESP_RFACTOR: %.4f\n", param->SRESP_RFACTOR);
    fprintf(LOG_DEST, "\tSRESP_TAULITTER: %.4f\n", param->SRESP_TAULITTER);
    fprintf(LOG_DEST, "\tSRESP_TAUINTER: %.4f\n", param->SRESP_TAUINTER);
    fprintf(LOG_DEST, "\tSRESP_TAUSLOW: %.4f\n", param->SRESP_TAUSLOW);
    fprintf(LOG_DEST, "\tSRESP_FAIR: %.4f\n", param->SRESP_FAIR);
    fprintf(LOG_DEST, "\tSRESP_FINTER: %.4f\n", param->SRESP_FINTER);
    fprintf(LOG_DEST, "\tSNOW_MAX_SURFACE_SWE: %.4f\n",
            param->SNOW_MAX_SURFACE_SWE);
    fprintf(LOG_DEST, "\tSNOW_LIQUID_WATER_CAPACITY: %.4f\n",
            param->SNOW_LIQUID_WATER_CAPACITY);
    fprintf(LOG_DEST, "\tSNOW_NEW_SNOW_DENSITY: %.4f\n",
            param->SNOW_NEW_SNOW_DENSITY);
    fprintf(LOG_DEST, "\tSNOW_DENS_DMLIMIT: %.4f\n", param->SNOW_DENS_DMLIMIT);
    fprintf(LOG_DEST, "\tSNOW_DENS_MAX_CHANGE: %.4f\n",
            param->SNOW_DENS_MAX_CHANGE);
    fprintf(LOG_DEST, "\tSNOW_DENS_ETA0: %.4f\n", param->SNOW_DENS_ETA0);
    fprintf(LOG_DEST, "\tSNOW_DENS_C1: %.4f\n", param->SNOW_DENS_C1);
    fprintf(LOG_DEST, "\tSNOW_DENS_C2: %.4f\n", param->SNOW_DENS_C2);
    fprintf(LOG_DEST, "\tSNOW_DENS_C5: %.4f\n", param->SNOW_DENS_C5);
    fprintf(LOG_DEST, "\tSNOW_DENS_C6: %.4f\n", param->SNOW_DENS_C6);
    fprintf(LOG_DEST, "\tSNOW_DENS_F: %.4f\n", param->SNOW_DENS_F);
    fprintf(LOG_DEST, "\tSNOW_MIN_SWQ_EB_THRES: %.4f\n",
            param->SNOW_MIN_SWQ_EB_THRES);
    fprintf(LOG_DEST, "\tSNOW_A1: %.4f\n", param->SNOW_A1);
    fprintf(LOG_DEST, "\tSNOW_A2: %.4f\n", param->SNOW_A2);
    fprintf(LOG_DEST, "\tSNOW_L1: %.4f\n", param->SNOW_L1);
    fprintf(LOG_DEST, "\tSNOW_L2: %.4f\n", param->SNOW_L2);
    fprintf(LOG_DEST, "\tSNOW_NEW_SNOW_ALB: %.4f\n", param->SNOW_NEW_SNOW_ALB);
    fprintf(LOG_DEST, "\tSNOW_ALB_ACCUM_A: %.4f\n", param->SNOW_ALB_ACCUM_A);
    fprintf(LOG_DEST, "\tSNOW_ALB_ACCUM_B: %.4f\n", param->SNOW_ALB_ACCUM_B);
    fprintf(LOG_DEST, "\tSNOW_ALB_THAW_A: %.4f\n", param->SNOW_ALB_THAW_A);
    fprintf(LOG_DEST, "\tSNOW_ALB_THAW_B: %.4f\n", param->SNOW_ALB_THAW_B);
    fprintf(LOG_DEST, "\tSNOW_TRACESNOW: %.4f\n", param->SNOW_TRACESNOW);
    fprintf(LOG_DEST, "\tSNOW_CONDUCT: %.4f\n", param->SNOW_CONDUCT);
    fprintf(LOG_DEST, "\tSNOW_MAX_SNOW_TEMP: %.4f\n",
            param->SNOW_MAX_SNOW_TEMP);
    fprintf(LOG_DEST, "\tSNOW_MIN_RAIN_TEMP: %.4f\n",
            param->SNOW_MIN_RAIN_TEMP);
    fprintf(LOG_DEST, "\tBLOWING_KA: %.4f\n", param->BLOWING_KA);
    fprintf(LOG_DEST, "\tBLOWING_CSALT: %.4f\n", param->BLOWING_CSALT);
    fprintf(LOG_DEST, "\tBLOWING_UTHRESH: %.4f\n", param->BLOWING_UTHRESH);
    fprintf(LOG_DEST, "\tBLOWING_KIN_VIS: %.4f\n", param->BLOWING_KIN_VIS);
    fprintf(LOG_DEST, "\tBLOWING_MAX_ITER: %d\n", param->BLOWING_MAX_ITER);
    fprintf(LOG_DEST, "\tBLOWING_K: %d\n", param->BLOWING_K);
    fprintf(LOG_DEST, "\tBLOWING_SETTLING: %.4f\n", param->BLOWING_SETTLING);
    fprintf(LOG_DEST, "\tBLOWING_NUMINCS: %d\n", param->BLOWING_NUMINCS);
    fprintf(LOG_DEST, "\tTREELINE_TEMPERATURE: %.4f\n",
            param->TREELINE_TEMPERATURE);
    fprintf(LOG_DEST, "\tSNOW_DT: %.4f\n", param->SNOW_DT);
    fprintf(LOG_DEST, "\tSURF_DT: %.4f\n", param->SURF_DT);
    fprintf(LOG_DEST, "\tSOIL_DT: %.4f\n", param->SOIL_DT);
    fprintf(LOG_DEST, "\tCANOPY_DT: %.4f\n", param->CANOPY_DT);
    fprintf(LOG_DEST, "\tCANOPY_VP: %.4f\n", param->CANOPY_VP);
    fprintf(LOG_DEST, "\tTOL_GRND: %.4f\n", param->TOL_GRND);
    fprintf(LOG_DEST, "\tTOL_OVER: %.4f\n", param->TOL_OVER);
    fprintf(LOG_DEST, "\tFROZEN_MAXITER: %d\n", param->FROZEN_MAXITER);
    fprintf(LOG_DEST, "\tNEWT_RAPH_MAXTRIAL: %d\n", param->NEWT_RAPH_MAXTRIAL);
    fprintf(LOG_DEST, "\tNEWT_RAPH_TOLX: %.4f\n", param->NEWT_RAPH_TOLX);
    fprintf(LOG_DEST, "\tNEWT_RAPH_TOLF: %.4f\n", param->NEWT_RAPH_TOLF);
    fprintf(LOG_DEST, "\tNEWT_RAPH_R_MAX: %.4f\n", param->NEWT_RAPH_R_MAX);
    fprintf(LOG_DEST, "\tNEWT_RAPH_R_MIN: %.4f\n", param->NEWT_RAPH_R_MIN);
    fprintf(LOG_DEST, "\tNEWT_RAPH_RELAX1: %.4f\n", param->NEWT_RAPH_RELAX1);
    fprintf(LOG_DEST, "\tNEWT_RAPH_RELAX2: %.4f\n", param->NEWT_RAPH_RELAX2);
    fprintf(LOG_DEST, "\tNEWT_RAPH_RELAX3: %.4f\n", param->NEWT_RAPH_RELAX3);
    fprintf(LOG_DEST, "\tNEWT_RAPH_EPS2: %.4f\n", param->NEWT_RAPH_EPS2);
    fprintf(LOG_DEST, "\tROOT_BRENT_MAXTRIES: %d\n",
            param->ROOT_BRENT_MAXTRIES);
    fprintf(LOG_DEST, "\tROOT_BRENT_MAXITER: %d\n", param->ROOT_BRENT_MAXITER);
    fprintf(LOG_DEST, "\tROOT_BRENT_TSTEP: %.4f\n", param->ROOT_BRENT_TSTEP);
    fprintf(LOG_DEST, "\tROOT_BRENT_T: %.4f\n", param->ROOT_BRENT_T);
    fprintf(LOG_DEST, "\tFROZEN_MAXITER: %d\n", param->FROZEN_MAXITER);
}

/******************************************************************************
 * @brief    Print save data structure.
 *****************************************************************************/
void
print_save_data(save_data_struct *save)
{
    fprintf(LOG_DEST, "save_data:\n");
    fprintf(LOG_DEST, "\ttotal_soil_moist: %.4f\n", save->total_soil_moist);
    fprintf(LOG_DEST, "\tsurfstor: %.4f\n", save->surfstor);
    fprintf(LOG_DEST, "\tswe: %.4f\n", save->swe);
    fprintf(LOG_DEST, "\twdew: %.4f\n", save->wdew);
}

/******************************************************************************
 * @brief     Print snow data structure.
 *****************************************************************************/
void
print_snow_data(snow_data_struct *snow)
{
    fprintf(LOG_DEST, "snow_data:\n");
    fprintf(LOG_DEST, "\talbedo            : %.4f\n", snow->albedo);
    fprintf(LOG_DEST, "\tcanopy_albedo     : %.4f\n", snow->canopy_albedo);
    fprintf(LOG_DEST, "\tcoldcontent       : %.4f\n", snow->coldcontent);
    fprintf(LOG_DEST, "\tcoverage          : %.4f\n", snow->coverage);
    fprintf(LOG_DEST, "\tdensity           : %.4f\n", snow->density);
    fprintf(LOG_DEST, "\tdepth             : %.4f\n", snow->depth);
    fprintf(LOG_DEST, "\tlast_snow         : %d\n", snow->last_snow);
    fprintf(LOG_DEST, "\tmax_snow_depth    : %.4f\n", snow->max_snow_depth);
    fprintf(LOG_DEST, "\tMELTING           : %d\n", snow->MELTING);
    fprintf(LOG_DEST, "\tpack_temp         : %.4f\n", snow->pack_temp);
    fprintf(LOG_DEST, "\tpack_water        : %.4f\n", snow->pack_water);
    fprintf(LOG_DEST, "\tsnow              : %d\n", snow->snow);
    fprintf(LOG_DEST, "\tsnow_canopy       : %.4f\n", snow->snow_canopy);
    fprintf(LOG_DEST, "\tstore_coverage    : %.4f\n", snow->store_coverage);
    fprintf(LOG_DEST, "\tstore_snow        : %d\n", snow->store_snow);
    fprintf(LOG_DEST, "\tstore_swq         : %.4f\n", snow->store_swq);
    fprintf(LOG_DEST, "\tsurf_temp         : %.4f\n", snow->surf_temp);
    fprintf(LOG_DEST, "\tsurf_temp_fbcount : %u\n", snow->surf_temp_fbcount);
    fprintf(LOG_DEST, "\tsurf_temp_fbflag  : %d\n", snow->surf_temp_fbflag);
    fprintf(LOG_DEST, "\tsurf_water        : %.4f\n", snow->surf_water);
    fprintf(LOG_DEST, "\tswq               : %.4f\n", snow->swq);
    fprintf(LOG_DEST, "\tsnow_distrib_slope: %.4f\n",
            snow->snow_distrib_slope);
    fprintf(LOG_DEST, "\ttmp_int_storage   : %.4f\n", snow->tmp_int_storage);
    fprintf(LOG_DEST, "\tblowing_flux      : %.4f\n", snow->blowing_flux);
    fprintf(LOG_DEST, "\tcanopy_vapor_flux : %.4f\n", snow->canopy_vapor_flux);
    fprintf(LOG_DEST, "\tmass_error        : %.4f\n", snow->mass_error);
    fprintf(LOG_DEST, "\tmelt              : %.4f\n", snow->melt);
    fprintf(LOG_DEST, "\tQnet              : %.4f\n", snow->Qnet);
    fprintf(LOG_DEST, "\tsurface_flux      : %.4f\n", snow->surface_flux);
    fprintf(LOG_DEST, "\ttransport         : %.4f\n", snow->transport);
    fprintf(LOG_DEST, "\tvapor_flux        : %.4f\n", snow->vapor_flux);
}

/******************************************************************************
 * @brief    Print soil_con_struct.
 *****************************************************************************/
void
print_soil_con(soil_con_struct *scon,
               size_t           nlayers,
               size_t           nnodes,
               size_t           nfrost,
               size_t           nbands,
               size_t           nzwt)
{
    size_t i;
    size_t j;

    fprintf(LOG_DEST, "soil_con:\n");
    fprintf(LOG_DEST, "\tFS_ACTIVE             : %d\n", scon->FS_ACTIVE);
    fprintf(LOG_DEST, "\tDs                    : %.4f\n", scon->Ds);
    fprintf(LOG_DEST, "\tDsmax                 : %.4f\n", scon->Dsmax);
    fprintf(LOG_DEST, "\tKsat                  :");
    for (i = 0; i < nlayers; i++) {
        fprintf(LOG_DEST, "\t%.4f", scon->Ksat[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tWcr                   :");
    for (i = 0; i < nlayers; i++) {
        fprintf(LOG_DEST, "\t%.4f", scon->Wcr[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tWpwp                  :");
    for (i = 0; i < nlayers; i++) {
        fprintf(LOG_DEST, "\t%.4f", scon->Wpwp[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tWs                    : %.4f\n", scon->Ws);
    fprintf(LOG_DEST, "\tAlbedoPar             : %.4f\n", scon->AlbedoPar);
    fprintf(LOG_DEST, "\talpha                 :");
    for (i = 0; i < nnodes; i++) {
        fprintf(LOG_DEST, "\t%.4f", scon->alpha[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tannual_prec           : %.4f\n", scon->annual_prec);
    fprintf(LOG_DEST, "\tavg_temp              : %.4f\n", scon->avg_temp);
    fprintf(LOG_DEST, "\tavgJulyAirTemp        : %.4f\n",
            scon->avgJulyAirTemp);
    fprintf(LOG_DEST, "\tb_infilt              : %.4f\n", scon->b_infilt);
    fprintf(LOG_DEST, "\tbeta                  :");
    for (i = 0; i < nnodes; i++) {
        fprintf(LOG_DEST, "\t%.4f", scon->beta[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tbubble                :");
    for (i = 0; i < nlayers; i++) {
        fprintf(LOG_DEST, "\t%.4f", scon->bubble[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tbubble_node           :");
    for (i = 0; i < nnodes; i++) {
        fprintf(LOG_DEST, "\t%.4f", scon->bubble_node[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tbulk_density          :");
    for (i = 0; i < nlayers; i++) {
        fprintf(LOG_DEST, "\t%.4f", scon->bulk_density[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tbulk_dens_min         :");
    for (i = 0; i < nlayers; i++) {
        fprintf(LOG_DEST, "\t%.4f", scon->bulk_dens_min[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tbulk_dens_org       :");
    for (i = 0; i < nlayers; i++) {
        fprintf(LOG_DEST, "\t%.4f", scon->bulk_dens_org[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tc                     : %.4f\n", scon->c);
    fprintf(LOG_DEST, "\tdepth                 :");
    for (i = 0; i < nlayers; i++) {
        fprintf(LOG_DEST, "\t%.4f", scon->depth[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tdp                    : %.4f\n", scon->dp);
    fprintf(LOG_DEST, "\tdz_node               :");
    for (i = 0; i < nnodes; i++) {
        fprintf(LOG_DEST, "\t%.4f", scon->dz_node[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tZsum_node             :");
    for (i = 0; i < nnodes; i++) {
        fprintf(LOG_DEST, "\t%.4f", scon->Zsum_node[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\texpt                  :");
    for (i = 0; i < nlayers; i++) {
        fprintf(LOG_DEST, "\t%.4f", scon->expt[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\texpt_node             :");
    for (i = 0; i < nnodes; i++) {
        fprintf(LOG_DEST, "\t%.4f", scon->expt_node[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tfrost_fract           :");
    for (i = 0; i < nfrost; i++) {
        fprintf(LOG_DEST, "\t%.4f", scon->frost_fract[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tfrost_slope           : %.4f\n", scon->frost_slope);
    fprintf(LOG_DEST, "\tgamma                 :");
    for (i = 0; i < nnodes; i++) {
        fprintf(LOG_DEST, "\t%.4f", scon->gamma[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tinit_moist            :");
    for (i = 0; i < nlayers; i++) {
        fprintf(LOG_DEST, "\t%.4f", scon->init_moist[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tmax_infil             : %.4f\n", scon->max_infil);
    fprintf(LOG_DEST, "\tmax_moist             :");
    for (i = 0; i < nlayers; i++) {
        fprintf(LOG_DEST, "\t%.4f", scon->max_moist[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tmax_moist_node        :");
    for (i = 0; i < nnodes; i++) {
        fprintf(LOG_DEST, "\t%.4f", scon->max_moist_node[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tmax_snow_distrib_slope: %.4f\n",
            scon->max_snow_distrib_slope);
    fprintf(LOG_DEST, "\tphi_s                 :");
    for (i = 0; i < nlayers; i++) {
        fprintf(LOG_DEST, "\t%.4f", scon->phi_s[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tporosity              :");
    for (i = 0; i < nlayers; i++) {
        fprintf(LOG_DEST, "\t%.4f", scon->porosity[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tquartz              :");
    for (i = 0; i < nlayers; i++) {
        fprintf(LOG_DEST, "\t%.4f", scon->quartz[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\torganic               :");
    for (i = 0; i < nlayers; i++) {
        fprintf(LOG_DEST, "\t%.4f", scon->organic[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tresid_moist           :");
    for (i = 0; i < nlayers; i++) {
        fprintf(LOG_DEST, "\t%.4f", scon->resid_moist[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\trough                 : %.4f\n", scon->rough);
    fprintf(LOG_DEST, "\tsnow_rough            : %.4f\n", scon->snow_rough);
    fprintf(LOG_DEST, "\tsoil_density          :");
    for (i = 0; i < nlayers; i++) {
        fprintf(LOG_DEST, "\t%.4f", scon->soil_density[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tsoil_dens_min         :");
    for (i = 0; i < nlayers; i++) {
        fprintf(LOG_DEST, "\t%.4f", scon->soil_dens_min[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tsoil_dens_org         :");
    for (i = 0; i < nlayers; i++) {
        fprintf(LOG_DEST, "\t%.4f", scon->soil_dens_org[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "BandElev                :");
    for (i = 0; i < nbands; i++) {
        fprintf(LOG_DEST, "\t%.4f", scon->BandElev[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "AreaFract               :");
    for (i = 0; i < nbands; i++) {
        fprintf(LOG_DEST, "\t%.4f", scon->AreaFract[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "Pfactor               :");
    for (i = 0; i < nbands; i++) {
        fprintf(LOG_DEST, "\t%.4f", scon->Pfactor[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "Tfactor               :");
    for (i = 0; i < nbands; i++) {
        fprintf(LOG_DEST, "\t%.4f", scon->Tfactor[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "AboveTreeLine         :");
    for (i = 0; i < nbands; i++) {
        fprintf(LOG_DEST, "\t%d", scon->AboveTreeLine[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\televation             : %.4f\n", scon->elevation);
    fprintf(LOG_DEST, "\tlat                   : %.4f\n", scon->lat);
    fprintf(LOG_DEST, "\tlng                   : %.4f\n", scon->lng);
    fprintf(LOG_DEST, "\tcell_area             : %.4f\n", scon->cell_area);
    fprintf(LOG_DEST, "\ttime_zone_lng         : %.4f\n", scon->time_zone_lng);
    fprintf(LOG_DEST, "\tgridcel               : %d\n", scon->gridcel);
    fprintf(LOG_DEST, "\tzwtvmoist_zwt         :");
    for (i = 0; i < nlayers + 2; i++) {
        for (j = 0; j < nzwt; j++) {
            fprintf(LOG_DEST, "\t%.4f", scon->zwtvmoist_zwt[i][j]);
        }
        fprintf(LOG_DEST, "\n\t\t\t");
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tzwtvmoist_moist       :");
    for (i = 0; i < nlayers + 2; i++) {
        for (j = 0; j < nzwt; j++) {
            fprintf(LOG_DEST, "\t%.4f", scon->zwtvmoist_moist[i][j]);
        }
        fprintf(LOG_DEST, "\n\t\t\t");
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tslope                 : %.4f\n", scon->slope);
    fprintf(LOG_DEST, "\taspect                : %.4f\n", scon->aspect);
    fprintf(LOG_DEST, "\tehoriz                : %.4f\n", scon->ehoriz);
    fprintf(LOG_DEST, "\twhoriz                : %.4f\n", scon->whoriz);
}

/******************************************************************************
 * @brief    Print veg_con structure.
 *****************************************************************************/
void
print_veg_con(veg_con_struct *vcon,
              size_t          nroots,
              char            blowing,
              char            lake,
              char            carbon,
              size_t          ncanopy)
{
    size_t i;

    fprintf(LOG_DEST, "veg_con:\n");
    fprintf(LOG_DEST, "\tCv              : %.4f\n", vcon->Cv);
    fprintf(LOG_DEST, "\tCv_sum          : %.4f\n", vcon->Cv_sum);
    fprintf(LOG_DEST, "\troot            :");
    for (i = 0; i < nroots; i++) {
        fprintf(LOG_DEST, "\t%.2f", vcon->root[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tzone_depth      :");
    for (i = 0; i < nroots; i++) {
        fprintf(LOG_DEST, "\t%.2f", vcon->zone_depth[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tzone_fract      :");
    for (i = 0; i < nroots; i++) {
        fprintf(LOG_DEST, "\t%.2f", vcon->zone_fract[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tveg_class       : %d\n", vcon->veg_class);
    fprintf(LOG_DEST, "\tvegetat_type_num: %zu\n", vcon->vegetat_type_num);
    if (blowing) {
        fprintf(LOG_DEST, "\tsigma_slope     : %.4f\n", vcon->sigma_slope);
        fprintf(LOG_DEST, "\tlag_one         : %.4f\n", vcon->lag_one);
        fprintf(LOG_DEST, "\tfetch           : %.4f\n", vcon->fetch);
    }
    if (lake) {
        fprintf(LOG_DEST, "\tLAKE            : %d\n", vcon->LAKE);
    }
    if (carbon) {
        fprintf(LOG_DEST, "\tCanopLayerBnd   :");
        for (i = 0; i < ncanopy; i++) {
            fprintf(LOG_DEST, "\t%.2f", vcon->CanopLayerBnd[i]);
        }
    }
}

/******************************************************************************
 * @brief    Print vegetation library variables.
 *****************************************************************************/
void
print_veg_lib(veg_lib_struct *vlib,
              char            carbon)
{
    size_t i;

    fprintf(LOG_DEST, "veg_lib:\n");
    fprintf(LOG_DEST, "\toverstory     : %d\n", vlib->overstory);
    fprintf(LOG_DEST, "\tLAI           :");
    for (i = 0; i < MONTHS_PER_YEAR; i++) {
        fprintf(LOG_DEST, "\t%.2f", vlib->LAI[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tWdmax         :");
    for (i = 0; i < MONTHS_PER_YEAR; i++) {
        fprintf(LOG_DEST, "\t%.2f", vlib->Wdmax[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\talbedo        :");
    for (i = 0; i < MONTHS_PER_YEAR; i++) {
        fprintf(LOG_DEST, "\t%.2f", vlib->albedo[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tdisplacement  :");
    for (i = 0; i < MONTHS_PER_YEAR; i++) {
        fprintf(LOG_DEST, "\t%.2f", vlib->displacement[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\temissivity    :");
    for (i = 0; i < MONTHS_PER_YEAR; i++) {
        fprintf(LOG_DEST, "\t%.2f", vlib->emissivity[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tNVegLibTypes  : %zu\n", vlib->NVegLibTypes);
    fprintf(LOG_DEST, "\trad_atten     : %.4f\n", vlib->rad_atten);
    fprintf(LOG_DEST, "\trarc          : %.4f\n", vlib->rarc);
    fprintf(LOG_DEST, "\trmin          : %.4f\n", vlib->rmin);
    fprintf(LOG_DEST, "\troughness     :");
    for (i = 0; i < MONTHS_PER_YEAR; i++) {
        fprintf(LOG_DEST, "\t%.2f", vlib->roughness[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\ttrunk_ratio   : %.4f\n", vlib->trunk_ratio);
    fprintf(LOG_DEST, "\twind_atten    : %.4f\n", vlib->wind_atten);
    fprintf(LOG_DEST, "\twind_h        : %.4f\n", vlib->wind_h);
    fprintf(LOG_DEST, "\tRGL           : %.4f\n", vlib->RGL);
    fprintf(LOG_DEST, "\tveg_class     : %d\n", vlib->veg_class);
    if (carbon) {
        fprintf(LOG_DEST, "\tCtype         : %d\n", vlib->Ctype);
        fprintf(LOG_DEST, "\tMaxCarboxRate : %.4f\n", vlib->MaxCarboxRate);
        fprintf(LOG_DEST, "\tMaxETransport : %.4f\n", vlib->MaxETransport);
        fprintf(LOG_DEST, "\tCO2Specificity: %.4f\n", vlib->CO2Specificity);
        fprintf(LOG_DEST, "\tLightUseEff   : %.4f\n", vlib->LightUseEff);
        fprintf(LOG_DEST, "\tNscaleFlag    : %d\n", vlib->NscaleFlag);
        fprintf(LOG_DEST, "\tWnpp_inhib    : %.4f\n", vlib->Wnpp_inhib);
        fprintf(LOG_DEST, "\tNPPfactor_sat : %.4f\n", vlib->NPPfactor_sat);
    }
}

/******************************************************************************
 * @brief    Print vegetation variables.
 *****************************************************************************/
void
print_veg_var(veg_var_struct *vvar,
              size_t          ncanopy)
{
    size_t i;

    fprintf(LOG_DEST, "veg_var:\n");
    fprintf(LOG_DEST, "\tcanopyevap   : %.4f\n", vvar->canopyevap);
    fprintf(LOG_DEST, "\tthroughfall  : %.4f\n", vvar->throughfall);
    fprintf(LOG_DEST, "\tWdew         : %.4f\n", vvar->Wdew);
    fprintf(LOG_DEST, "\tNscaleFactor :");
    for (i = 0; i < ncanopy; i++) {
        fprintf(LOG_DEST, "\t%.4f", vvar->NscaleFactor[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\taPARLayer    :");
    for (i = 0; i < ncanopy; i++) {
        fprintf(LOG_DEST, "\t%.4f", vvar->aPARLayer[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tCiLayer      :");
    for (i = 0; i < ncanopy; i++) {
        fprintf(LOG_DEST, "\t%.4f", vvar->CiLayer[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\trsLayer      :");
    for (i = 0; i < ncanopy; i++) {
        fprintf(LOG_DEST, "\t%.4f", vvar->rsLayer[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\taPAR         : %.4f\n", vvar->aPAR);
    fprintf(LOG_DEST, "\tCi           : %.4f\n", vvar->Ci);
    fprintf(LOG_DEST, "\trc           : %.4f\n", vvar->rc);
    fprintf(LOG_DEST, "\tNPPfactor    : %.4f\n", vvar->NPPfactor);
    fprintf(LOG_DEST, "\tGPP          : %.4f\n", vvar->GPP);
    fprintf(LOG_DEST, "\tRphoto       : %.4f\n", vvar->Rphoto);
    fprintf(LOG_DEST, "\tRdark        : %.4f\n", vvar->Rdark);
    fprintf(LOG_DEST, "\tRmaint       : %.4f\n", vvar->Rmaint);
    fprintf(LOG_DEST, "\tRgrowth      : %.4f\n", vvar->Rgrowth);
    fprintf(LOG_DEST, "\tRaut         : %.4f\n", vvar->Raut);
    fprintf(LOG_DEST, "\tNPP          : %.4f\n", vvar->NPP);
    fprintf(LOG_DEST, "\tLitterfall   : %.4f\n", vvar->Litterfall);
    fprintf(LOG_DEST, "\tAnnualNPP    : %.4f\n", vvar->AnnualNPP);
    fprintf(LOG_DEST, "\tAnnualNPPPrev: %.4f\n", vvar->AnnualNPPPrev);
}
