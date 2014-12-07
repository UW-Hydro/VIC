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
                size_t            nfrost,
                size_t            npet)
{
    size_t i;

    fprintf(LOG_DEST, "cell_data:\n");
    fprintf(LOG_DEST, "\taero_resist :");
    for (i = 0; i < 2; i++) {
        fprintf(LOG_DEST, "\t%.4lf", cell->aero_resist[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tasat        : %.4lf\n", cell->asat);
    fprintf(LOG_DEST, "\tbaseflow    : %.4lf\n", cell->baseflow);
    fprintf(LOG_DEST, "\tCLitter     : %.4lf\n", cell->CLitter);
    fprintf(LOG_DEST, "\tCInter      : %.4lf\n", cell->CInter);
    fprintf(LOG_DEST, "\tCSlow       : %.4lf\n", cell->CSlow);
    fprintf(LOG_DEST, "\tinflow      : %.4lf\n", cell->inflow);
    fprintf(LOG_DEST, "\tpot_evap    :");
    for (i = 0; i < npet; i++) {
        fprintf(LOG_DEST, "\t%.4lf", cell->pot_evap[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\trunoff      : %.4lf\n", cell->runoff);
    for (i = 0; i < nlayers; i++) {
        fprintf(LOG_DEST, "\tlayer %zd   :\n", i);
        print_layer_data(&(cell->layer[i]), nfrost);
    }
    fprintf(LOG_DEST, "\tRhLitter    : %.4lf\n", cell->RhLitter);
    fprintf(LOG_DEST, "\tRhLitter2Atm: %.4lf\n", cell->RhLitter2Atm);
    fprintf(LOG_DEST, "\tRhInter     : %.4lf\n", cell->RhInter);
    fprintf(LOG_DEST, "\tRhSlow      : %.4lf\n", cell->RhSlow);
    fprintf(LOG_DEST, "\tRhTot       : %.4lf\n", cell->RhTot);
    fprintf(LOG_DEST, "\trootmoist   : %.4lf\n", cell->rootmoist);
    fprintf(LOG_DEST, "\twetness     : %.4lf\n", cell->wetness);
    fprintf(LOG_DEST, "\tzwt         : %.4lf\n", cell->zwt);
    fprintf(LOG_DEST, "\tzwt_lumped  : %.4lf\n", cell->zwt_lumped);
}

/******************************************************************************
 * @brief    Print day-month-year structure.
 *****************************************************************************/
void
print_dmy(dmy_struct *dmy)
{
    fprintf(LOG_DEST, "dmy:\n");
    fprintf(LOG_DEST, "\tday        : %d\n", dmy->day);
    fprintf(LOG_DEST, "\tday_in_year: %d\n", dmy->day_in_year);
    fprintf(LOG_DEST, "\thour       : %d\n", dmy->hour);
    fprintf(LOG_DEST, "\tmonth      : %d\n", dmy->month);
    fprintf(LOG_DEST, "\tyear       : %d\n", dmy->year);
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
    fprintf(LOG_DEST, "\tAlbedoLake       : %.4lf\n", eb->AlbedoLake);
    fprintf(LOG_DEST, "\tAlbedoOver       : %.4lf\n", eb->AlbedoOver);
    fprintf(LOG_DEST, "\tAlbedoUnder      : %.4lf\n", eb->AlbedoUnder);
    fprintf(LOG_DEST, "\tCs               :");
    for (i = 0; i < 2; i++) {
        fprintf(LOG_DEST, "\t%.4lf", eb->Cs[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tCs_node          :");
    for (i = 0; i < nnodes; i++) {
        fprintf(LOG_DEST, "\t%.4lf", eb->Cs_node[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tfdepth           :");
    for (i = 0; i < nfronts; i++) {
        fprintf(LOG_DEST, "\t%.4lf", eb->fdepth[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tfrozen           : %d\n", eb->frozen);
    fprintf(LOG_DEST, "\tice              :");
    for (i = 0; i < nnodes; i++) {
        fprintf(LOG_DEST, "\t%.4lf", eb->ice[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tkappa            :");
    for (i = 0; i < 2; i++) {
        fprintf(LOG_DEST, "\t%.4lf", eb->kappa[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tkappa_node       :");
    for (i = 0; i < nnodes; i++) {
        fprintf(LOG_DEST, "\t%.4lf", eb->kappa_node[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tmoist            :");
    for (i = 0; i < nnodes; i++) {
        fprintf(LOG_DEST, "\t%.4lf", eb->moist[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tNfrost           : %zu\n", eb->Nfrost);
    fprintf(LOG_DEST, "\tNthaw            : %zu\n", eb->Nthaw);
    fprintf(LOG_DEST, "\tT                :");
    for (i = 0; i < nnodes; i++) {
        fprintf(LOG_DEST, "\t%.4lf", eb->T[i]);
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
    fprintf(LOG_DEST, "\tTcanopy          : %.4lf\n", eb->Tcanopy);
    fprintf(LOG_DEST, "\tTcanopy_fbflag   : %d\n", eb->Tcanopy_fbflag);
    fprintf(LOG_DEST, "\tTcanopy_fbcount  : %d\n", eb->Tcanopy_fbcount);
    fprintf(LOG_DEST, "\ttdepth           :");
    for (i = 0; i < nfronts; i++) {
        fprintf(LOG_DEST, "\t%.4lf", eb->tdepth[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tTfoliage         : %.4lf\n", eb->Tfoliage);
    fprintf(LOG_DEST, "\tTfoliage_fbflag  : %d\n", eb->Tfoliage_fbflag);
    fprintf(LOG_DEST, "\tTfoliage_fbcount : %d\n", eb->Tfoliage_fbcount);
    fprintf(LOG_DEST, "\tTsurf            : %.4lf\n", eb->Tsurf);
    fprintf(LOG_DEST, "\tTsurf_fbflag     : %d\n", eb->Tsurf_fbflag);
    fprintf(LOG_DEST, "\tTsurf_fbcount    : %d\n", eb->Tsurf_fbcount);
    fprintf(LOG_DEST, "\tunfrozen         : %.4lf\n", eb->unfrozen);
    fprintf(LOG_DEST, "\tadvected_sensible: %.4lf\n", eb->advected_sensible);
    fprintf(LOG_DEST, "\tadvection        : %.4lf\n", eb->advection);
    fprintf(LOG_DEST, "\tAtmosError       : %.4lf\n", eb->AtmosError);
    fprintf(LOG_DEST, "\tAtmosLatent      : %.4lf\n", eb->AtmosLatent);
    fprintf(LOG_DEST, "\tAtmosLatentSub   : %.4lf\n", eb->AtmosLatentSub);
    fprintf(LOG_DEST, "\tAtmosSensible    : %.4lf\n", eb->AtmosSensible);
    fprintf(LOG_DEST, "\tcanopy_advection : %.4lf\n", eb->canopy_advection);
    fprintf(LOG_DEST, "\tcanopy_latent    : %.4lf\n", eb->canopy_latent);
    fprintf(LOG_DEST, "\tcanopy_latent_sub: %.4lf\n", eb->canopy_latent_sub);
    fprintf(LOG_DEST, "\tcanopy_refreeze  : %.4lf\n", eb->canopy_refreeze);
    fprintf(LOG_DEST, "\tcanopy_sensible  : %.4lf\n", eb->canopy_sensible);
    fprintf(LOG_DEST, "\tdeltaCC          : %.4lf\n", eb->deltaCC);
    fprintf(LOG_DEST, "\tdeltaH           : %.4lf\n", eb->deltaH);
    fprintf(LOG_DEST, "\terror            : %.4lf\n", eb->error);
    fprintf(LOG_DEST, "\tfusion           : %.4lf\n", eb->fusion);
    fprintf(LOG_DEST, "\tgrnd_flux        : %.4lf\n", eb->grnd_flux);
    fprintf(LOG_DEST, "\tlatent           : %.4lf\n", eb->latent);
    fprintf(LOG_DEST, "\tlatent_sub       : %.4lf\n", eb->latent_sub);
    fprintf(LOG_DEST, "\tlongwave         : %.4lf\n", eb->longwave);
    fprintf(LOG_DEST, "\tLongOverIn       : %.4lf\n", eb->LongOverIn);
    fprintf(LOG_DEST, "\tLongUnderIn      : %.4lf\n", eb->LongUnderIn);
    fprintf(LOG_DEST, "\tLongUnderOut     : %.4lf\n", eb->LongUnderOut);
    fprintf(LOG_DEST, "\tmelt_energy      : %.4lf\n", eb->melt_energy);
    fprintf(LOG_DEST, "\tNetLongAtmos     : %.4lf\n", eb->NetLongAtmos);
    fprintf(LOG_DEST, "\tNetLongOver      : %.4lf\n", eb->NetLongOver);
    fprintf(LOG_DEST, "\tNetLongUnder     : %.4lf\n", eb->NetLongUnder);
    fprintf(LOG_DEST, "\tNetShortAtmos    : %.4lf\n", eb->NetShortAtmos);
    fprintf(LOG_DEST, "\tNetShortGrnd     : %.4lf\n", eb->NetShortGrnd);
    fprintf(LOG_DEST, "\tNetShortOver     : %.4lf\n", eb->NetShortOver);
    fprintf(LOG_DEST, "\tNetShortUnder    : %.4lf\n", eb->NetShortUnder);
    fprintf(LOG_DEST, "\tout_long_canopy  : %.4lf\n", eb->out_long_canopy);
    fprintf(LOG_DEST, "\tout_long_surface : %.4lf\n", eb->out_long_surface);
    fprintf(LOG_DEST, "\trefreeze_energy  : %.4lf\n", eb->refreeze_energy);
    fprintf(LOG_DEST, "\tsensible         : %.4lf\n", eb->sensible);
    fprintf(LOG_DEST, "\tshortwave        : %.4lf\n", eb->shortwave);
    fprintf(LOG_DEST, "\tShortOverIn      : %.4lf\n", eb->ShortOverIn);
    fprintf(LOG_DEST, "\tShortUnderIn     : %.4lf\n", eb->ShortUnderIn);
    fprintf(LOG_DEST, "\tsnow_flux        : %.4lf\n", eb->snow_flux);
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
    fprintf(LOG_DEST, "\tmultiplier: %lf\n", force_type->multiplier);
}

/******************************************************************************
 * @brief    Print global parameters structure.
 *****************************************************************************/
void
print_global_param(global_param_struct *gp)
{
    size_t i;

    fprintf(LOG_DEST, "global_param:\n");
    fprintf(LOG_DEST, "\tmeasure_h    : %.4lf\n", gp->measure_h);
    fprintf(LOG_DEST, "\twind_h       : %.4lf\n", gp->wind_h);
    fprintf(LOG_DEST, "\tresolution   : %.4f\n", gp->resolution);
    fprintf(LOG_DEST, "\tdt           : %d\n", gp->dt);
    fprintf(LOG_DEST, "\tout_dt       : %d\n", gp->out_dt);
    fprintf(LOG_DEST, "\tendday       : %d\n", gp->endday);
    fprintf(LOG_DEST, "\tendmonth     : %d\n", gp->endmonth);
    fprintf(LOG_DEST, "\tendyear      : %d\n", gp->endyear);
    for (i = 0; i < 2; i++) {
        fprintf(LOG_DEST, "\tforceday[%zd]   : %d\n", i, gp->forceday[i]);
        fprintf(LOG_DEST, "\tforcehour[%zd]  : %d\n", i, gp->forcehour[i]);
        fprintf(LOG_DEST, "\tforcemonth[%zd] : %d\n", i, gp->forcemonth[i]);
        fprintf(LOG_DEST, "\tforceoffset[%zd]: %d\n", i, gp->forceoffset[i]);
        fprintf(LOG_DEST, "\tforceskip[%zd]  : %d\n", i, gp->forceskip[i]);
        fprintf(LOG_DEST, "\tforceyear[%zd]  : %d\n", i, gp->forceyear[i]);
    }
    fprintf(LOG_DEST, "\tnrecs        : %d\n", gp->nrecs);
    fprintf(LOG_DEST, "\tskipyear     : %d\n", gp->skipyear);
    fprintf(LOG_DEST, "\tstartday     : %d\n", gp->startday);
    fprintf(LOG_DEST, "\tstarthour    : %d\n", gp->starthour);
    fprintf(LOG_DEST, "\tstartmonth   : %d\n", gp->startmonth);
    fprintf(LOG_DEST, "\tstartyear    : %d\n", gp->startyear);
    fprintf(LOG_DEST, "\tstateday     : %d\n", gp->stateday);
    fprintf(LOG_DEST, "\tstatemonth   : %d\n", gp->statemonth);
    fprintf(LOG_DEST, "\tstateyear    : %d\n", gp->stateyear);
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
        fprintf(LOG_DEST, "\t%.4lf", lcon->z[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tbasin    :");
    for (i = 0; i < nlnodes; i++) {
        fprintf(LOG_DEST, "\t%.4lf", lcon->basin[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tCl       :");
    for (i = 0; i < nlnodes; i++) {
        fprintf(LOG_DEST, "\t%.4lf", lcon->Cl[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tb        : %.4lf\n", lcon->b);
    fprintf(LOG_DEST, "\tmaxdepth : %.4lf\n", lcon->maxdepth);
    fprintf(LOG_DEST, "\tmindepth : %.4lf\n", lcon->mindepth);
    fprintf(LOG_DEST, "\tmaxvolume: %.4lf\n", lcon->maxvolume);
    fprintf(LOG_DEST, "\tminvolume: %.4lf\n", lcon->minvolume);
    fprintf(LOG_DEST, "\tbpercent : %.4f\n", lcon->bpercent);
    fprintf(LOG_DEST, "\trpercent : %.4f\n", lcon->rpercent);
    fprintf(LOG_DEST, "\twfrac    : %.4lf\n", lcon->wfrac);
    fprintf(LOG_DEST, "\tdepth_in : %.4lf\n", lcon->depth_in);
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
               size_t           nfrost,
               size_t           npet)
{
    size_t i;

    fprintf(LOG_DEST, "lake_var:\n");
    fprintf(LOG_DEST, "\tactivenod      : %d\n", lvar->activenod);
    fprintf(LOG_DEST, "\tdz             : %.4lf\n", lvar->dz);
    fprintf(LOG_DEST, "\tsurfdz         : %.4lf\n", lvar->surfdz);
    fprintf(LOG_DEST, "\tldepth         : %.4lf\n", lvar->ldepth);
    fprintf(LOG_DEST, "\tsurface        :");
    for (i = 0; i < nlnodes + 1; i++) {
        fprintf(LOG_DEST, "\t%.4lf", lvar->surface[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tsarea          : %.4lf\n", lvar->sarea);
    fprintf(LOG_DEST, "\tsarea_save     : %.4lf\n", lvar->sarea_save);
    fprintf(LOG_DEST, "\tvolume         : %.4lf\n", lvar->volume);
    fprintf(LOG_DEST, "\tvolume_save    : %.4lf\n", lvar->volume_save);
    fprintf(LOG_DEST, "\ttemp           :");
    for (i = 0; i < nlnodes; i++) {
        fprintf(LOG_DEST, "\t%.4lf", lvar->temp[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\ttempavg        : %.4lf\n", lvar->tempavg);
    fprintf(LOG_DEST, "\tareai          : %.4lf\n", lvar->areai);
    fprintf(LOG_DEST, "\tnew_ice_area   : %.4lf\n", lvar->new_ice_area);
    fprintf(LOG_DEST, "\tice_water_eq   : %.4lf\n", lvar->ice_water_eq);
    fprintf(LOG_DEST, "\thice           : %.4lf\n", lvar->hice);
    fprintf(LOG_DEST, "\ttempi          : %.4lf\n", lvar->tempi);
    fprintf(LOG_DEST, "\tswe            : %.4lf\n", lvar->swe);
    fprintf(LOG_DEST, "\tswe_save       : %.4lf\n", lvar->swe_save);
    fprintf(LOG_DEST, "\tsurf_temp      : %.4lf\n", lvar->surf_temp);
    fprintf(LOG_DEST, "\tpack_temp      : %.4lf\n", lvar->pack_temp);
    fprintf(LOG_DEST, "\tcoldcontent    : %.4lf\n", lvar->coldcontent);
    fprintf(LOG_DEST, "\tsurf_water     : %.4lf\n", lvar->surf_water);
    fprintf(LOG_DEST, "\tpack_water     : %.4lf\n", lvar->pack_water);
    fprintf(LOG_DEST, "\tSAlbedo        : %.4lf\n", lvar->SAlbedo);
    fprintf(LOG_DEST, "\tsdepth         : %.4lf\n", lvar->sdepth);
    fprintf(LOG_DEST, "\taero_resist    : %.4lf\n", lvar->aero_resist);
    fprintf(LOG_DEST, "\tdensity        :");
    for (i = 0; i < nlnodes; i++) {
        fprintf(LOG_DEST, "\t%.4lf", lvar->density[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tbaseflow_in    : %.4lf\n", lvar->baseflow_in);
    fprintf(LOG_DEST, "\tbaseflow_out   : %.4lf\n", lvar->baseflow_out);
    fprintf(LOG_DEST, "\tchannel_in     : %.4lf\n", lvar->channel_in);
    fprintf(LOG_DEST, "\tevapw          : %.4lf\n", lvar->evapw);
    fprintf(LOG_DEST, "\tice_throughfall: %.4lf\n", lvar->ice_throughfall);
    fprintf(LOG_DEST, "\tprec           : %.4lf\n", lvar->prec);
    fprintf(LOG_DEST, "\trecharge       : %.4lf\n", lvar->recharge);
    fprintf(LOG_DEST, "\trunoff_in      : %.4lf\n", lvar->runoff_in);
    fprintf(LOG_DEST, "\trunoff_out     : %.4lf\n", lvar->runoff_out);
    fprintf(LOG_DEST, "\tsnowmlt        : %.4lf\n", lvar->snowmlt);
    fprintf(LOG_DEST, "\tvapor_flux     : %.4lf\n", lvar->vapor_flux);
    print_snow_data(&(lvar->snow));
    print_energy_bal(&(lvar->energy), nnodes, nfronts);
    print_cell_data(&(lvar->soil), nlayers, nfrost, npet);
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
    fprintf(LOG_DEST, "\tCs   : %.4lf\n", ldata->Cs);
    fprintf(LOG_DEST, "\tT    : %.4lf\n", ldata->T);
    fprintf(LOG_DEST, "\tevap : %.4lf\n", ldata->evap);
    fprintf(LOG_DEST, "\tice  :");
    for (i = 0; i < nfrost; i++) {
        fprintf(LOG_DEST, "\t%.4lf", ldata->ice[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tkappa: %.4lf\n", ldata->kappa);
    fprintf(LOG_DEST, "\tmoist: %.4lf\n", ldata->moist);
    fprintf(LOG_DEST, "\tphi  : %.4lf\n", ldata->phi);
    fprintf(LOG_DEST, "\tzwt  : %.4lf\n", ldata->zwt);
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
    fprintf(LOG_DEST, "\tLW_CLOUD             : %d\n", option->LW_CLOUD);
    fprintf(LOG_DEST, "\tLW_TYPE              : %d\n", option->LW_TYPE);
    fprintf(LOG_DEST, "\tMTCLIM_SWE_CORR      : %d\n", option->MTCLIM_SWE_CORR);
    fprintf(LOG_DEST, "\tNcanopy              : %zu\n", option->Ncanopy);
    fprintf(LOG_DEST, "\tNfrost               : %zu\n", option->Nfrost);
    fprintf(LOG_DEST, "\tNlakenode            : %zu\n", option->Nlakenode);
    fprintf(LOG_DEST, "\tNlayer               : %zu\n", option->Nlayer);
    fprintf(LOG_DEST, "\tNnode                : %zu\n", option->Nnode);
    fprintf(LOG_DEST, "\tNOFLUX               : %d\n", option->NOFLUX);
    fprintf(LOG_DEST, "\tNVEGTYPES            : %zu\n", option->NVEGTYPES);
    fprintf(LOG_DEST, "\tPLAPSE               : %d\n", option->PLAPSE);
    fprintf(LOG_DEST, "\tRC_MODE              : %d\n", option->RC_MODE);
    fprintf(LOG_DEST, "\tROOT_ZONES           : %zu\n", option->ROOT_ZONES);
    fprintf(LOG_DEST, "\tQUICK_FLUX           : %d\n", option->QUICK_FLUX);
    fprintf(LOG_DEST, "\tQUICK_SOLVE          : %d\n", option->QUICK_SOLVE);
    fprintf(LOG_DEST, "\tSHARE_LAYER_MOIST    : %d\n",
            option->SHARE_LAYER_MOIST);
    fprintf(LOG_DEST, "\tSNOW_DENSITY         : %d\n", option->SNOW_DENSITY);
    fprintf(LOG_DEST, "\tSNOW_BAND            : %zu\n", option->SNOW_BAND);
    fprintf(LOG_DEST, "\tSNOW_STEP            : %d\n", option->SNOW_STEP);
    fprintf(LOG_DEST, "\tSPATIAL_FROST        : %d\n", option->SPATIAL_FROST);
    fprintf(LOG_DEST, "\tSPATIAL_SNOW         : %d\n", option->SPATIAL_SNOW);
    fprintf(LOG_DEST, "\tTFALLBACK            : %d\n", option->TFALLBACK);
    fprintf(LOG_DEST, "\tVP_INTERP            : %d\n", option->VP_INTERP);
    fprintf(LOG_DEST, "\tVP_ITER              : %d\n", option->VP_ITER);
    fprintf(LOG_DEST, "\tALMA_INPUT           : %d\n", option->ALMA_INPUT);
    fprintf(LOG_DEST, "\tBASEFLOW             : %d\n", option->BASEFLOW);
    fprintf(LOG_DEST, "\tGRID_DECIMAL         : %d\n", option->GRID_DECIMAL);
    fprintf(LOG_DEST, "\tVEGLIB_PHOTO         : %d\n", option->VEGLIB_PHOTO);
    fprintf(LOG_DEST, "\tVEGPARAM_LAI         : %d\n", option->VEGPARAM_LAI);
    fprintf(LOG_DEST, "\tLAI_SRC              : %d\n", option->LAI_SRC);
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
    fprintf(LOG_DEST, "\tOUTPUT_FORCE         : %d\n", option->OUTPUT_FORCE);
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
        fprintf(LOG_DEST, "\t%.4lf", out->data[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\taggdata:");
    for (i = 0; i < nelem; i++) {
        fprintf(LOG_DEST, "\t%.4lf", out->aggdata[i]);
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
    fprintf(LOG_DEST, "\tFORCE_DT    : %d %d\n", param_set->FORCE_DT[0],
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
    fprintf(LOG_DEST, "\tN_TYPES     : %d %d\n", param_set->N_TYPES[0],
            param_set->N_TYPES[1]);
}

/******************************************************************************
 * @brief    Print model parameters.
 *****************************************************************************/
void
print_parameters(parameters_struct *param)
{
    fprintf(LOG_DEST, "parameters:\n");
    fprintf(LOG_DEST, "\tLAPSE_RATE: %.4lf\n", param->LAPSE_RATE);
    fprintf(LOG_DEST, "\tGAUGE_HEIGHT: %.4lf\n", param->GAUGE_HEIGHT);
    fprintf(LOG_DEST, "\tWIND_SPEED_DEFAULT: %.4lf\n",
            param->WIND_SPEED_DEFAULT);
    fprintf(LOG_DEST, "\tWIND_SPEED_MIN: %.4lf\n", param->WIND_SPEED_MIN);
    fprintf(LOG_DEST, "\tHUGE_RESIST: %.4lf\n", param->HUGE_RESIST);
    fprintf(LOG_DEST, "\tALBEDO_BARE_SOIL: %.4lf\n", param->ALBEDO_BARE_SOIL);
    fprintf(LOG_DEST, "\tALBEDO_H20_SURF: %.4lf\n", param->ALBEDO_H20_SURF);
    fprintf(LOG_DEST, "\tEMISS_GRND: %.4lf\n", param->EMISS_GRND);
    fprintf(LOG_DEST, "\tEMISS_VEG: %.4lf\n", param->EMISS_VEG);
    fprintf(LOG_DEST, "\tEMISS_ICE: %.4lf\n", param->EMISS_ICE);
    fprintf(LOG_DEST, "\tEMISS_SNOW: %.4lf\n", param->EMISS_SNOW);
    fprintf(LOG_DEST, "\tEMISS_H2O: %.4lf\n", param->EMISS_H2O);
    fprintf(LOG_DEST, "\tSOIL_RESID_MOIST: %.4lf\n", param->SOIL_RESID_MOIST);
    fprintf(LOG_DEST, "\tSOIL_SLAB_MOIST_FRACT: %.4lf\n",
            param->SOIL_SLAB_MOIST_FRACT);
    fprintf(LOG_DEST, "\tVEG_LAI_SNOW_MULTIPLIER: %.4lf\n",
            param->VEG_LAI_SNOW_MULTIPLIER);
    fprintf(LOG_DEST, "\tVEG_MIN_INTERCEPTION_STORAGE: %.4lf\n",
            param->VEG_MIN_INTERCEPTION_STORAGE);
    fprintf(LOG_DEST, "\tVEG_LAI_WATER_FACTOR: %.4lf\n",
            param->VEG_LAI_WATER_FACTOR);
    fprintf(LOG_DEST, "\tCANOPY_CLOSURE: %.4lf\n", param->CANOPY_CLOSURE);
    fprintf(LOG_DEST, "\tCANOPY_RSMAX: %.4lf\n", param->CANOPY_RSMAX);
    fprintf(LOG_DEST, "\tCANOPY_VPDMINFACTOR: %.4lf\n",
            param->CANOPY_VPDMINFACTOR);
    fprintf(LOG_DEST, "\tMTCLIM_SOLAR_CONSTANT: %.4lf\n",
            param->MTCLIM_SOLAR_CONSTANT);
    fprintf(LOG_DEST, "\tMTCLIM_TDAYCOEF: %.4lf\n", param->MTCLIM_TDAYCOEF);
    fprintf(LOG_DEST, "\tMTCLIM_SNOW_TCRIT: %.4lf\n", param->MTCLIM_SNOW_TCRIT);
    fprintf(LOG_DEST, "\tMTCLIM_SNOW_TRATE: %.4lf\n", param->MTCLIM_SNOW_TRATE);
    fprintf(LOG_DEST, "\tMTCLIM_TBASE: %.4lf\n", param->MTCLIM_TBASE);
    fprintf(LOG_DEST, "\tMTCLIM_ABASE: %.4lf\n", param->MTCLIM_ABASE);
    fprintf(LOG_DEST, "\tMTCLIM_C: %.4lf\n", param->MTCLIM_C);
    fprintf(LOG_DEST, "\tMTCLIM_B0: %.4lf\n", param->MTCLIM_B0);
    fprintf(LOG_DEST, "\tMTCLIM_B1: %.4lf\n", param->MTCLIM_B1);
    fprintf(LOG_DEST, "\tMTCLIM_B2: %.4lf\n", param->MTCLIM_B2);
    fprintf(LOG_DEST, "\tMTCLIM_RAIN_SCALAR: %.4lf\n",
            param->MTCLIM_RAIN_SCALAR);
    fprintf(LOG_DEST, "\tMTCLIM_DIF_ALB: %.4lf\n", param->MTCLIM_DIF_ALB);
    fprintf(LOG_DEST, "\tMTCLIM_SC_INT: %.4lf\n", param->MTCLIM_SC_INT);
    fprintf(LOG_DEST, "\tMTCLIM_SC_SLOPE: %.4lf\n", param->MTCLIM_SC_SLOPE);
    fprintf(LOG_DEST, "\tMTCLIM_SRADDT: %.4lf\n", param->MTCLIM_SRADDT);
    fprintf(LOG_DEST, "\tMTCLIM_SW_PREC_THRESH: %.4lf\n",
            param->MTCLIM_SW_PREC_THRESH);
    fprintf(LOG_DEST, "\tLAKE_TMELT: %.4lf\n", param->LAKE_TMELT);
    fprintf(LOG_DEST, "\tLAKE_MAX_SURFACE: %.4lf\n", param->LAKE_MAX_SURFACE);
    fprintf(LOG_DEST, "\tLAKE_BETA: %.4lf\n", param->LAKE_BETA);
    fprintf(LOG_DEST, "\tLAKE_FRACMIN: %.4lf\n", param->LAKE_FRACMIN);
    fprintf(LOG_DEST, "\tLAKE_FRACLIM: %.4lf\n", param->LAKE_FRACLIM);
    fprintf(LOG_DEST, "\tLAKE_DM: %.4lf\n", param->LAKE_DM);
    fprintf(LOG_DEST, "\tLAKE_SNOWCRIT: %.4lf\n", param->LAKE_SNOWCRIT);
    fprintf(LOG_DEST, "\tLAKE_ZWATER: %.4lf\n", param->LAKE_ZWATER);
    fprintf(LOG_DEST, "\tLAKE_ZSNOW: %.4lf\n", param->LAKE_ZSNOW);
    fprintf(LOG_DEST, "\tLAKE_RHOSNOW: %.4lf\n", param->LAKE_RHOSNOW);
    fprintf(LOG_DEST, "\tLAKE_CONDI: %.4lf\n", param->LAKE_CONDI);
    fprintf(LOG_DEST, "\tLAKE_CONDS: %.4lf\n", param->LAKE_CONDS);
    fprintf(LOG_DEST, "\tLAKE_LAMISW: %.4lf\n", param->LAKE_LAMISW);
    fprintf(LOG_DEST, "\tLAKE_LAMILW: %.4lf\n", param->LAKE_LAMILW);
    fprintf(LOG_DEST, "\tLAKE_LAMSSW: %.4lf\n", param->LAKE_LAMSSW);
    fprintf(LOG_DEST, "\tLAKE_LAMSLW: %.4lf\n", param->LAKE_LAMSLW);
    fprintf(LOG_DEST, "\tLAKE_LAMWSW: %.4lf\n", param->LAKE_LAMWSW);
    fprintf(LOG_DEST, "\tLAKE_LAMWLW: %.4lf\n", param->LAKE_LAMWLW);
    fprintf(LOG_DEST, "\tLAKE_A1: %.4lf\n", param->LAKE_A1);
    fprintf(LOG_DEST, "\tLAKE_A2: %.4lf\n", param->LAKE_A2);
    fprintf(LOG_DEST, "\tLAKE_QWTAU: %.4lf\n", param->LAKE_QWTAU);
    fprintf(LOG_DEST, "\tLAKE_MAX_ITER: %d\n", param->LAKE_MAX_ITER);
    fprintf(LOG_DEST, "\tSVP_A: %.4lf\n", param->SVP_A);
    fprintf(LOG_DEST, "\tSVP_B: %.4lf\n", param->SVP_B);
    fprintf(LOG_DEST, "\tSVP_C: %.4lf\n", param->SVP_C);
    fprintf(LOG_DEST, "\tCARBON_CATMCURRENT: %.4lf\n",
            param->CARBON_CATMCURRENT);
    fprintf(LOG_DEST, "\tCARBON_SW2PAR: %.4lf\n", param->CARBON_SW2PAR);
    fprintf(LOG_DEST, "\tPHOTO_OMEGA: %.4lf\n", param->PHOTO_OMEGA);
    fprintf(LOG_DEST, "\tPHOTO_LAIMAX: %.4lf\n", param->PHOTO_LAIMAX);
    fprintf(LOG_DEST, "\tPHOTO_LAILIMIT: %.4lf\n", param->PHOTO_LAILIMIT);
    fprintf(LOG_DEST, "\tPHOTO_LAIMIN: %.4lf\n", param->PHOTO_LAIMIN);
    fprintf(LOG_DEST, "\tPHOTO_EPAR: %.4lf\n", param->PHOTO_EPAR);
    fprintf(LOG_DEST, "\tPHOTO_FCMAX: %.4lf\n", param->PHOTO_FCMAX);
    fprintf(LOG_DEST, "\tPHOTO_FCMIN: %.4lf\n", param->PHOTO_FCMIN);
    fprintf(LOG_DEST, "\tPHOTO_ZENITHMIN: %.4lf\n", param->PHOTO_ZENITHMIN);
    fprintf(LOG_DEST, "\tPHOTO_ZENITHMINPAR: %.4lf\n",
            param->PHOTO_ZENITHMINPAR);
    fprintf(LOG_DEST, "\tPHOTO_ALBSOIPARMIN: %.4lf\n",
            param->PHOTO_ALBSOIPARMIN);
    fprintf(LOG_DEST, "\tPHOTO_MINMAXETRANS: %.4lf\n",
            param->PHOTO_MINMAXETRANS);
    fprintf(LOG_DEST, "\tPHOTO_MINSTOMCOND: %.4lf\n", param->PHOTO_MINSTOMCOND);
    fprintf(LOG_DEST, "\tPHOTO_FCI1C3: %.4lf\n", param->PHOTO_FCI1C3);
    fprintf(LOG_DEST, "\tPHOTO_FCI1C4: %.4lf\n", param->PHOTO_FCI1C4);
    fprintf(LOG_DEST, "\tPHOTO_OX: %.4lf\n", param->PHOTO_OX);
    fprintf(LOG_DEST, "\tPHOTO_KC: %.4lf\n", param->PHOTO_KC);
    fprintf(LOG_DEST, "\tPHOTO_KO: %.4lf\n", param->PHOTO_KO);
    fprintf(LOG_DEST, "\tPHOTO_EC: %.4lf\n", param->PHOTO_EC);
    fprintf(LOG_DEST, "\tPHOTO_EO: %.4lf\n", param->PHOTO_EO);
    fprintf(LOG_DEST, "\tPHOTO_EV: %.4lf\n", param->PHOTO_EV);
    fprintf(LOG_DEST, "\tPHOTO_ER: %.4lf\n", param->PHOTO_ER);
    fprintf(LOG_DEST, "\tPHOTO_ALC3: %.4lf\n", param->PHOTO_ALC3);
    fprintf(LOG_DEST, "\tPHOTO_FRDC3: %.4lf\n", param->PHOTO_FRDC3);
    fprintf(LOG_DEST, "\tPHOTO_EK: %.4lf\n", param->PHOTO_EK);
    fprintf(LOG_DEST, "\tPHOTO_ALC4: %.4lf\n", param->PHOTO_ALC4);
    fprintf(LOG_DEST, "\tPHOTO_FRDC4: %.4lf\n", param->PHOTO_FRDC4);
    fprintf(LOG_DEST, "\tPHOTO_THETA: %.4lf\n", param->PHOTO_THETA);
    fprintf(LOG_DEST, "\tPHOTO_FRLEAF: %.4lf\n", param->PHOTO_FRLEAF);
    fprintf(LOG_DEST, "\tPHOTO_FRGROWTH: %.4lf\n", param->PHOTO_FRGROWTH);
    fprintf(LOG_DEST, "\tSRESP_E0_LT: %.4lf\n", param->SRESP_E0_LT);
    fprintf(LOG_DEST, "\tSRESP_T0_LT: %.4lf\n", param->SRESP_T0_LT);
    fprintf(LOG_DEST, "\tSRESP_WMINFM: %.4lf\n", param->SRESP_WMINFM);
    fprintf(LOG_DEST, "\tSRESP_WMAXFM: %.4lf\n", param->SRESP_WMAXFM);
    fprintf(LOG_DEST, "\tSRESP_WOPTFM: %.4lf\n", param->SRESP_WOPTFM);
    fprintf(LOG_DEST, "\tSRESP_RHSAT: %.4lf\n", param->SRESP_RHSAT);
    fprintf(LOG_DEST, "\tSRESP_RFACTOR: %.4lf\n", param->SRESP_RFACTOR);
    fprintf(LOG_DEST, "\tSRESP_TAULITTER: %.4lf\n", param->SRESP_TAULITTER);
    fprintf(LOG_DEST, "\tSRESP_TAUINTER: %.4lf\n", param->SRESP_TAUINTER);
    fprintf(LOG_DEST, "\tSRESP_TAUSLOW: %.4lf\n", param->SRESP_TAUSLOW);
    fprintf(LOG_DEST, "\tSRESP_FAIR: %.4lf\n", param->SRESP_FAIR);
    fprintf(LOG_DEST, "\tSRESP_FINTER: %.4lf\n", param->SRESP_FINTER);
    fprintf(LOG_DEST, "\tSNOW_MAX_SURFACE_SWE: %.4lf\n",
            param->SNOW_MAX_SURFACE_SWE);
    fprintf(LOG_DEST, "\tSNOW_LIQUID_WATER_CAPACITY: %.4lf\n",
            param->SNOW_LIQUID_WATER_CAPACITY);
    fprintf(LOG_DEST, "\tSNOW_NEW_SNOW_DENSITY: %.4lf\n",
            param->SNOW_NEW_SNOW_DENSITY);
    fprintf(LOG_DEST, "\tSNOW_DENS_DMLIMIT: %.4lf\n", param->SNOW_DENS_DMLIMIT);
    fprintf(LOG_DEST, "\tSNOW_DENS_MAX_CHANGE: %.4lf\n",
            param->SNOW_DENS_MAX_CHANGE);
    fprintf(LOG_DEST, "\tSNOW_DENS_ETA0: %.4lf\n", param->SNOW_DENS_ETA0);
    fprintf(LOG_DEST, "\tSNOW_DENS_C1: %.4lf\n", param->SNOW_DENS_C1);
    fprintf(LOG_DEST, "\tSNOW_DENS_C2: %.4lf\n", param->SNOW_DENS_C2);
    fprintf(LOG_DEST, "\tSNOW_DENS_C5: %.4lf\n", param->SNOW_DENS_C5);
    fprintf(LOG_DEST, "\tSNOW_DENS_C6: %.4lf\n", param->SNOW_DENS_C6);
    fprintf(LOG_DEST, "\tSNOW_DENS_F: %.4lf\n", param->SNOW_DENS_F);
    fprintf(LOG_DEST, "\tSNOW_MIN_SWQ_EB_THRES: %.4lf\n",
            param->SNOW_MIN_SWQ_EB_THRES);
    fprintf(LOG_DEST, "\tSNOW_A1: %.4lf\n", param->SNOW_A1);
    fprintf(LOG_DEST, "\tSNOW_A2: %.4lf\n", param->SNOW_A2);
    fprintf(LOG_DEST, "\tSNOW_L1: %.4lf\n", param->SNOW_L1);
    fprintf(LOG_DEST, "\tSNOW_L2: %.4lf\n", param->SNOW_L2);
    fprintf(LOG_DEST, "\tSNOW_NEW_SNOW_ALB: %.4lf\n", param->SNOW_NEW_SNOW_ALB);
    fprintf(LOG_DEST, "\tSNOW_ALB_ACCUM_A: %.4lf\n", param->SNOW_ALB_ACCUM_A);
    fprintf(LOG_DEST, "\tSNOW_ALB_ACCUM_B: %.4lf\n", param->SNOW_ALB_ACCUM_B);
    fprintf(LOG_DEST, "\tSNOW_ALB_THAW_A: %.4lf\n", param->SNOW_ALB_THAW_A);
    fprintf(LOG_DEST, "\tSNOW_ALB_THAW_B: %.4lf\n", param->SNOW_ALB_THAW_B);
    fprintf(LOG_DEST, "\tSNOW_TRACESNOW: %.4lf\n", param->SNOW_TRACESNOW);
    fprintf(LOG_DEST, "\tSNOW_CONDUCT: %.4lf\n", param->SNOW_CONDUCT);
    fprintf(LOG_DEST, "\tSNOW_MAX_SNOW_TEMP: %.4lf\n",
            param->SNOW_MAX_SNOW_TEMP);
    fprintf(LOG_DEST, "\tSNOW_MIN_RAIN_TEMP: %.4lf\n",
            param->SNOW_MIN_RAIN_TEMP);
    fprintf(LOG_DEST, "\tBLOWING_KA: %.4lf\n", param->BLOWING_KA);
    fprintf(LOG_DEST, "\tBLOWING_CSALT: %.4lf\n", param->BLOWING_CSALT);
    fprintf(LOG_DEST, "\tBLOWING_UTHRESH: %.4lf\n", param->BLOWING_UTHRESH);
    fprintf(LOG_DEST, "\tBLOWING_KIN_VIS: %.4lf\n", param->BLOWING_KIN_VIS);
    fprintf(LOG_DEST, "\tBLOWING_MAX_ITER: %d\n", param->BLOWING_MAX_ITER);
    fprintf(LOG_DEST, "\tBLOWING_K: %d\n", param->BLOWING_K);
    fprintf(LOG_DEST, "\tBLOWING_SETTLING: %.4lf\n", param->BLOWING_SETTLING);
    fprintf(LOG_DEST, "\tBLOWING_NUMINCS: %d\n", param->BLOWING_NUMINCS);
    fprintf(LOG_DEST, "\tTREELINE_TEMPERATURE: %.4lf\n",
            param->TREELINE_TEMPERATURE);
    fprintf(LOG_DEST, "\tSNOW_DT: %.4lf\n", param->SNOW_DT);
    fprintf(LOG_DEST, "\tSURF_DT: %.4lf\n", param->SURF_DT);
    fprintf(LOG_DEST, "\tSOIL_DT: %.4lf\n", param->SOIL_DT);
    fprintf(LOG_DEST, "\tCANOPY_DT: %.4lf\n", param->CANOPY_DT);
    fprintf(LOG_DEST, "\tCANOPY_VP: %.4lf\n", param->CANOPY_VP);
    fprintf(LOG_DEST, "\tTOL_GRND: %.4lf\n", param->TOL_GRND);
    fprintf(LOG_DEST, "\tTOL_OVER: %.4lf\n", param->TOL_OVER);
    fprintf(LOG_DEST, "\tFROZEN_MAXITER: %d\n", param->FROZEN_MAXITER);
    fprintf(LOG_DEST, "\tNEWT_RAPH_MAXTRIAL: %d\n", param->NEWT_RAPH_MAXTRIAL);
    fprintf(LOG_DEST, "\tNEWT_RAPH_TOLX: %.4lf\n", param->NEWT_RAPH_TOLX);
    fprintf(LOG_DEST, "\tNEWT_RAPH_TOLF: %.4lf\n", param->NEWT_RAPH_TOLF);
    fprintf(LOG_DEST, "\tNEWT_RAPH_R_MAX: %.4lf\n", param->NEWT_RAPH_R_MAX);
    fprintf(LOG_DEST, "\tNEWT_RAPH_R_MIN: %.4lf\n", param->NEWT_RAPH_R_MIN);
    fprintf(LOG_DEST, "\tNEWT_RAPH_RELAX1: %.4lf\n", param->NEWT_RAPH_RELAX1);
    fprintf(LOG_DEST, "\tNEWT_RAPH_RELAX2: %.4lf\n", param->NEWT_RAPH_RELAX2);
    fprintf(LOG_DEST, "\tNEWT_RAPH_RELAX3: %.4lf\n", param->NEWT_RAPH_RELAX3);
    fprintf(LOG_DEST, "\tNEWT_RAPH_EPS2: %.4lf\n", param->NEWT_RAPH_EPS2);
    fprintf(LOG_DEST, "\tROOT_BRENT_MAXTRIES: %d\n",
            param->ROOT_BRENT_MAXTRIES);
    fprintf(LOG_DEST, "\tROOT_BRENT_MAXITER: %d\n", param->ROOT_BRENT_MAXITER);
    fprintf(LOG_DEST, "\tROOT_BRENT_TSTEP: %.4lf\n", param->ROOT_BRENT_TSTEP);
    fprintf(LOG_DEST, "\tROOT_BRENT_T: %.4lf\n", param->ROOT_BRENT_T);
    fprintf(LOG_DEST, "\tFROZEN_MAXITER: %d\n", param->FROZEN_MAXITER);
}

/******************************************************************************
 * @brief    Print save data structure.
 *****************************************************************************/
void
print_save_data(save_data_struct *save)
{
    fprintf(LOG_DEST, "save_data:\n");
    fprintf(LOG_DEST, "\ttotal_soil_moist: %.4lf\n", save->total_soil_moist);
    fprintf(LOG_DEST, "\tsurfstor: %.4lf\n", save->surfstor);
    fprintf(LOG_DEST, "\tswe: %.4lf\n", save->swe);
    fprintf(LOG_DEST, "\twdew: %.4lf\n", save->wdew);
}

/******************************************************************************
 * @brief     Print snow data structure.
 *****************************************************************************/
void
print_snow_data(snow_data_struct *snow)
{
    fprintf(LOG_DEST, "snow_data:\n");
    fprintf(LOG_DEST, "\talbedo            : %.4lf\n", snow->albedo);
    fprintf(LOG_DEST, "\tcanopy_albedo     : %.4lf\n", snow->canopy_albedo);
    fprintf(LOG_DEST, "\tcoldcontent       : %.4lf\n", snow->coldcontent);
    fprintf(LOG_DEST, "\tcoverage          : %.4lf\n", snow->coverage);
    fprintf(LOG_DEST, "\tdensity           : %.4lf\n", snow->density);
    fprintf(LOG_DEST, "\tdepth             : %.4lf\n", snow->depth);
    fprintf(LOG_DEST, "\tlast_snow         : %d\n", snow->last_snow);
    fprintf(LOG_DEST, "\tmax_snow_depth    : %.4lf\n", snow->max_snow_depth);
    fprintf(LOG_DEST, "\tMELTING           : %d\n", snow->MELTING);
    fprintf(LOG_DEST, "\tpack_temp         : %.4lf\n", snow->pack_temp);
    fprintf(LOG_DEST, "\tpack_water        : %.4lf\n", snow->pack_water);
    fprintf(LOG_DEST, "\tsnow              : %d\n", snow->snow);
    fprintf(LOG_DEST, "\tsnow_canopy       : %.4lf\n", snow->snow_canopy);
    fprintf(LOG_DEST, "\tstore_coverage    : %.4lf\n", snow->store_coverage);
    fprintf(LOG_DEST, "\tstore_snow        : %d\n", snow->store_snow);
    fprintf(LOG_DEST, "\tstore_swq         : %.4lf\n", snow->store_swq);
    fprintf(LOG_DEST, "\tsurf_temp         : %.4lf\n", snow->surf_temp);
    fprintf(LOG_DEST, "\tsurf_temp_fbcount : %u\n", snow->surf_temp_fbcount);
    fprintf(LOG_DEST, "\tsurf_temp_fbflag  : %d\n", snow->surf_temp_fbflag);
    fprintf(LOG_DEST, "\tsurf_water        : %.4lf\n", snow->surf_water);
    fprintf(LOG_DEST, "\tswq               : %.4lf\n", snow->swq);
    fprintf(LOG_DEST, "\tsnow_distrib_slope: %.4lf\n",
            snow->snow_distrib_slope);
    fprintf(LOG_DEST, "\ttmp_int_storage   : %.4lf\n", snow->tmp_int_storage);
    fprintf(LOG_DEST, "\tblowing_flux      : %.4lf\n", snow->blowing_flux);
    fprintf(LOG_DEST, "\tcanopy_vapor_flux : %.4lf\n", snow->canopy_vapor_flux);
    fprintf(LOG_DEST, "\tmass_error        : %.4lf\n", snow->mass_error);
    fprintf(LOG_DEST, "\tmelt              : %.4lf\n", snow->melt);
    fprintf(LOG_DEST, "\tQnet              : %.4lf\n", snow->Qnet);
    fprintf(LOG_DEST, "\tsurface_flux      : %.4lf\n", snow->surface_flux);
    fprintf(LOG_DEST, "\ttransport         : %.4lf\n", snow->transport);
    fprintf(LOG_DEST, "\tvapor_flux        : %.4lf\n", snow->vapor_flux);
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
    fprintf(LOG_DEST, "\tDs                    : %.4lf\n", scon->Ds);
    fprintf(LOG_DEST, "\tDsmax                 : %.4lf\n", scon->Dsmax);
    fprintf(LOG_DEST, "\tKsat                  :");
    for (i = 0; i < nlayers; i++) {
        fprintf(LOG_DEST, "\t%.4lf", scon->Ksat[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tWcr                   :");
    for (i = 0; i < nlayers; i++) {
        fprintf(LOG_DEST, "\t%.4lf", scon->Wcr[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tWpwp                  :");
    for (i = 0; i < nlayers; i++) {
        fprintf(LOG_DEST, "\t%.4lf", scon->Wpwp[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tWs                    : %.4lf\n", scon->Ws);
    fprintf(LOG_DEST, "\tAlbedoPar             : %.4f\n", scon->AlbedoPar);
    fprintf(LOG_DEST, "\talpha                 :");
    for (i = 0; i < nnodes; i++) {
        fprintf(LOG_DEST, "\t%.4lf", scon->alpha[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tannual_prec           : %.4lf\n", scon->annual_prec);
    fprintf(LOG_DEST, "\tavg_temp              : %.4lf\n", scon->avg_temp);
    fprintf(LOG_DEST, "\tavgJulyAirTemp        : %.4lf\n",
            scon->avgJulyAirTemp);
    fprintf(LOG_DEST, "\tb_infilt              : %.4lf\n", scon->b_infilt);
    fprintf(LOG_DEST, "\tbeta                  :");
    for (i = 0; i < nnodes; i++) {
        fprintf(LOG_DEST, "\t%.4lf", scon->beta[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tbubble                :");
    for (i = 0; i < nlayers; i++) {
        fprintf(LOG_DEST, "\t%.4lf", scon->bubble[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tbubble_node           :");
    for (i = 0; i < nnodes; i++) {
        fprintf(LOG_DEST, "\t%.4lf", scon->bubble_node[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tbulk_density          :");
    for (i = 0; i < nlayers; i++) {
        fprintf(LOG_DEST, "\t%.4lf", scon->bulk_density[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tbulk_dens_min         :");
    for (i = 0; i < nlayers; i++) {
        fprintf(LOG_DEST, "\t%.4lf", scon->bulk_dens_min[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tbulk_dens_org       :");
    for (i = 0; i < nlayers; i++) {
        fprintf(LOG_DEST, "\t%.4lf", scon->bulk_dens_org[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tc                     : %.4lf\n", scon->c);
    fprintf(LOG_DEST, "\tdepth                 :");
    for (i = 0; i < nlayers; i++) {
        fprintf(LOG_DEST, "\t%.4lf", scon->depth[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tdp                    : %.4lf\n", scon->dp);
    fprintf(LOG_DEST, "\tdz_node               :");
    for (i = 0; i < nnodes; i++) {
        fprintf(LOG_DEST, "\t%.4lf", scon->dz_node[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tZsum_node             :");
    for (i = 0; i < nnodes; i++) {
        fprintf(LOG_DEST, "\t%.4lf", scon->Zsum_node[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\texpt                  :");
    for (i = 0; i < nlayers; i++) {
        fprintf(LOG_DEST, "\t%.4lf", scon->expt[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\texpt_node             :");
    for (i = 0; i < nnodes; i++) {
        fprintf(LOG_DEST, "\t%.4lf", scon->expt_node[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tfrost_fract           :");
    for (i = 0; i < nfrost; i++) {
        fprintf(LOG_DEST, "\t%.4lf", scon->frost_fract[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tfrost_slope           : %.4lf\n", scon->frost_slope);
    fprintf(LOG_DEST, "\tgamma                 :");
    for (i = 0; i < nnodes; i++) {
        fprintf(LOG_DEST, "\t%.4lf", scon->gamma[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tinit_moist            :");
    for (i = 0; i < nlayers; i++) {
        fprintf(LOG_DEST, "\t%.4lf", scon->init_moist[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tmax_infil             : %.4lf\n", scon->max_infil);
    fprintf(LOG_DEST, "\tmax_moist             :");
    for (i = 0; i < nlayers; i++) {
        fprintf(LOG_DEST, "\t%.4lf", scon->max_moist[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tmax_moist_node        :");
    for (i = 0; i < nnodes; i++) {
        fprintf(LOG_DEST, "\t%.4lf", scon->max_moist_node[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tmax_snow_distrib_slope: %.4lf\n",
            scon->max_snow_distrib_slope);
    fprintf(LOG_DEST, "\tphi_s                 :");
    for (i = 0; i < nlayers; i++) {
        fprintf(LOG_DEST, "\t%.4lf", scon->phi_s[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tporosity              :");
    for (i = 0; i < nlayers; i++) {
        fprintf(LOG_DEST, "\t%.4lf", scon->porosity[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tquartz              :");
    for (i = 0; i < nlayers; i++) {
        fprintf(LOG_DEST, "\t%.4lf", scon->quartz[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\torganic               :");
    for (i = 0; i < nlayers; i++) {
        fprintf(LOG_DEST, "\t%.4lf", scon->organic[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tresid_moist           :");
    for (i = 0; i < nlayers; i++) {
        fprintf(LOG_DEST, "\t%.4lf", scon->resid_moist[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\trough                 : %.4lf\n", scon->rough);
    fprintf(LOG_DEST, "\tsnow_rough            : %.4lf\n", scon->snow_rough);
    fprintf(LOG_DEST, "\tsoil_density          :");
    for (i = 0; i < nlayers; i++) {
        fprintf(LOG_DEST, "\t%.4lf", scon->soil_density[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tsoil_dens_min         :");
    for (i = 0; i < nlayers; i++) {
        fprintf(LOG_DEST, "\t%.4lf", scon->soil_dens_min[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tsoil_dens_org         :");
    for (i = 0; i < nlayers; i++) {
        fprintf(LOG_DEST, "\t%.4lf", scon->soil_dens_org[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "BandElev                :");
    for (i = 0; i < nbands; i++) {
        fprintf(LOG_DEST, "\t%.4lf", scon->BandElev[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "AreaFract               :");
    for (i = 0; i < nbands; i++) {
        fprintf(LOG_DEST, "\t%.4lf", scon->AreaFract[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "Pfactor               :");
    for (i = 0; i < nbands; i++) {
        fprintf(LOG_DEST, "\t%.4lf", scon->Pfactor[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "Tfactor               :");
    for (i = 0; i < nbands; i++) {
        fprintf(LOG_DEST, "\t%.4lf", scon->Tfactor[i]);
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
    fprintf(LOG_DEST, "\tcell_area             : %.4lf\n", scon->cell_area);
    fprintf(LOG_DEST, "\ttime_zone_lng         : %.4f\n", scon->time_zone_lng);
    fprintf(LOG_DEST, "\tgridcel               : %d\n", scon->gridcel);
    fprintf(LOG_DEST, "\tzwtvmoist_zwt         :");
    for (i = 0; i < nlayers + 2; i++) {
        for (j = 0; j < nzwt; j++) {
            fprintf(LOG_DEST, "\t%.4lf", scon->zwtvmoist_zwt[i][j]);
        }
        fprintf(LOG_DEST, "\n\t\t\t");
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tzwtvmoist_moist       :");
    for (i = 0; i < nlayers + 2; i++) {
        for (j = 0; j < nzwt; j++) {
            fprintf(LOG_DEST, "\t%.4lf", scon->zwtvmoist_moist[i][j]);
        }
        fprintf(LOG_DEST, "\n\t\t\t");
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tslope                 : %.4lf\n", scon->slope);
    fprintf(LOG_DEST, "\taspect                : %.4lf\n", scon->aspect);
    fprintf(LOG_DEST, "\tehoriz                : %.4lf\n", scon->ehoriz);
    fprintf(LOG_DEST, "\twhoriz                : %.4lf\n", scon->whoriz);
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
    fprintf(LOG_DEST, "\tCv              : %.4lf\n", vcon->Cv);
    fprintf(LOG_DEST, "\tCv_sum          : %.4lf\n", vcon->Cv_sum);
    fprintf(LOG_DEST, "\troot            :");
    for (i = 0; i < nroots; i++) {
        fprintf(LOG_DEST, "\t%.2lf", vcon->root[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tzone_depth      :");
    for (i = 0; i < nroots; i++) {
        fprintf(LOG_DEST, "\t%.2lf", vcon->zone_depth[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tzone_fract      :");
    for (i = 0; i < nroots; i++) {
        fprintf(LOG_DEST, "\t%.2lf", vcon->zone_fract[i]);
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
            fprintf(LOG_DEST, "\t%.2lf", vcon->CanopLayerBnd[i]);
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
        fprintf(LOG_DEST, "\t%.2lf", vlib->LAI[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tWdmax         :");
    for (i = 0; i < MONTHS_PER_YEAR; i++) {
        fprintf(LOG_DEST, "\t%.2lf", vlib->Wdmax[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\talbedo        :");
    for (i = 0; i < MONTHS_PER_YEAR; i++) {
        fprintf(LOG_DEST, "\t%.2lf", vlib->albedo[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tdisplacement  :");
    for (i = 0; i < MONTHS_PER_YEAR; i++) {
        fprintf(LOG_DEST, "\t%.2lf", vlib->displacement[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\temissivity    :");
    for (i = 0; i < MONTHS_PER_YEAR; i++) {
        fprintf(LOG_DEST, "\t%.2lf", vlib->emissivity[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tNVegLibTypes  : %zu\n", vlib->NVegLibTypes);
    fprintf(LOG_DEST, "\trad_atten     : %.4lf\n", vlib->rad_atten);
    fprintf(LOG_DEST, "\trarc          : %.4lf\n", vlib->rarc);
    fprintf(LOG_DEST, "\trmin          : %.4f\n", vlib->rmin);
    fprintf(LOG_DEST, "\troughness     :");
    for (i = 0; i < MONTHS_PER_YEAR; i++) {
        fprintf(LOG_DEST, "\t%.2f", vlib->roughness[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\ttrunk_ratio   : %.4lf\n", vlib->trunk_ratio);
    fprintf(LOG_DEST, "\twind_atten    : %.4lf\n", vlib->wind_atten);
    fprintf(LOG_DEST, "\twind_h        : %.4lf\n", vlib->wind_h);
    fprintf(LOG_DEST, "\tRGL           : %.4f\n", vlib->RGL);
    fprintf(LOG_DEST, "\tveg_class     : %d\n", vlib->veg_class);
    if (carbon) {
        fprintf(LOG_DEST, "\tCtype         : %d\n", vlib->Ctype);
        fprintf(LOG_DEST, "\tMaxCarboxRate : %.4lf\n", vlib->MaxCarboxRate);
        fprintf(LOG_DEST, "\tMaxETransport : %.4lf\n", vlib->MaxETransport);
        fprintf(LOG_DEST, "\tCO2Specificity: %.4lf\n", vlib->CO2Specificity);
        fprintf(LOG_DEST, "\tLightUseEff   : %.4lf\n", vlib->LightUseEff);
        fprintf(LOG_DEST, "\tNscaleFlag    : %d\n", vlib->NscaleFlag);
        fprintf(LOG_DEST, "\tWnpp_inhib    : %.4lf\n", vlib->Wnpp_inhib);
        fprintf(LOG_DEST, "\tNPPfactor_sat : %.4lf\n", vlib->NPPfactor_sat);
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
    fprintf(LOG_DEST, "\tcanopyevap   : %.4lf\n", vvar->canopyevap);
    fprintf(LOG_DEST, "\tthroughfall  : %.4lf\n", vvar->throughfall);
    fprintf(LOG_DEST, "\tWdew         : %.4lf\n", vvar->Wdew);
    fprintf(LOG_DEST, "\tNscaleFactor :");
    for (i = 0; i < ncanopy; i++) {
        fprintf(LOG_DEST, "\t%.4lf", vvar->NscaleFactor[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\taPARLayer    :");
    for (i = 0; i < ncanopy; i++) {
        fprintf(LOG_DEST, "\t%.4lf", vvar->aPARLayer[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tCiLayer      :");
    for (i = 0; i < ncanopy; i++) {
        fprintf(LOG_DEST, "\t%.4lf", vvar->CiLayer[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\trsLayer      :");
    for (i = 0; i < ncanopy; i++) {
        fprintf(LOG_DEST, "\t%.4lf", vvar->rsLayer[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\taPAR         : %.4lf\n", vvar->aPAR);
    fprintf(LOG_DEST, "\tCi           : %.4lf\n", vvar->Ci);
    fprintf(LOG_DEST, "\trc           : %.4lf\n", vvar->rc);
    fprintf(LOG_DEST, "\tNPPfactor    : %.4lf\n", vvar->NPPfactor);
    fprintf(LOG_DEST, "\tGPP          : %.4lf\n", vvar->GPP);
    fprintf(LOG_DEST, "\tRphoto       : %.4lf\n", vvar->Rphoto);
    fprintf(LOG_DEST, "\tRdark        : %.4lf\n", vvar->Rdark);
    fprintf(LOG_DEST, "\tRmaint       : %.4lf\n", vvar->Rmaint);
    fprintf(LOG_DEST, "\tRgrowth      : %.4lf\n", vvar->Rgrowth);
    fprintf(LOG_DEST, "\tRaut         : %.4lf\n", vvar->Raut);
    fprintf(LOG_DEST, "\tNPP          : %.4lf\n", vvar->NPP);
    fprintf(LOG_DEST, "\tLitterfall   : %.4lf\n", vvar->Litterfall);
    fprintf(LOG_DEST, "\tAnnualNPP    : %.4lf\n", vvar->AnnualNPP);
    fprintf(LOG_DEST, "\tAnnualNPPPrev: %.4lf\n", vvar->AnnualNPPPrev);
}
