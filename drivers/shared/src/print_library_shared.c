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

    printf("cell_data:\n");
    printf("\taero_resist :");
    for (i = 0; i < 2; i++) {
        printf("\t%.4lf", cell->aero_resist[i]);
    }
    printf("\n");
    printf("\tasat        : %.4lf\n", cell->asat);
    printf("\tbaseflow    : %.4lf\n", cell->baseflow);
    printf("\tCLitter     : %.4lf\n", cell->CLitter);
    printf("\tCInter      : %.4lf\n", cell->CInter);
    printf("\tCSlow       : %.4lf\n", cell->CSlow);
    printf("\tinflow      : %.4lf\n", cell->inflow);
    printf("\tpot_evap    :");
    for (i = 0; i < npet; i++) {
        printf("\t%.4lf", cell->pot_evap[i]);
    }
    printf("\n");
    printf("\trunoff      : %.4lf\n", cell->runoff);
    for (i = 0; i < nlayers; i++) {
        printf("\tlayer %zd   :\n", i);
        print_layer_data(&(cell->layer[i]), nfrost);
    }
    printf("\tRhLitter    : %.4lf\n", cell->RhLitter);
    printf("\tRhLitter2Atm: %.4lf\n", cell->RhLitter2Atm);
    printf("\tRhInter     : %.4lf\n", cell->RhInter);
    printf("\tRhSlow      : %.4lf\n", cell->RhSlow);
    printf("\tRhTot       : %.4lf\n", cell->RhTot);
    printf("\trootmoist   : %.4lf\n", cell->rootmoist);
    printf("\twetness     : %.4lf\n", cell->wetness);
    printf("\tzwt         : %.4lf\n", cell->zwt);
    printf("\tzwt_lumped  : %.4lf\n", cell->zwt_lumped);
}

/******************************************************************************
 * @brief    Print day-month-year structure.
 *****************************************************************************/
void
print_dmy(dmy_struct *dmy)
{
    printf("dmy:\n");
    printf("\tday        : %d\n", dmy->day);
    printf("\tday_in_year: %d\n", dmy->day_in_year);
    printf("\thour       : %d\n", dmy->hour);
    printf("\tmonth      : %d\n", dmy->month);
    printf("\tyear       : %d\n", dmy->year);
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

    printf("energy_bal:\n");
    printf("\tAlbedoLake       : %.4lf\n", eb->AlbedoLake);
    printf("\tAlbedoOver       : %.4lf\n", eb->AlbedoOver);
    printf("\tAlbedoUnder      : %.4lf\n", eb->AlbedoUnder);
    printf("\tCs               :");
    for (i = 0; i < 2; i++) {
        printf("\t%.4lf", eb->Cs[i]);
    }
    printf("\n");
    printf("\tCs_node          :");
    for (i = 0; i < nnodes; i++) {
        printf("\t%.4lf", eb->Cs_node[i]);
    }
    printf("\n");
    printf("\tfdepth           :");
    for (i = 0; i < nfronts; i++) {
        printf("\t%.4lf", eb->fdepth[i]);
    }
    printf("\n");
    printf("\tfrozen           : %d\n", eb->frozen);
    printf("\tice              :");
    for (i = 0; i < nnodes; i++) {
        printf("\t%.4lf", eb->ice[i]);
    }
    printf("\n");
    printf("\tkappa            :");
    for (i = 0; i < 2; i++) {
        printf("\t%.4lf", eb->kappa[i]);
    }
    printf("\n");
    printf("\tkappa_node       :");
    for (i = 0; i < nnodes; i++) {
        printf("\t%.4lf", eb->kappa_node[i]);
    }
    printf("\n");
    printf("\tmoist            :");
    for (i = 0; i < nnodes; i++) {
        printf("\t%.4lf", eb->moist[i]);
    }
    printf("\n");
    printf("\tNfrost           : %zu\n", eb->Nfrost);
    printf("\tNthaw            : %zu\n", eb->Nthaw);
    printf("\tT                :");
    for (i = 0; i < nnodes; i++) {
        printf("\t%.4lf", eb->T[i]);
    }
    printf("\n");
    printf("\tT_fbflag         :");
    for (i = 0; i < nnodes; i++) {
        printf("\t%d", eb->T_fbflag[i]);
    }
    printf("\n");
    printf("\tT_fbcount        :");
    for (i = 0; i < nnodes; i++) {
        printf("\t%d", eb->T_fbcount[i]);
    }
    printf("\n");
    printf("\tT1_index         : %d\n", eb->T1_index);
    printf("\tTcanopy          : %.4lf\n", eb->Tcanopy);
    printf("\tTcanopy_fbflag   : %d\n", eb->Tcanopy_fbflag);
    printf("\tTcanopy_fbcount  : %d\n", eb->Tcanopy_fbcount);
    printf("\ttdepth           :");
    for (i = 0; i < nfronts; i++) {
        printf("\t%.4lf", eb->tdepth[i]);
    }
    printf("\n");
    printf("\tTfoliage         : %.4lf\n", eb->Tfoliage);
    printf("\tTfoliage_fbflag  : %d\n", eb->Tfoliage_fbflag);
    printf("\tTfoliage_fbcount : %d\n", eb->Tfoliage_fbcount);
    printf("\tTsurf            : %.4lf\n", eb->Tsurf);
    printf("\tTsurf_fbflag     : %d\n", eb->Tsurf_fbflag);
    printf("\tTsurf_fbcount    : %d\n", eb->Tsurf_fbcount);
    printf("\tunfrozen         : %.4lf\n", eb->unfrozen);
    printf("\tadvected_sensible: %.4lf\n", eb->advected_sensible);
    printf("\tadvection        : %.4lf\n", eb->advection);
    printf("\tAtmosError       : %.4lf\n", eb->AtmosError);
    printf("\tAtmosLatent      : %.4lf\n", eb->AtmosLatent);
    printf("\tAtmosLatentSub   : %.4lf\n", eb->AtmosLatentSub);
    printf("\tAtmosSensible    : %.4lf\n", eb->AtmosSensible);
    printf("\tcanopy_advection : %.4lf\n", eb->canopy_advection);
    printf("\tcanopy_latent    : %.4lf\n", eb->canopy_latent);
    printf("\tcanopy_latent_sub: %.4lf\n", eb->canopy_latent_sub);
    printf("\tcanopy_refreeze  : %.4lf\n", eb->canopy_refreeze);
    printf("\tcanopy_sensible  : %.4lf\n", eb->canopy_sensible);
    printf("\tdeltaCC          : %.4lf\n", eb->deltaCC);
    printf("\tdeltaH           : %.4lf\n", eb->deltaH);
    printf("\terror            : %.4lf\n", eb->error);
    printf("\tfusion           : %.4lf\n", eb->fusion);
    printf("\tgrnd_flux        : %.4lf\n", eb->grnd_flux);
    printf("\tlatent           : %.4lf\n", eb->latent);
    printf("\tlatent_sub       : %.4lf\n", eb->latent_sub);
    printf("\tlongwave         : %.4lf\n", eb->longwave);
    printf("\tLongOverIn       : %.4lf\n", eb->LongOverIn);
    printf("\tLongUnderIn      : %.4lf\n", eb->LongUnderIn);
    printf("\tLongUnderOut     : %.4lf\n", eb->LongUnderOut);
    printf("\tmelt_energy      : %.4lf\n", eb->melt_energy);
    printf("\tNetLongAtmos     : %.4lf\n", eb->NetLongAtmos);
    printf("\tNetLongOver      : %.4lf\n", eb->NetLongOver);
    printf("\tNetLongUnder     : %.4lf\n", eb->NetLongUnder);
    printf("\tNetShortAtmos    : %.4lf\n", eb->NetShortAtmos);
    printf("\tNetShortGrnd     : %.4lf\n", eb->NetShortGrnd);
    printf("\tNetShortOver     : %.4lf\n", eb->NetShortOver);
    printf("\tNetShortUnder    : %.4lf\n", eb->NetShortUnder);
    printf("\tout_long_canopy  : %.4lf\n", eb->out_long_canopy);
    printf("\tout_long_surface : %.4lf\n", eb->out_long_surface);
    printf("\trefreeze_energy  : %.4lf\n", eb->refreeze_energy);
    printf("\tsensible         : %.4lf\n", eb->sensible);
    printf("\tshortwave        : %.4lf\n", eb->shortwave);
    printf("\tShortOverIn      : %.4lf\n", eb->ShortOverIn);
    printf("\tShortUnderIn     : %.4lf\n", eb->ShortUnderIn);
    printf("\tsnow_flux        : %.4lf\n", eb->snow_flux);
}

/******************************************************************************
 * @brief    Print filenames structure.
 *****************************************************************************/
void
print_filenames(filenames_struct *fnames)
{
    printf("filenames:\n");
    printf("\tforcing[0]   : %s\n", fnames->forcing[0]);
    printf("\tforcing[1]   : %s\n", fnames->forcing[1]);
    printf("\tf_path_pfx[0]: %s\n", fnames->f_path_pfx[0]);
    printf("\tf_path_pfx[1]: %s\n", fnames->f_path_pfx[1]);
    printf("\tglobal       : %s\n", fnames->global);
    printf("\tconstants    : %s\n", fnames->constants);
    printf("\tdomain       : %s\n", fnames->domain);
    printf("\tinit_state   : %s\n", fnames->init_state);
    printf("\tlakeparam    : %s\n", fnames->lakeparam);
    printf("\tresult_dir   : %s\n", fnames->result_dir);
    printf("\tsnowband     : %s\n", fnames->snowband);
    printf("\tsoil         : %s\n", fnames->soil);
    printf("\tstatefile    : %s\n", fnames->statefile);
    printf("\tveg          : %s\n", fnames->veg);
    printf("\tveglib       : %s\n", fnames->veglib);
    printf("\tlog_path     : %s\n", fnames->log_path);
}

/******************************************************************************
 * @brief    Print file path structure.
 *****************************************************************************/
void
print_filep(filep_struct *fp)
{
    printf("filep:\n");
    printf("\tforcing[0] : %p\n", fp->forcing[0]);
    printf("\tforcing[1] : %p\n", fp->forcing[1]);
    printf("\tglobalparam: %p\n", fp->globalparam);
    printf("\tconstants  : %p\n", fp->constants);
    printf("\tdomain     : %p\n", fp->domain);
    printf("\tinit_state : %p\n", fp->init_state);
    printf("\tlakeparam  : %p\n", fp->lakeparam);
    printf("\tsnowband   : %p\n", fp->snowband);
    printf("\tsoilparam  : %p\n", fp->soilparam);
    printf("\tstatefile  : %p\n", fp->statefile);
    printf("\tveglib     : %p\n", fp->veglib);
    printf("\tvegparam   : %p\n", fp->vegparam);
    printf("\tlogfile    : %p\n", fp->logfile);
}

/******************************************************************************
 * @brief    Print forcing type structure.
 *****************************************************************************/
void
print_force_type(force_type_struct *force_type)
{
    printf("force_type:\n");
    printf("\tSIGNED    : %d\n", force_type->SIGNED);
    printf("\tSUPPLIED  : %d\n", force_type->SUPPLIED);
    printf("\tmultiplier: %lf\n", force_type->multiplier);
}

/******************************************************************************
 * @brief    Print global parameters structure.
 *****************************************************************************/
void
print_global_param(global_param_struct *gp)
{
    size_t i;

    printf("global_param:\n");
    printf("\tmeasure_h    : %.4lf\n", gp->measure_h);
    printf("\twind_h       : %.4lf\n", gp->wind_h);
    printf("\tresolution   : %.4f\n", gp->resolution);
    printf("\tdt           : %d\n", gp->dt);
    printf("\tout_dt       : %d\n", gp->out_dt);
    printf("\tendday       : %d\n", gp->endday);
    printf("\tendmonth     : %d\n", gp->endmonth);
    printf("\tendyear      : %d\n", gp->endyear);
    for (i = 0; i < 2; i++) {
        printf("\tforceday[%zd]   : %d\n", i, gp->forceday[i]);
        printf("\tforcehour[%zd]  : %d\n", i, gp->forcehour[i]);
        printf("\tforcemonth[%zd] : %d\n", i, gp->forcemonth[i]);
        printf("\tforceoffset[%zd]: %d\n", i, gp->forceoffset[i]);
        printf("\tforceskip[%zd]  : %d\n", i, gp->forceskip[i]);
        printf("\tforceyear[%zd]  : %d\n", i, gp->forceyear[i]);
    }
    printf("\tnrecs        : %d\n", gp->nrecs);
    printf("\tskipyear     : %d\n", gp->skipyear);
    printf("\tstartday     : %d\n", gp->startday);
    printf("\tstarthour    : %d\n", gp->starthour);
    printf("\tstartmonth   : %d\n", gp->startmonth);
    printf("\tstartyear    : %d\n", gp->startyear);
    printf("\tstateday     : %d\n", gp->stateday);
    printf("\tstatemonth   : %d\n", gp->statemonth);
    printf("\tstateyear    : %d\n", gp->stateyear);
}

/******************************************************************************
 * @brief    Print lake_con_structure.
 *****************************************************************************/
void
print_lake_con(lake_con_struct *lcon,
               size_t           nlnodes)
{
    size_t i;

    printf("lake_con:\n");
    printf("\tnumnod   : %zu\n", lcon->numnod);
    printf("\tz        :");
    for (i = 0; i < nlnodes; i++) {
        printf("\t%.4lf", lcon->z[i]);
    }
    printf("\n");
    printf("\tbasin    :");
    for (i = 0; i < nlnodes; i++) {
        printf("\t%.4lf", lcon->basin[i]);
    }
    printf("\n");
    printf("\tCl       :");
    for (i = 0; i < nlnodes; i++) {
        printf("\t%.4lf", lcon->Cl[i]);
    }
    printf("\n");
    printf("\tb        : %.4lf\n", lcon->b);
    printf("\tmaxdepth : %.4lf\n", lcon->maxdepth);
    printf("\tmindepth : %.4lf\n", lcon->mindepth);
    printf("\tmaxvolume: %.4lf\n", lcon->maxvolume);
    printf("\tminvolume: %.4lf\n", lcon->minvolume);
    printf("\tbpercent : %.4f\n", lcon->bpercent);
    printf("\trpercent : %.4f\n", lcon->rpercent);
    printf("\twfrac    : %.4lf\n", lcon->wfrac);
    printf("\tdepth_in : %.4lf\n", lcon->depth_in);
    printf("\tlake_idx : %d\n", lcon->lake_idx);
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

    printf("lake_var:\n");
    printf("\tactivenod      : %d\n", lvar->activenod);
    printf("\tdz             : %.4lf\n", lvar->dz);
    printf("\tsurfdz         : %.4lf\n", lvar->surfdz);
    printf("\tldepth         : %.4lf\n", lvar->ldepth);
    printf("\tsurface        :");
    for (i = 0; i < nlnodes + 1; i++) {
        printf("\t%.4lf", lvar->surface[i]);
    }
    printf("\n");
    printf("\tsarea          : %.4lf\n", lvar->sarea);
    printf("\tsarea_save     : %.4lf\n", lvar->sarea_save);
    printf("\tvolume         : %.4lf\n", lvar->volume);
    printf("\tvolume_save    : %.4lf\n", lvar->volume_save);
    printf("\ttemp           :");
    for (i = 0; i < nlnodes; i++) {
        printf("\t%.4lf", lvar->temp[i]);
    }
    printf("\n");
    printf("\ttempavg        : %.4lf\n", lvar->tempavg);
    printf("\tareai          : %.4lf\n", lvar->areai);
    printf("\tnew_ice_area   : %.4lf\n", lvar->new_ice_area);
    printf("\tice_water_eq   : %.4lf\n", lvar->ice_water_eq);
    printf("\thice           : %.4lf\n", lvar->hice);
    printf("\ttempi          : %.4lf\n", lvar->tempi);
    printf("\tswe            : %.4lf\n", lvar->swe);
    printf("\tswe_save       : %.4lf\n", lvar->swe_save);
    printf("\tsurf_temp      : %.4lf\n", lvar->surf_temp);
    printf("\tpack_temp      : %.4lf\n", lvar->pack_temp);
    printf("\tcoldcontent    : %.4lf\n", lvar->coldcontent);
    printf("\tsurf_water     : %.4lf\n", lvar->surf_water);
    printf("\tpack_water     : %.4lf\n", lvar->pack_water);
    printf("\tSAlbedo        : %.4lf\n", lvar->SAlbedo);
    printf("\tsdepth         : %.4lf\n", lvar->sdepth);
    printf("\taero_resist    : %.4lf\n", lvar->aero_resist);
    printf("\tdensity        :");
    for (i = 0; i < nlnodes; i++) {
        printf("\t%.4lf", lvar->density[i]);
    }
    printf("\n");
    printf("\tbaseflow_in    : %.4lf\n", lvar->baseflow_in);
    printf("\tbaseflow_out   : %.4lf\n", lvar->baseflow_out);
    printf("\tchannel_in     : %.4lf\n", lvar->channel_in);
    printf("\tevapw          : %.4lf\n", lvar->evapw);
    printf("\tice_throughfall: %.4lf\n", lvar->ice_throughfall);
    printf("\tprec           : %.4lf\n", lvar->prec);
    printf("\trecharge       : %.4lf\n", lvar->recharge);
    printf("\trunoff_in      : %.4lf\n", lvar->runoff_in);
    printf("\trunoff_out     : %.4lf\n", lvar->runoff_out);
    printf("\tsnowmlt        : %.4lf\n", lvar->snowmlt);
    printf("\tvapor_flux     : %.4lf\n", lvar->vapor_flux);
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

    printf("layer_data:\n");
    printf("\tCs   : %.4lf\n", ldata->Cs);
    printf("\tT    : %.4lf\n", ldata->T);
    printf("\tevap : %.4lf\n", ldata->evap);
    printf("\tice  :");
    for (i = 0; i < nfrost; i++) {
        printf("\t%.4lf", ldata->ice[i]);
    }
    printf("\n");
    printf("\tkappa: %.4lf\n", ldata->kappa);
    printf("\tmoist: %.4lf\n", ldata->moist);
    printf("\tphi  : %.4lf\n", ldata->phi);
    printf("\tzwt  : %.4lf\n", ldata->zwt);
}

/******************************************************************************
 * @brief    Print options structure.
 *****************************************************************************/
void
print_option(option_struct *option)
{
    printf("option:\n");
    printf("\tAboveTreelineVeg     : %d\n", option->AboveTreelineVeg);
    printf("\tAERO_RESIST_CANSNOW  : %d\n", option->AERO_RESIST_CANSNOW);
    printf("\tBLOWING              : %d\n", option->BLOWING);
    printf("\tBLOWING_VAR_THRESHOLD: %d\n", option->BLOWING_VAR_THRESHOLD);
    printf("\tBLOWING_CALC_PROB    : %d\n", option->BLOWING_CALC_PROB);
    printf("\tBLOWING_SIMPLE       : %d\n", option->BLOWING_SIMPLE);
    printf("\tBLOWING_FETCH        : %d\n", option->BLOWING_FETCH);
    printf("\tBLOWING_SPATIAL_WIND : %d\n", option->BLOWING_SPATIAL_WIND);
    printf("\tCARBON               : %d\n", option->CARBON);
    printf("\tCLOSE_ENERGY         : %d\n", option->CLOSE_ENERGY);
    printf("\tCOMPUTE_TREELINE     : %d\n", option->COMPUTE_TREELINE);
    printf("\tCONTINUEONERROR      : %d\n", option->CONTINUEONERROR);
    printf("\tCORRPREC             : %d\n", option->CORRPREC);
    printf("\tEQUAL_AREA           : %d\n", option->EQUAL_AREA);
    printf("\tEXP_TRANS            : %d\n", option->EXP_TRANS);
    printf("\tFROZEN_SOIL          : %d\n", option->FROZEN_SOIL);
    printf("\tFULL_ENERGY          : %d\n", option->FULL_ENERGY);
    printf("\tGRND_FLUX_TYPE       : %d\n", option->GRND_FLUX_TYPE);
    printf("\tIMPLICIT             : %d\n", option->IMPLICIT);
    printf("\tJULY_TAVG_SUPPLIED   : %d\n", option->JULY_TAVG_SUPPLIED);
    printf("\tLAKES                : %d\n", option->LAKES);
    printf("\tLW_CLOUD             : %d\n", option->LW_CLOUD);
    printf("\tLW_TYPE              : %d\n", option->LW_TYPE);
    printf("\tMTCLIM_SWE_CORR      : %d\n", option->MTCLIM_SWE_CORR);
    printf("\tNcanopy              : %zu\n", option->Ncanopy);
    printf("\tNfrost               : %zu\n", option->Nfrost);
    printf("\tNlakenode            : %zu\n", option->Nlakenode);
    printf("\tNlayer               : %zu\n", option->Nlayer);
    printf("\tNnode                : %zu\n", option->Nnode);
    printf("\tNOFLUX               : %d\n", option->NOFLUX);
    printf("\tNVEGTYPES            : %zu\n", option->NVEGTYPES);
    printf("\tPLAPSE               : %d\n", option->PLAPSE);
    printf("\tRC_MODE              : %d\n", option->RC_MODE);
    printf("\tROOT_ZONES           : %zu\n", option->ROOT_ZONES);
    printf("\tQUICK_FLUX           : %d\n", option->QUICK_FLUX);
    printf("\tQUICK_SOLVE          : %d\n", option->QUICK_SOLVE);
    printf("\tSHARE_LAYER_MOIST    : %d\n", option->SHARE_LAYER_MOIST);
    printf("\tSNOW_DENSITY         : %d\n", option->SNOW_DENSITY);
    printf("\tSNOW_BAND            : %zu\n", option->SNOW_BAND);
    printf("\tSNOW_STEP            : %d\n", option->SNOW_STEP);
    printf("\tSPATIAL_FROST        : %d\n", option->SPATIAL_FROST);
    printf("\tSPATIAL_SNOW         : %d\n", option->SPATIAL_SNOW);
    printf("\tTFALLBACK            : %d\n", option->TFALLBACK);
    printf("\tVP_INTERP            : %d\n", option->VP_INTERP);
    printf("\tVP_ITER              : %d\n", option->VP_ITER);
    printf("\tALMA_INPUT           : %d\n", option->ALMA_INPUT);
    printf("\tBASEFLOW             : %d\n", option->BASEFLOW);
    printf("\tGRID_DECIMAL         : %d\n", option->GRID_DECIMAL);
    printf("\tVEGLIB_PHOTO         : %d\n", option->VEGLIB_PHOTO);
    printf("\tVEGPARAM_LAI         : %d\n", option->VEGPARAM_LAI);
    printf("\tLAI_SRC              : %d\n", option->LAI_SRC);
    printf("\tLAKE_PROFILE         : %d\n", option->LAKE_PROFILE);
    printf("\tORGANIC_FRACT        : %d\n", option->ORGANIC_FRACT);
    printf("\tBINARY_STATE_FILE    : %d\n", option->BINARY_STATE_FILE);
    printf("\tINIT_STATE           : %d\n", option->INIT_STATE);
    printf("\tSAVE_STATE           : %d\n", option->SAVE_STATE);
    printf("\tALMA_OUTPUT          : %d\n", option->ALMA_OUTPUT);
    printf("\tBINARY_OUTPUT        : %d\n", option->BINARY_OUTPUT);
    printf("\tCOMPRESS             : %d\n", option->COMPRESS);
    printf("\tMOISTFRACT           : %d\n", option->MOISTFRACT);
    printf("\tNoutfiles            : %zu\n", option->Noutfiles);
    printf("\tOUTPUT_FORCE         : %d\n", option->OUTPUT_FORCE);
    printf("\tPRT_HEADER           : %d\n", option->PRT_HEADER);
    printf("\tPRT_SNOW_BAND        : %d\n", option->PRT_SNOW_BAND);
}

/******************************************************************************
 * @brief    Print out data structure.
 *****************************************************************************/
void
print_out_data(out_data_struct *out,
               size_t           nelem)
{
    size_t i;

    printf("out_data:\n");
    printf("\tvarname: %s\n", out->varname);
    printf("\twrite: %d\n", out->write);
    printf("\tformat: %s\n", out->format);
    printf("\ttype: %d\n", out->type);
    printf("\tmult: %.4f\n", out->mult);
    printf("\tnelem: %d\n", out->nelem);
    printf("\tdata:");
    for (i = 0; i < nelem; i++) {
        printf("\t%.4lf", out->data[i]);
    }
    printf("\n");
    printf("\taggdata:");
    for (i = 0; i < nelem; i++) {
        printf("\t%.4lf", out->aggdata[i]);
    }
    printf("\n");
}

/******************************************************************************
 * @brief    Print out data file structure.
 *****************************************************************************/
void
print_out_data_file(out_data_file_struct *outf)
{
    printf("\tprefix: %s\n", outf->prefix);
    printf("\tfilename: %s\n", outf->filename);
    printf("\tfh: %p\n", outf->fh);
    printf("\tnvars: %zu\n", outf->nvars);
    printf("\tvarid: %p\n", outf->varid);
}

/******************************************************************************
 * @brief    print param set structure.
 *****************************************************************************/
void
print_param_set(param_set_struct *param_set)
{
    size_t i;

    printf("param_set:\n");
    for (i = 0; i < N_FORCING_TYPES; i++) {
        print_force_type(&(param_set->TYPE[i]));
    }
    printf("\tFORCE_DT    : %d %d\n", param_set->FORCE_DT[0],
           param_set->FORCE_DT[1]);
    printf("\tFORCE_ENDIAN: %d %d\n", param_set->FORCE_ENDIAN[0],
           param_set->FORCE_ENDIAN[1]);
    printf("\tFORCE_FORMAT: %d %d\n", param_set->FORCE_FORMAT[0],
           param_set->FORCE_FORMAT[1]);
    printf("\tFORCE_INDEX :\n");
    for (i = 0; i < N_FORCING_TYPES; i++) {
        printf("\t\t%zd: %d %d\n", i, param_set->FORCE_INDEX[0][i],
               param_set->FORCE_INDEX[1][i]);
    }
    printf("\tN_TYPES     : %d %d\n", param_set->N_TYPES[0],
           param_set->N_TYPES[1]);
}

/******************************************************************************
 * @brief    Print model parameters.
 *****************************************************************************/
void
print_parameters(parameters_struct *param)
{
    printf("parameters:\n");
    printf("\tLAPSE_RATE: %.4lf\n", param->LAPSE_RATE);
    printf("\tGAUGE_HEIGHT: %.4lf\n", param->GAUGE_HEIGHT);
    printf("\tWIND_SPEED_DEFAULT: %.4lf\n", param->WIND_SPEED_DEFAULT);
    printf("\tWIND_SPEED_MIN: %.4lf\n", param->WIND_SPEED_MIN);
    printf("\tHUGE_RESIST: %.4lf\n", param->HUGE_RESIST);
    printf("\tALBEDO_BARE_SOIL: %.4lf\n", param->ALBEDO_BARE_SOIL);
    printf("\tALBEDO_H20_SURF: %.4lf\n", param->ALBEDO_H20_SURF);
    printf("\tEMISS_GRND: %.4lf\n", param->EMISS_GRND);
    printf("\tEMISS_VEG: %.4lf\n", param->EMISS_VEG);
    printf("\tEMISS_ICE: %.4lf\n", param->EMISS_ICE);
    printf("\tEMISS_SNOW: %.4lf\n", param->EMISS_SNOW);
    printf("\tEMISS_H2O: %.4lf\n", param->EMISS_H2O);
    printf("\tSOIL_RESID_MOIST: %.4lf\n", param->SOIL_RESID_MOIST);
    printf("\tSOIL_SLAB_MOIST_FRACT: %.4lf\n", param->SOIL_SLAB_MOIST_FRACT);
    printf("\tVEG_LAI_SNOW_MULTIPLIER: %.4lf\n",
           param->VEG_LAI_SNOW_MULTIPLIER);
    printf("\tVEG_MIN_INTERCEPTION_STORAGE: %.4lf\n",
           param->VEG_MIN_INTERCEPTION_STORAGE);
    printf("\tVEG_LAI_WATER_FACTOR: %.4lf\n", param->VEG_LAI_WATER_FACTOR);
    printf("\tCANOPY_CLOSURE: %.4lf\n", param->CANOPY_CLOSURE);
    printf("\tCANOPY_RSMAX: %.4lf\n", param->CANOPY_RSMAX);
    printf("\tCANOPY_VPDMINFACTOR: %.4lf\n", param->CANOPY_VPDMINFACTOR);
    printf("\tMTCLIM_SOLAR_CONSTANT: %.4lf\n", param->MTCLIM_SOLAR_CONSTANT);
    printf("\tMTCLIM_TDAYCOEF: %.4lf\n", param->MTCLIM_TDAYCOEF);
    printf("\tMTCLIM_SNOW_TCRIT: %.4lf\n", param->MTCLIM_SNOW_TCRIT);
    printf("\tMTCLIM_SNOW_TRATE: %.4lf\n", param->MTCLIM_SNOW_TRATE);
    printf("\tMTCLIM_TBASE: %.4lf\n", param->MTCLIM_TBASE);
    printf("\tMTCLIM_ABASE: %.4lf\n", param->MTCLIM_ABASE);
    printf("\tMTCLIM_C: %.4lf\n", param->MTCLIM_C);
    printf("\tMTCLIM_B0: %.4lf\n", param->MTCLIM_B0);
    printf("\tMTCLIM_B1: %.4lf\n", param->MTCLIM_B1);
    printf("\tMTCLIM_B2: %.4lf\n", param->MTCLIM_B2);
    printf("\tMTCLIM_RAIN_SCALAR: %.4lf\n", param->MTCLIM_RAIN_SCALAR);
    printf("\tMTCLIM_DIF_ALB: %.4lf\n", param->MTCLIM_DIF_ALB);
    printf("\tMTCLIM_SC_INT: %.4lf\n", param->MTCLIM_SC_INT);
    printf("\tMTCLIM_SC_SLOPE: %.4lf\n", param->MTCLIM_SC_SLOPE);
    printf("\tMTCLIM_SRADDT: %.4lf\n", param->MTCLIM_SRADDT);
    printf("\tMTCLIM_SW_PREC_THRESH: %.4lf\n", param->MTCLIM_SW_PREC_THRESH);
    printf("\tLAKE_TMELT: %.4lf\n", param->LAKE_TMELT);
    printf("\tLAKE_MAX_SURFACE: %.4lf\n", param->LAKE_MAX_SURFACE);
    printf("\tLAKE_BETA: %.4lf\n", param->LAKE_BETA);
    printf("\tLAKE_FRACMIN: %.4lf\n", param->LAKE_FRACMIN);
    printf("\tLAKE_FRACLIM: %.4lf\n", param->LAKE_FRACLIM);
    printf("\tLAKE_DM: %.4lf\n", param->LAKE_DM);
    printf("\tLAKE_SNOWCRIT: %.4lf\n", param->LAKE_SNOWCRIT);
    printf("\tLAKE_ZWATER: %.4lf\n", param->LAKE_ZWATER);
    printf("\tLAKE_ZSNOW: %.4lf\n", param->LAKE_ZSNOW);
    printf("\tLAKE_RHOSNOW: %.4lf\n", param->LAKE_RHOSNOW);
    printf("\tLAKE_CONDI: %.4lf\n", param->LAKE_CONDI);
    printf("\tLAKE_CONDS: %.4lf\n", param->LAKE_CONDS);
    printf("\tLAKE_LAMISW: %.4lf\n", param->LAKE_LAMISW);
    printf("\tLAKE_LAMILW: %.4lf\n", param->LAKE_LAMILW);
    printf("\tLAKE_LAMSSW: %.4lf\n", param->LAKE_LAMSSW);
    printf("\tLAKE_LAMSLW: %.4lf\n", param->LAKE_LAMSLW);
    printf("\tLAKE_LAMWSW: %.4lf\n", param->LAKE_LAMWSW);
    printf("\tLAKE_LAMWLW: %.4lf\n", param->LAKE_LAMWLW);
    printf("\tLAKE_A1: %.4lf\n", param->LAKE_A1);
    printf("\tLAKE_A2: %.4lf\n", param->LAKE_A2);
    printf("\tLAKE_QWTAU: %.4lf\n", param->LAKE_QWTAU);
    printf("\tLAKE_MAX_ITER: %d\n", param->LAKE_MAX_ITER);
    printf("\tSVP_A: %.4lf\n", param->SVP_A);
    printf("\tSVP_B: %.4lf\n", param->SVP_B);
    printf("\tSVP_C: %.4lf\n", param->SVP_C);
    printf("\tCARBON_CATMCURRENT: %.4lf\n", param->CARBON_CATMCURRENT);
    printf("\tCARBON_SW2PAR: %.4lf\n", param->CARBON_SW2PAR);
    printf("\tPHOTO_OMEGA: %.4lf\n", param->PHOTO_OMEGA);
    printf("\tPHOTO_LAIMAX: %.4lf\n", param->PHOTO_LAIMAX);
    printf("\tPHOTO_LAILIMIT: %.4lf\n", param->PHOTO_LAILIMIT);
    printf("\tPHOTO_LAIMIN: %.4lf\n", param->PHOTO_LAIMIN);
    printf("\tPHOTO_EPAR: %.4lf\n", param->PHOTO_EPAR);
    printf("\tPHOTO_FCMAX: %.4lf\n", param->PHOTO_FCMAX);
    printf("\tPHOTO_FCMIN: %.4lf\n", param->PHOTO_FCMIN);
    printf("\tPHOTO_ZENITHMIN: %.4lf\n", param->PHOTO_ZENITHMIN);
    printf("\tPHOTO_ZENITHMINPAR: %.4lf\n", param->PHOTO_ZENITHMINPAR);
    printf("\tPHOTO_ALBSOIPARMIN: %.4lf\n", param->PHOTO_ALBSOIPARMIN);
    printf("\tPHOTO_MINMAXETRANS: %.4lf\n", param->PHOTO_MINMAXETRANS);
    printf("\tPHOTO_MINSTOMCOND: %.4lf\n", param->PHOTO_MINSTOMCOND);
    printf("\tPHOTO_FCI1C3: %.4lf\n", param->PHOTO_FCI1C3);
    printf("\tPHOTO_FCI1C4: %.4lf\n", param->PHOTO_FCI1C4);
    printf("\tPHOTO_OX: %.4lf\n", param->PHOTO_OX);
    printf("\tPHOTO_KC: %.4lf\n", param->PHOTO_KC);
    printf("\tPHOTO_KO: %.4lf\n", param->PHOTO_KO);
    printf("\tPHOTO_EC: %.4lf\n", param->PHOTO_EC);
    printf("\tPHOTO_EO: %.4lf\n", param->PHOTO_EO);
    printf("\tPHOTO_EV: %.4lf\n", param->PHOTO_EV);
    printf("\tPHOTO_ER: %.4lf\n", param->PHOTO_ER);
    printf("\tPHOTO_ALC3: %.4lf\n", param->PHOTO_ALC3);
    printf("\tPHOTO_FRDC3: %.4lf\n", param->PHOTO_FRDC3);
    printf("\tPHOTO_EK: %.4lf\n", param->PHOTO_EK);
    printf("\tPHOTO_ALC4: %.4lf\n", param->PHOTO_ALC4);
    printf("\tPHOTO_FRDC4: %.4lf\n", param->PHOTO_FRDC4);
    printf("\tPHOTO_THETA: %.4lf\n", param->PHOTO_THETA);
    printf("\tPHOTO_FRLEAF: %.4lf\n", param->PHOTO_FRLEAF);
    printf("\tPHOTO_FRGROWTH: %.4lf\n", param->PHOTO_FRGROWTH);
    printf("\tSRESP_E0_LT: %.4lf\n", param->SRESP_E0_LT);
    printf("\tSRESP_T0_LT: %.4lf\n", param->SRESP_T0_LT);
    printf("\tSRESP_WMINFM: %.4lf\n", param->SRESP_WMINFM);
    printf("\tSRESP_WMAXFM: %.4lf\n", param->SRESP_WMAXFM);
    printf("\tSRESP_WOPTFM: %.4lf\n", param->SRESP_WOPTFM);
    printf("\tSRESP_RHSAT: %.4lf\n", param->SRESP_RHSAT);
    printf("\tSRESP_RFACTOR: %.4lf\n", param->SRESP_RFACTOR);
    printf("\tSRESP_TAULITTER: %.4lf\n", param->SRESP_TAULITTER);
    printf("\tSRESP_TAUINTER: %.4lf\n", param->SRESP_TAUINTER);
    printf("\tSRESP_TAUSLOW: %.4lf\n", param->SRESP_TAUSLOW);
    printf("\tSRESP_FAIR: %.4lf\n", param->SRESP_FAIR);
    printf("\tSRESP_FINTER: %.4lf\n", param->SRESP_FINTER);
    printf("\tSNOW_MAX_SURFACE_SWE: %.4lf\n", param->SNOW_MAX_SURFACE_SWE);
    printf("\tSNOW_LIQUID_WATER_CAPACITY: %.4lf\n",
           param->SNOW_LIQUID_WATER_CAPACITY);
    printf("\tSNOW_NEW_SNOW_DENSITY: %.4lf\n", param->SNOW_NEW_SNOW_DENSITY);
    printf("\tSNOW_DENS_DMLIMIT: %.4lf\n", param->SNOW_DENS_DMLIMIT);
    printf("\tSNOW_DENS_MAX_CHANGE: %.4lf\n", param->SNOW_DENS_MAX_CHANGE);
    printf("\tSNOW_DENS_ETA0: %.4lf\n", param->SNOW_DENS_ETA0);
    printf("\tSNOW_DENS_C1: %.4lf\n", param->SNOW_DENS_C1);
    printf("\tSNOW_DENS_C2: %.4lf\n", param->SNOW_DENS_C2);
    printf("\tSNOW_DENS_C5: %.4lf\n", param->SNOW_DENS_C5);
    printf("\tSNOW_DENS_C6: %.4lf\n", param->SNOW_DENS_C6);
    printf("\tSNOW_DENS_F: %.4lf\n", param->SNOW_DENS_F);
    printf("\tSNOW_MIN_SWQ_EB_THRES: %.4lf\n", param->SNOW_MIN_SWQ_EB_THRES);
    printf("\tSNOW_A1: %.4lf\n", param->SNOW_A1);
    printf("\tSNOW_A2: %.4lf\n", param->SNOW_A2);
    printf("\tSNOW_L1: %.4lf\n", param->SNOW_L1);
    printf("\tSNOW_L2: %.4lf\n", param->SNOW_L2);
    printf("\tSNOW_NEW_SNOW_ALB: %.4lf\n", param->SNOW_NEW_SNOW_ALB);
    printf("\tSNOW_ALB_ACCUM_A: %.4lf\n", param->SNOW_ALB_ACCUM_A);
    printf("\tSNOW_ALB_ACCUM_B: %.4lf\n", param->SNOW_ALB_ACCUM_B);
    printf("\tSNOW_ALB_THAW_A: %.4lf\n", param->SNOW_ALB_THAW_A);
    printf("\tSNOW_ALB_THAW_B: %.4lf\n", param->SNOW_ALB_THAW_B);
    printf("\tSNOW_TRACESNOW: %.4lf\n", param->SNOW_TRACESNOW);
    printf("\tSNOW_CONDUCT: %.4lf\n", param->SNOW_CONDUCT);
    printf("\tSNOW_MAX_SNOW_TEMP: %.4lf\n", param->SNOW_MAX_SNOW_TEMP);
    printf("\tSNOW_MIN_RAIN_TEMP: %.4lf\n", param->SNOW_MIN_RAIN_TEMP);
    printf("\tBLOWING_KA: %.4lf\n", param->BLOWING_KA);
    printf("\tBLOWING_CSALT: %.4lf\n", param->BLOWING_CSALT);
    printf("\tBLOWING_UTHRESH: %.4lf\n", param->BLOWING_UTHRESH);
    printf("\tBLOWING_KIN_VIS: %.4lf\n", param->BLOWING_KIN_VIS);
    printf("\tBLOWING_MAX_ITER: %d\n", param->BLOWING_MAX_ITER);
    printf("\tBLOWING_K: %d\n", param->BLOWING_K);
    printf("\tBLOWING_SETTLING: %.4lf\n", param->BLOWING_SETTLING);
    printf("\tBLOWING_NUMINCS: %d\n", param->BLOWING_NUMINCS);
    printf("\tTREELINE_TEMPERATURE: %.4lf\n", param->TREELINE_TEMPERATURE);
    printf("\tSNOW_DT: %.4lf\n", param->SNOW_DT);
    printf("\tSURF_DT: %.4lf\n", param->SURF_DT);
    printf("\tSOIL_DT: %.4lf\n", param->SOIL_DT);
    printf("\tCANOPY_DT: %.4lf\n", param->CANOPY_DT);
    printf("\tCANOPY_VP: %.4lf\n", param->CANOPY_VP);
    printf("\tTOL_GRND: %.4lf\n", param->TOL_GRND);
    printf("\tTOL_OVER: %.4lf\n", param->TOL_OVER);
    printf("\tFROZEN_MAXITER: %d\n", param->FROZEN_MAXITER);
    printf("\tNEWT_RAPH_MAXTRIAL: %d\n", param->NEWT_RAPH_MAXTRIAL);
    printf("\tNEWT_RAPH_TOLX: %.4lf\n", param->NEWT_RAPH_TOLX);
    printf("\tNEWT_RAPH_TOLF: %.4lf\n", param->NEWT_RAPH_TOLF);
    printf("\tNEWT_RAPH_R_MAX: %.4lf\n", param->NEWT_RAPH_R_MAX);
    printf("\tNEWT_RAPH_R_MIN: %.4lf\n", param->NEWT_RAPH_R_MIN);
    printf("\tNEWT_RAPH_RELAX1: %.4lf\n", param->NEWT_RAPH_RELAX1);
    printf("\tNEWT_RAPH_RELAX2: %.4lf\n", param->NEWT_RAPH_RELAX2);
    printf("\tNEWT_RAPH_RELAX3: %.4lf\n", param->NEWT_RAPH_RELAX3);
    printf("\tNEWT_RAPH_EPS2: %.4lf\n", param->NEWT_RAPH_EPS2);
    printf("\tROOT_BRENT_MAXTRIES: %d\n", param->ROOT_BRENT_MAXTRIES);
    printf("\tROOT_BRENT_MAXITER: %d\n", param->ROOT_BRENT_MAXITER);
    printf("\tROOT_BRENT_TSTEP: %.4lf\n", param->ROOT_BRENT_TSTEP);
    printf("\tROOT_BRENT_T: %.4lf\n", param->ROOT_BRENT_T);
    printf("\tFROZEN_MAXITER: %d\n", param->FROZEN_MAXITER);
}

/******************************************************************************
 * @brief    Print save data structure.
 *****************************************************************************/
void
print_save_data(save_data_struct *save)
{
    printf("save_data:\n");
    printf("\ttotal_soil_moist: %.4lf\n", save->total_soil_moist);
    printf("\tsurfstor: %.4lf\n", save->surfstor);
    printf("\tswe: %.4lf\n", save->swe);
    printf("\twdew: %.4lf\n", save->wdew);
}

/******************************************************************************
 * @brief     Print snow data structure.
 *****************************************************************************/
void
print_snow_data(snow_data_struct *snow)
{
    printf("snow_data:\n");
    printf("\talbedo            : %.4lf\n", snow->albedo);
    printf("\tcanopy_albedo     : %.4lf\n", snow->canopy_albedo);
    printf("\tcoldcontent       : %.4lf\n", snow->coldcontent);
    printf("\tcoverage          : %.4lf\n", snow->coverage);
    printf("\tdensity           : %.4lf\n", snow->density);
    printf("\tdepth             : %.4lf\n", snow->depth);
    printf("\tlast_snow         : %d\n", snow->last_snow);
    printf("\tmax_snow_depth    : %.4lf\n", snow->max_snow_depth);
    printf("\tMELTING           : %d\n", snow->MELTING);
    printf("\tpack_temp         : %.4lf\n", snow->pack_temp);
    printf("\tpack_water        : %.4lf\n", snow->pack_water);
    printf("\tsnow              : %d\n", snow->snow);
    printf("\tsnow_canopy       : %.4lf\n", snow->snow_canopy);
    printf("\tstore_coverage    : %.4lf\n", snow->store_coverage);
    printf("\tstore_snow        : %d\n", snow->store_snow);
    printf("\tstore_swq         : %.4lf\n", snow->store_swq);
    printf("\tsurf_temp         : %.4lf\n", snow->surf_temp);
    printf("\tsurf_temp_fbcount : %u\n", snow->surf_temp_fbcount);
    printf("\tsurf_temp_fbflag  : %d\n", snow->surf_temp_fbflag);
    printf("\tsurf_water        : %.4lf\n", snow->surf_water);
    printf("\tswq               : %.4lf\n", snow->swq);
    printf("\tsnow_distrib_slope: %.4lf\n", snow->snow_distrib_slope);
    printf("\ttmp_int_storage   : %.4lf\n", snow->tmp_int_storage);
    printf("\tblowing_flux      : %.4lf\n", snow->blowing_flux);
    printf("\tcanopy_vapor_flux : %.4lf\n", snow->canopy_vapor_flux);
    printf("\tmass_error        : %.4lf\n", snow->mass_error);
    printf("\tmelt              : %.4lf\n", snow->melt);
    printf("\tQnet              : %.4lf\n", snow->Qnet);
    printf("\tsurface_flux      : %.4lf\n", snow->surface_flux);
    printf("\ttransport         : %.4lf\n", snow->transport);
    printf("\tvapor_flux        : %.4lf\n", snow->vapor_flux);
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

    printf("soil_con:\n");
    printf("\tFS_ACTIVE             : %d\n", scon->FS_ACTIVE);
    printf("\tDs                    : %.4lf\n", scon->Ds);
    printf("\tDsmax                 : %.4lf\n", scon->Dsmax);
    printf("\tKsat                  :");
    for (i = 0; i < nlayers; i++) {
        printf("\t%.4lf", scon->Ksat[i]);
    }
    printf("\n");
    printf("\tWcr                   :");
    for (i = 0; i < nlayers; i++) {
        printf("\t%.4lf", scon->Wcr[i]);
    }
    printf("\n");
    printf("\tWpwp                  :");
    for (i = 0; i < nlayers; i++) {
        printf("\t%.4lf", scon->Wpwp[i]);
    }
    printf("\n");
    printf("\tWs                    : %.4lf\n", scon->Ws);
    printf("\tAlbedoPar             : %.4f\n", scon->AlbedoPar);
    printf("\talpha                 :");
    for (i = 0; i < nnodes; i++) {
        printf("\t%.4lf", scon->alpha[i]);
    }
    printf("\n");
    printf("\tannual_prec           : %.4lf\n", scon->annual_prec);
    printf("\tavg_temp              : %.4lf\n", scon->avg_temp);
    printf("\tavgJulyAirTemp        : %.4lf\n", scon->avgJulyAirTemp);
    printf("\tb_infilt              : %.4lf\n", scon->b_infilt);
    printf("\tbeta                  :");
    for (i = 0; i < nnodes; i++) {
        printf("\t%.4lf", scon->beta[i]);
    }
    printf("\n");
    printf("\tbubble                :");
    for (i = 0; i < nlayers; i++) {
        printf("\t%.4lf", scon->bubble[i]);
    }
    printf("\n");
    printf("\tbubble_node           :");
    for (i = 0; i < nnodes; i++) {
        printf("\t%.4lf", scon->bubble_node[i]);
    }
    printf("\n");
    printf("\tbulk_density          :");
    for (i = 0; i < nlayers; i++) {
        printf("\t%.4lf", scon->bulk_density[i]);
    }
    printf("\n");
    printf("\tbulk_dens_min         :");
    for (i = 0; i < nlayers; i++) {
        printf("\t%.4lf", scon->bulk_dens_min[i]);
    }
    printf("\n");
    printf("\tbulk_dens_org       :");
    for (i = 0; i < nlayers; i++) {
        printf("\t%.4lf", scon->bulk_dens_org[i]);
    }
    printf("\n");
    printf("\tc                     : %.4lf\n", scon->c);
    printf("\tdepth                 :");
    for (i = 0; i < nlayers; i++) {
        printf("\t%.4lf", scon->depth[i]);
    }
    printf("\n");
    printf("\tdp                    : %.4lf\n", scon->dp);
    printf("\tdz_node               :");
    for (i = 0; i < nnodes; i++) {
        printf("\t%.4lf", scon->dz_node[i]);
    }
    printf("\n");
    printf("\tZsum_node             :");
    for (i = 0; i < nnodes; i++) {
        printf("\t%.4lf", scon->Zsum_node[i]);
    }
    printf("\n");
    printf("\texpt                  :");
    for (i = 0; i < nlayers; i++) {
        printf("\t%.4lf", scon->expt[i]);
    }
    printf("\n");
    printf("\texpt_node             :");
    for (i = 0; i < nnodes; i++) {
        printf("\t%.4lf", scon->expt_node[i]);
    }
    printf("\n");
    printf("\tfrost_fract           :");
    for (i = 0; i < nfrost; i++) {
        printf("\t%.4lf", scon->frost_fract[i]);
    }
    printf("\n");
    printf("\tfrost_slope           : %.4lf\n", scon->frost_slope);
    printf("\tgamma                 :");
    for (i = 0; i < nnodes; i++) {
        printf("\t%.4lf", scon->gamma[i]);
    }
    printf("\n");
    printf("\tinit_moist            :");
    for (i = 0; i < nlayers; i++) {
        printf("\t%.4lf", scon->init_moist[i]);
    }
    printf("\n");
    printf("\tmax_infil             : %.4lf\n", scon->max_infil);
    printf("\tmax_moist             :");
    for (i = 0; i < nlayers; i++) {
        printf("\t%.4lf", scon->max_moist[i]);
    }
    printf("\n");
    printf("\tmax_moist_node        :");
    for (i = 0; i < nnodes; i++) {
        printf("\t%.4lf", scon->max_moist_node[i]);
    }
    printf("\n");
    printf("\tmax_snow_distrib_slope: %.4lf\n", scon->max_snow_distrib_slope);
    printf("\tphi_s                 :");
    for (i = 0; i < nlayers; i++) {
        printf("\t%.4lf", scon->phi_s[i]);
    }
    printf("\n");
    printf("\tporosity              :");
    for (i = 0; i < nlayers; i++) {
        printf("\t%.4lf", scon->porosity[i]);
    }
    printf("\n");
    printf("\tquartz              :");
    for (i = 0; i < nlayers; i++) {
        printf("\t%.4lf", scon->quartz[i]);
    }
    printf("\n");
    printf("\torganic               :");
    for (i = 0; i < nlayers; i++) {
        printf("\t%.4lf", scon->organic[i]);
    }
    printf("\n");
    printf("\tresid_moist           :");
    for (i = 0; i < nlayers; i++) {
        printf("\t%.4lf", scon->resid_moist[i]);
    }
    printf("\n");
    printf("\trough                 : %.4lf\n", scon->rough);
    printf("\tsnow_rough            : %.4lf\n", scon->snow_rough);
    printf("\tsoil_density          :");
    for (i = 0; i < nlayers; i++) {
        printf("\t%.4lf", scon->soil_density[i]);
    }
    printf("\n");
    printf("\tsoil_dens_min         :");
    for (i = 0; i < nlayers; i++) {
        printf("\t%.4lf", scon->soil_dens_min[i]);
    }
    printf("\n");
    printf("\tsoil_dens_org         :");
    for (i = 0; i < nlayers; i++) {
        printf("\t%.4lf", scon->soil_dens_org[i]);
    }
    printf("\n");
    printf("BandElev                :");
    for (i = 0; i < nbands; i++) {
        printf("\t%.4lf", scon->BandElev[i]);
    }
    printf("\n");
    printf("AreaFract               :");
    for (i = 0; i < nbands; i++) {
        printf("\t%.4lf", scon->AreaFract[i]);
    }
    printf("\n");
    printf("Pfactor               :");
    for (i = 0; i < nbands; i++) {
        printf("\t%.4lf", scon->Pfactor[i]);
    }
    printf("\n");
    printf("Tfactor               :");
    for (i = 0; i < nbands; i++) {
        printf("\t%.4lf", scon->Tfactor[i]);
    }
    printf("\n");
    printf("AboveTreeLine         :");
    for (i = 0; i < nbands; i++) {
        printf("\t%d", scon->AboveTreeLine[i]);
    }
    printf("\n");
    printf("\televation             : %.4f\n", scon->elevation);
    printf("\tlat                   : %.4f\n", scon->lat);
    printf("\tlng                   : %.4f\n", scon->lng);
    printf("\tcell_area             : %.4lf\n", scon->cell_area);
    printf("\ttime_zone_lng         : %.4f\n", scon->time_zone_lng);
    printf("\tgridcel               : %d\n", scon->gridcel);
    printf("\tzwtvmoist_zwt         :");
    for (i = 0; i < nlayers + 2; i++) {
        for (j = 0; j < nzwt; j++) {
            printf("\t%.4lf", scon->zwtvmoist_zwt[i][j]);
        }
        printf("\n\t\t\t");
    }
    printf("\n");
    printf("\tzwtvmoist_moist       :");
    for (i = 0; i < nlayers + 2; i++) {
        for (j = 0; j < nzwt; j++) {
            printf("\t%.4lf", scon->zwtvmoist_moist[i][j]);
        }
        printf("\n\t\t\t");
    }
    printf("\n");
    printf("\tslope                 : %.4lf\n", scon->slope);
    printf("\taspect                : %.4lf\n", scon->aspect);
    printf("\tehoriz                : %.4lf\n", scon->ehoriz);
    printf("\twhoriz                : %.4lf\n", scon->whoriz);
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

    printf("veg_con:\n");
    printf("\tCv              : %.4lf\n", vcon->Cv);
    printf("\tCv_sum          : %.4lf\n", vcon->Cv_sum);
    printf("\troot            :");
    for (i = 0; i < nroots; i++) {
        printf("\t%.2lf", vcon->root[i]);
    }
    printf("\n");
    printf("\tzone_depth      :");
    for (i = 0; i < nroots; i++) {
        printf("\t%.2lf", vcon->zone_depth[i]);
    }
    printf("\n");
    printf("\tzone_fract      :");
    for (i = 0; i < nroots; i++) {
        printf("\t%.2lf", vcon->zone_fract[i]);
    }
    printf("\n");
    printf("\tveg_class       : %d\n", vcon->veg_class);
    printf("\tvegetat_type_num: %zu\n", vcon->vegetat_type_num);
    if (blowing) {
        printf("\tsigma_slope     : %.4f\n", vcon->sigma_slope);
        printf("\tlag_one         : %.4f\n", vcon->lag_one);
        printf("\tfetch           : %.4f\n", vcon->fetch);
    }
    if (lake) {
        printf("\tLAKE            : %d\n", vcon->LAKE);
    }
    if (carbon) {
        printf("\tCanopLayerBnd   :");
        for (i = 0; i < ncanopy; i++) {
            printf("\t%.2lf", vcon->CanopLayerBnd[i]);
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

    printf("veg_lib:\n");
    printf("\toverstory     : %d\n", vlib->overstory);
    printf("\tLAI           :");
    for (i = 0; i < MONTHS_PER_YEAR; i++) {
        printf("\t%.2lf", vlib->LAI[i]);
    }
    printf("\n");
    printf("\tWdmax         :");
    for (i = 0; i < MONTHS_PER_YEAR; i++) {
        printf("\t%.2lf", vlib->Wdmax[i]);
    }
    printf("\n");
    printf("\talbedo        :");
    for (i = 0; i < MONTHS_PER_YEAR; i++) {
        printf("\t%.2lf", vlib->albedo[i]);
    }
    printf("\n");
    printf("\tdisplacement  :");
    for (i = 0; i < MONTHS_PER_YEAR; i++) {
        printf("\t%.2lf", vlib->displacement[i]);
    }
    printf("\n");
    printf("\temissivity    :");
    for (i = 0; i < MONTHS_PER_YEAR; i++) {
        printf("\t%.2lf", vlib->emissivity[i]);
    }
    printf("\n");
    printf("\tNVegLibTypes  : %zu\n", vlib->NVegLibTypes);
    printf("\trad_atten     : %.4lf\n", vlib->rad_atten);
    printf("\trarc          : %.4lf\n", vlib->rarc);
    printf("\trmin          : %.4f\n", vlib->rmin);
    printf("\troughness     :");
    for (i = 0; i < MONTHS_PER_YEAR; i++) {
        printf("\t%.2f", vlib->roughness[i]);
    }
    printf("\n");
    printf("\ttrunk_ratio   : %.4lf\n", vlib->trunk_ratio);
    printf("\twind_atten    : %.4lf\n", vlib->wind_atten);
    printf("\twind_h        : %.4lf\n", vlib->wind_h);
    printf("\tRGL           : %.4f\n", vlib->RGL);
    printf("\tveg_class     : %d\n", vlib->veg_class);
    if (carbon) {
        printf("\tCtype         : %d\n", vlib->Ctype);
        printf("\tMaxCarboxRate : %.4lf\n", vlib->MaxCarboxRate);
        printf("\tMaxETransport : %.4lf\n", vlib->MaxETransport);
        printf("\tCO2Specificity: %.4lf\n", vlib->CO2Specificity);
        printf("\tLightUseEff   : %.4lf\n", vlib->LightUseEff);
        printf("\tNscaleFlag    : %d\n", vlib->NscaleFlag);
        printf("\tWnpp_inhib    : %.4lf\n", vlib->Wnpp_inhib);
        printf("\tNPPfactor_sat : %.4lf\n", vlib->NPPfactor_sat);
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

    printf("veg_var:\n");
    printf("\tcanopyevap   : %.4lf\n", vvar->canopyevap);
    printf("\tthroughfall  : %.4lf\n", vvar->throughfall);
    printf("\tWdew         : %.4lf\n", vvar->Wdew);
    printf("\tNscaleFactor :");
    for (i = 0; i < ncanopy; i++) {
        printf("\t%.4lf", vvar->NscaleFactor[i]);
    }
    printf("\n");
    printf("\taPARLayer    :");
    for (i = 0; i < ncanopy; i++) {
        printf("\t%.4lf", vvar->aPARLayer[i]);
    }
    printf("\n");
    printf("\tCiLayer      :");
    for (i = 0; i < ncanopy; i++) {
        printf("\t%.4lf", vvar->CiLayer[i]);
    }
    printf("\n");
    printf("\trsLayer      :");
    for (i = 0; i < ncanopy; i++) {
        printf("\t%.4lf", vvar->rsLayer[i]);
    }
    printf("\n");
    printf("\taPAR         : %.4lf\n", vvar->aPAR);
    printf("\tCi           : %.4lf\n", vvar->Ci);
    printf("\trc           : %.4lf\n", vvar->rc);
    printf("\tNPPfactor    : %.4lf\n", vvar->NPPfactor);
    printf("\tGPP          : %.4lf\n", vvar->GPP);
    printf("\tRphoto       : %.4lf\n", vvar->Rphoto);
    printf("\tRdark        : %.4lf\n", vvar->Rdark);
    printf("\tRmaint       : %.4lf\n", vvar->Rmaint);
    printf("\tRgrowth      : %.4lf\n", vvar->Rgrowth);
    printf("\tRaut         : %.4lf\n", vvar->Raut);
    printf("\tNPP          : %.4lf\n", vvar->NPP);
    printf("\tLitterfall   : %.4lf\n", vvar->Litterfall);
    printf("\tAnnualNPP    : %.4lf\n", vvar->AnnualNPP);
    printf("\tAnnualNPPPrev: %.4lf\n", vvar->AnnualNPPPrev);
}
