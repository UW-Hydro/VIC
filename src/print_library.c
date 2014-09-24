#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>

void
print_all_vars(all_vars_struct *all)
{
    printf("all_vars:\n");
    printf("\tNot implemented: Print individual structures.\n");
}

void
print_atmos_data(atmos_data_struct *atmos, size_t nr)
{
    size_t i;
    
    printf("atmos_data  :\n");
    printf("\tair_temp  :");
    for (i = 0; i <= nr; i++) {
        printf("\t%.4lf", atmos->air_temp[i]);
    }
    printf("\n");
    printf("\tCatm      :");
    for (i = 0; i <= nr; i++) {
        printf("\t%.4lf", atmos->Catm[i]);
    }
    printf("\n");
    printf("\tchannel_in:");
    for (i = 0; i <= nr; i++) {
        printf("\t%.4lf", atmos->channel_in[i]);
    }
    printf("\n");
    printf("\tcoszen    :");
    for (i = 0; i <= nr; i++) {
        printf("\t%.4lf", atmos->coszen[i]);
    }
    printf("\n");
    printf("\tdensity   :");
    for (i = 0; i <= nr; i++) {
        printf("\t%.4lf", atmos->density[i]);
    }
    printf("\n");
    printf("\tfdir      :");
    for (i = 0; i <= nr; i++) {
        printf("\t%.4lf", atmos->fdir[i]);
    }
    printf("\n");
    printf("\tlongwave  :");
    for (i = 0; i <= nr; i++) {
        printf("\t%.4lf", atmos->longwave[i]);
    }
    printf("\n");
    printf("\tout_prec  :");
    for (i = 0; i <= nr; i++) {
        printf("\t%.4lf", atmos->out_prec);
    }
    printf("\n");
    printf("\tout_rain  :");
    for (i = 0; i <= nr; i++) {
        printf("\t%.4lf", atmos->out_rain);
    }
    printf("\n");
    printf("\tout_snow  :");
    for (i = 0; i <= nr; i++) {
        printf("\t%.4lf", atmos->out_snow);
    }
    printf("\n");
    printf("\tpar       :");
    for (i = 0; i <= nr; i++) {
        printf("\t%.4lf", atmos->par[i]);
    }
    printf("\n");
    printf("\tprec      :");
    for (i = 0; i <= nr; i++) {
        printf("\t%.4lf", atmos->prec[i]);
    }
    printf("\n");
    printf("\tpressure  :");
    for (i = 0; i <= nr; i++) {
        printf("\t%.4lf", atmos->pressure[i]);
    }
    printf("\n");
    printf("\tshortwave :");
    for (i = 0; i <= nr; i++) {
        printf("\t%.4lf", atmos->shortwave[i]);
    }
    printf("\n");
    printf("\tsnowflag  :");
    for (i = 0; i <= nr; i++) {
        printf("\t%d\n", atmos->snowflag[i]);
    }
    printf("\n");
    printf("\ttskc      :");
    for (i = 0; i <= nr; i++) {
        printf("\t%.4lf", atmos->tskc[i]);
    }
    printf("\n");
    printf("\tvp        :");
    for (i = 0; i <= nr; i++) {
        printf("\t%.4lf", atmos->vp[i]);
    }
    printf("\n");
    printf("\tvpd       :");
    for (i = 0; i <= nr; i++) {
        printf("\t%.4lf", atmos->vpd[i]);
    }
    printf("\n");
    printf("\twind      :");
    for (i = 0; i <= nr; i++) {
        printf("\t%.4lf", atmos->wind[i]);
    }
    printf("\n");
}

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
    printf("\tNfrost           : %d\n", eb->Nfrost);
    printf("\tNthaw            : %d\n", eb->Nthaw);
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

void
print_Error(Error_struct *error)
{
    printf("Error:\n");
    printf("\tNot implemented: Print individual structures.\n");
}

void
print_filenames(filenames_struct *fnames)
{
    printf("filenames:\n");
    printf("\tforcing[0]   : %s\n", fnames->forcing[0]);
    printf("\tforcing[1]   : %s\n", fnames->forcing[1]);
    printf("\tf_path_pfx[0]: %s\n", fnames->f_path_pfx[0]);
    printf("\tf_path_pfx[1]: %s\n", fnames->f_path_pfx[1]);
    printf("\tglobal       : %s\n", fnames->global);
    printf("\tinit_state   : %s\n", fnames->init_state);
    printf("\tlakeparam    : %s\n", fnames->lakeparam);
    printf("\tresult_dir   : %s\n", fnames->result_dir);
    printf("\tsnowband     : %s\n", fnames->snowband);
    printf("\tsoil         : %s\n", fnames->soil);
    printf("\tstatefile    : %s\n", fnames->statefile);
    printf("\tveg          : %s\n", fnames->veg);
    printf("\tveglib       : %s\n", fnames->veglib);
}

void
print_filep(filep_struct *fp)
{
    printf("filep:\n");
    printf("\tforcing[0] : %p\n", fp->forcing[0]);
    printf("\tforcing[1] : %p\n", fp->forcing[1]);
    printf("\tglobalparam: %p\n", fp->globalparam);
    printf("\tinit_state : %p\n", fp->init_state);
    printf("\tlakeparam  : %p\n", fp->lakeparam);
    printf("\tsnowband   : %p\n", fp->snowband);
    printf("\tsoilparam  : %p\n", fp->soilparam);
    printf("\tstatefile  : %p\n", fp->statefile);
    printf("\tveglib     : %p\n", fp->veglib);
    printf("\tvegparam   : %p\n", fp->vegparam);
}

void
print_force_type(force_type_struct *force_type)
{
    printf("force_type:\n");
    printf("\tSIGNED    : %d\n", force_type->SIGNED);
    printf("\tSUPPLIED  : %d\n", force_type->SUPPLIED);
    printf("\tmultiplier: %lf\n", force_type->multiplier);
}

void
print_global_param(global_param_struct *gp)
{
    printf("global_param:\n");
    printf("\tMAX_SNOW_TEMP: %.4lf\n", gp->MAX_SNOW_TEMP);
    printf("\tMIN_RAIN_TEMP: %.4lf\n", gp->MIN_RAIN_TEMP);
    printf("\tmeasure_h    : %.4lf\n", gp->measure_h);
    printf("\twind_h       : %.4lf\n", gp->wind_h);
    printf("\tresolution   : %.4f\n", gp->resolution);
    printf("\tdt           : %d\n", gp->dt);
    printf("\tout_dt       : %d\n", gp->out_dt);
    printf("\tendday       : %d\n", gp->endday);
    printf("\tendmonth     : %d\n", gp->endmonth);
    printf("\tendyear      : %d\n", gp->endyear);
    printf("\tforceday[0]  : %d\n", gp->forceday[0]);
    printf("\tforceday[1]  : %d\n", gp->forceday[1]);
    printf("\tforcehour[0] : %d\n", gp->forcehour[0]);
    printf("\tforcehour[1] : %d\n", gp->forcehour[1]);
    printf("\tforcemonth[0]: %d\n", gp->forcemonth[0]);
    printf("\tforcemonth[1]: %d\n", gp->forcemonth[1]);
    printf("\tforceskip[0] : %d\n", gp->forceskip[0]);
    printf("\tforceskip[1] : %d\n", gp->forceskip[1]);
    printf("\tforceyear[0] : %d\n", gp->forceyear[0]);
    printf("\tforceyear[1] : %d\n", gp->forceyear[1]);
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

void
print_lake_con(lake_con_struct *lcon,
               size_t           nlnodes)
{
    size_t i;

    printf("lake_con:\n");
    printf("\tnumnod   : %d\n", lcon->numnod);
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
    printf("\teta_a    : %.4lf\n", lcon->eta_a);
    printf("\twfrac    : %.4lf\n", lcon->wfrac);
    printf("\tdepth_in : %.4lf\n", lcon->depth_in);
    printf("\tlake_idx : %d\n", lcon->lake_idx);
}

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

void
print_option(option_struct *option)
{
    printf("option:\n");
    printf("\tAboveTreelineVeg   : %d\n", option->AboveTreelineVeg);
    printf("\tAERO_RESIST_CANSNOW: %d\n", option->AERO_RESIST_CANSNOW);
    printf("\tBLOWING            : %d\n", option->BLOWING);
    printf("\tCARBON             : %d\n", option->CARBON);
    printf("\tCLOSE_ENERGY       : %d\n", option->CLOSE_ENERGY);
    printf("\tCOMPUTE_TREELINE   : %d\n", option->COMPUTE_TREELINE);
    printf("\tCONTINUEONERROR    : %d\n", option->CONTINUEONERROR);
    printf("\tCORRPREC           : %d\n", option->CORRPREC);
    printf("\tEQUAL_AREA         : %d\n", option->EQUAL_AREA);
    printf("\tEXP_TRANS          : %d\n", option->EXP_TRANS);
    printf("\tFROZEN_SOIL        : %d\n", option->FROZEN_SOIL);
    printf("\tFULL_ENERGY        : %d\n", option->FULL_ENERGY);
    printf("\tGRND_FLUX_TYPE     : %d\n", option->GRND_FLUX_TYPE);
    printf("\tIMPLICIT           : %d\n", option->IMPLICIT);
    printf("\tJULY_TAVG_SUPPLIED : %d\n", option->JULY_TAVG_SUPPLIED);
    printf("\tLAKES              : %d\n", option->LAKES);
    printf("\tLW_CLOUD           : %d\n", option->LW_CLOUD);
    printf("\tLW_TYPE            : %d\n", option->LW_TYPE);
    printf("\tMIN_WIND_SPEED     : %.4f\n", option->MIN_WIND_SPEED);
    printf("\tMTCLIM_SWE_CORR    : %d\n", option->MTCLIM_SWE_CORR);
    printf("\tNcanopy            : %d\n", option->Ncanopy);
    printf("\tNfrost             : %d\n", option->Nfrost);
    printf("\tNlakenode          : %d\n", option->Nlakenode);
    printf("\tNlayer             : %d\n", option->Nlayer);
    printf("\tNnode              : %d\n", option->Nnode);
    printf("\tNOFLUX             : %d\n", option->NOFLUX);
    // printf("\tNVEGTYPES          : %d\n", option->NVEGTYPES);
    printf("\tPLAPSE             : %d\n", option->PLAPSE);
    printf("\tRC_MODE            : %d\n", option->RC_MODE);
    printf("\tROOT_ZONES         : %d\n", option->ROOT_ZONES);
    printf("\tQUICK_FLUX         : %d\n", option->QUICK_FLUX);
    printf("\tQUICK_SOLVE        : %d\n", option->QUICK_SOLVE);
    printf("\tSHARE_LAYER_MOIST  : %d\n", option->SHARE_LAYER_MOIST);
    printf("\tSNOW_DENSITY       : %d\n", option->SNOW_DENSITY);
    printf("\tSNOW_BAND          : %d\n", option->SNOW_BAND);
    printf("\tSNOW_STEP          : %d\n", option->SNOW_STEP);
    printf("\tSPATIAL_FROST      : %d\n", option->SPATIAL_FROST);
    printf("\tSPATIAL_SNOW       : %d\n", option->SPATIAL_SNOW);
    printf("\tSW_PREC_THRESH     : %.4f\n", option->SW_PREC_THRESH);
    printf("\tTFALLBACK          : %d\n", option->TFALLBACK);
    printf("\tVP_INTERP          : %d\n", option->VP_INTERP);
    printf("\tVP_ITER            : %d\n", option->VP_ITER);
    printf("\tALMA_INPUT         : %d\n", option->ALMA_INPUT);
    printf("\tBASEFLOW           : %d\n", option->BASEFLOW);
    printf("\tGRID_DECIMAL       : %d\n", option->GRID_DECIMAL);
    printf("\tVEGLIB_PHOTO       : %d\n", option->VEGLIB_PHOTO);
    printf("\tVEGPARAM_LAI       : %d\n", option->VEGPARAM_LAI);
    printf("\tLAI_SRC            : %d\n", option->LAI_SRC);
    printf("\tLAKE_PROFILE       : %d\n", option->LAKE_PROFILE);
    printf("\tORGANIC_FRACT      : %d\n", option->ORGANIC_FRACT);
    printf("\tBINARY_STATE_FILE  : %d\n", option->BINARY_STATE_FILE);
    printf("\tINIT_STATE         : %d\n", option->INIT_STATE);
    printf("\tSAVE_STATE         : %d\n", option->SAVE_STATE);
    printf("\tALMA_OUTPUT        : %d\n", option->ALMA_OUTPUT);
    printf("\tBINARY_OUTPUT      : %d\n", option->BINARY_OUTPUT);
    printf("\tCOMPRESS           : %d\n", option->COMPRESS);
    printf("\tMOISTFRACT         : %d\n", option->MOISTFRACT);
    printf("\tNoutfiles          : %d\n", option->Noutfiles);
    printf("\tOUTPUT_FORCE       : %d\n", option->OUTPUT_FORCE);
    printf("\tPRT_HEADER         : %d\n", option->PRT_HEADER);
    printf("\tPRT_SNOW_BAND      : %d\n", option->PRT_SNOW_BAND);
}

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

void
print_out_data_file(out_data_file_struct *outf)
{
    size_t i;

    printf("\tprefix: %s\n", outf->prefix);
    printf("\tfilename: %s\n", outf->filename);
    printf("\tfh: %p\n", outf->fh);
    printf("\tnvars: %d\n", outf->nvars);
    printf("\tvarid: %p\n", outf->varid);
}

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

void
print_save_data(save_data_struct *save)
{
    printf("save_data:\n");
    printf("\ttotal_soil_moist: %.4lf\n", save->total_soil_moist);
    printf("\tsurfstor: %.4lf\n", save->surfstor);
    printf("\tswe: %.4lf\n", save->swe);
    printf("\twdew: %.4lf\n", save->wdew);
}

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
    printf("\tsurf_temp_fbcount : %.4lf\n", snow->surf_temp_fbcount);
    printf("\tsurf_temp_fbflag  : %.4lf\n", snow->surf_temp_fbflag);
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
    printf("\tvegetat_type_num: %d\n", vcon->vegetat_type_num);
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

void
print_veg_lib(veg_lib_struct *vlib,
              char            carbon)
{
    size_t i;

    printf("veg_lib:\n");
    printf("\toverstory     : %d\n", vlib->overstory);
    printf("\tLAI           :");
    for (i = 0; i < MONTHSPERYEAR; i++) {
        printf("\t%.2lf", vlib->LAI[i]);
    }
    printf("\n");
    printf("\tWdmax         :");
    for (i = 0; i < MONTHSPERYEAR; i++) {
        printf("\t%.2lf", vlib->Wdmax[i]);
    }
    printf("\n");
    printf("\talbedo        :");
    for (i = 0; i < MONTHSPERYEAR; i++) {
        printf("\t%.2lf", vlib->albedo[i]);
    }
    printf("\n");
    printf("\tdisplacement  :");
    for (i = 0; i < MONTHSPERYEAR; i++) {
        printf("\t%.2lf", vlib->displacement[i]);
    }
    printf("\n");
    printf("\temissivity    :");
    for (i = 0; i < MONTHSPERYEAR; i++) {
        printf("\t%.2lf", vlib->emissivity[i]);
    }
    printf("\n");
    printf("\tNVegLibTypes  : %d\n", vlib->NVegLibTypes);
    printf("\trad_atten     : %.4lf\n", vlib->rad_atten);
    printf("\trarc          : %.4lf\n", vlib->rarc);
    printf("\trmin          : %.4f\n", vlib->rmin);
    printf("\troughness     :");
    for (i = 0; i < MONTHSPERYEAR; i++) {
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
