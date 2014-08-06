#include <vic_def.h>
#include <vic_run.h>
#include <vic_driver_image.h>

void
print_domain(domain_struct *domain,
             bool           print_loc)
{
    int i;

    printf("domain:\n");
    printf("\tncells_global: %zd\n", domain->ncells_global);
    printf("\tn_nx: %zd\n", domain->n_nx);
    printf("\tn_ny: %zd\n", domain->n_ny);
    printf("\tncells_local: %zd\n", domain->ncells_local);
    printf("\tlocations: %p\n", domain->locations);
    if (print_loc) {
        for (i = 0; i < domain->ncells_global; i++) {
            print_location(&(domain->locations[i]));
        }
    }
}

void
print_filep(filep_struct *fp)
{
    printf("filep:\n");
    printf("\tforcing[0]: %p\n", fp->forcing[0]);
    printf("\tforcing[1]: %p\n", fp->forcing[1]);
    printf("\tglobalparam: %p\n", fp->globalparam);
    printf("\tdomain: %p\n", fp->domain);
    printf("\tinit_state: %p\n", fp->init_state);
    printf("\tlakeparam: %p\n", fp->lakeparam);
    printf("\tsnowband: %p\n", fp->snowband);
    printf("\tsoilparam: %p\n", fp->soilparam);
    printf("\tstatefile: %p\n", fp->statefile);
    printf("\tveglib: %p\n", fp->veglib);
    printf("\tvegparam: %p\n", fp->vegparam);
}

void
print_force_type(force_type_struct *force_type)
{
    printf("force_type:\n");
    printf("\tSIGNED: %d, SUPPLIED: %d, multiplier: %lf\n",
           force_type->SIGNED, force_type->SUPPLIED,
           force_type->multiplier);
}

void
print_location(location_struct *location)
{
    printf("%zd: (%zd, %zd) {%zd: (%zd, %zd)}\n",
           location->global_cell_idx, location->global_x_idx,
           location->global_y_idx, location->local_cell_idx,
           location->local_x_idx, location->local_y_idx);
    printf("\t(%lf, %lf): %lf %lf\n",
           location->longitude, location->latitude,
           location->area, location->frac);
}

void
print_param_set(param_set_struct *param_set)
{
    size_t i;

    printf("param_set:\n");
    for (i = 0; i < N_FORCING_TYPES; i++) {
        print_force_type(&(param_set->TYPE[i]));
    }
    printf("\tFORCE_DT: %d %d\n", param_set->FORCE_DT[0],
           param_set->FORCE_DT[1]);
    printf("\tFORCE_ENDIAN: %d %d\n", param_set->FORCE_ENDIAN[0],
           param_set->FORCE_ENDIAN[1]);
    printf("\tFORCE_FORMAT: %d %d\n", param_set->FORCE_FORMAT[0],
           param_set->FORCE_FORMAT[1]);
    printf("\tFORCE_INDEX:\n");
    for (i = 0; i < N_FORCING_TYPES; i++) {
        printf("\t\t%zd: %d %d\n", i, param_set->FORCE_INDEX[0][i],
               param_set->FORCE_INDEX[1][i]);
    }
    printf("\tN_TYPES: %d %d\n", param_set->N_TYPES[0], param_set->N_TYPES[1]);
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
    printf("\tCv: %.4lf\n", vcon->Cv);
    printf("\tCv_sum: %.4lf\n", vcon->Cv_sum);
    printf("\troot:");
    for (i = 0; i < nroots; i++) {
        printf("\t%.2lf", vcon->root[i]);
    }
    printf("\n");
    printf("\tzone_depth:");
    for (i = 0; i < nroots; i++) {
        printf("\t%.2lf", vcon->zone_depth[i]);
    }
    printf("\n");
    printf("\tzone_fract:");
    for (i = 0; i < nroots; i++) {
        printf("\t%.2lf", vcon->zone_fract[i]);
    }
    printf("\n");
    printf("\tveg_class: %d\n", vcon->veg_class);
    printf("\tvegetat_type_num: %d\n", vcon->vegetat_type_num);
    if (blowing) {
        printf("\tsigma_slope: %.4f\n", vcon->sigma_slope);
        printf("\tlag_one: %.4f\n", vcon->lag_one);
        printf("\tfetch: %.4f\n", vcon->fetch);
    }
    if (lake) {
        printf("\tLAKE: %d\n", vcon->LAKE);
    }
    if (carbon) {
        printf("\tCanopLayerBnd:");
        for (i = 0; i < nroots; i++) {
            printf("\t%.2lf", vcon->CanopLayerBnd[i]);
        }
    }
}

void
print_veg_con_map(veg_con_map_struct *veg_con_map)
{
    size_t i;

    printf("veg_con_map:\n");
    printf("\tnv_types: %zd\n", veg_con_map->nv_types);
    printf("\tnv_active: %zd\n", veg_con_map->nv_active);
    for (i = 0; i < veg_con_map->nv_types; i++) {
        printf("\t%zd: %d (vidx) %lf (Cv)\n", i, veg_con_map->vidx[i],
               veg_con_map->Cv[i]);
    }
}

void
print_veg_lib(veg_lib_struct *vlib,
              char            carbon)
{
    size_t i;

    printf("veg_lib:\n");
    printf("\toverstory: %d\n", vlib->overstory);
    printf("\tLAI:");
    for (i = 0; i < MONTHSPERYEAR; i++) {
        printf("\t%.2lf", vlib->LAI[i]);
    }
    printf("\n");
    printf("\tWdmax:");
    for (i = 0; i < MONTHSPERYEAR; i++) {
        printf("\t%.2lf", vlib->Wdmax[i]);
    }
    printf("\n");
    printf("\talbedo:");
    for (i = 0; i < MONTHSPERYEAR; i++) {
        printf("\t%.2lf", vlib->albedo[i]);
    }
    printf("\n");
    printf("\tdisplacement:");
    for (i = 0; i < MONTHSPERYEAR; i++) {
        printf("\t%.2lf", vlib->displacement[i]);
    }
    printf("\n");
    printf("\temissivity:");
    for (i = 0; i < MONTHSPERYEAR; i++) {
        printf("\t%.2lf", vlib->emissivity[i]);
    }
    printf("\n");
    printf("\tNVegLibTypes: %d\n", vlib->NVegLibTypes);
    printf("\trad_atten: %.4lf\n", vlib->rad_atten);
    printf("\trarc: %.4lf\n", vlib->rarc);
    printf("\trmin: %.4f\n", vlib->rmin);
    printf("\troughness:");
    for (i = 0; i < MONTHSPERYEAR; i++) {
        printf("\t%.2f", vlib->roughness[i]);
    }
    printf("\n");
    printf("\ttrunk_ratio: %.4lf\n", vlib->trunk_ratio);
    printf("\twind_atten: %.4lf\n", vlib->wind_atten);
    printf("\twind_h: %.4lf\n", vlib->wind_h);
    printf("\tRGL: %.4f\n", vlib->RGL);
    printf("\tveg_class: %d\n", vlib->veg_class);
    if (carbon) {
        printf("\tCtype: %d\n", vlib->Ctype);
        printf("\tMaxCarboxRate: %.4lf\n", vlib->MaxCarboxRate);
        printf("\tMaxETransport: %.4lf\n", vlib->MaxETransport);
        printf("\tCO2Specificity: %.4lf\n", vlib->CO2Specificity);
        printf("\tLightUseEff: %.4lf\n", vlib->LightUseEff);
        printf("\tNscaleFlag: %d\n", vlib->NscaleFlag);
        printf("\tWnpp_inhib: %.4lf\n", vlib->Wnpp_inhib);
        printf("\tNPPfactor_sat: %.4lf\n", vlib->NPPfactor_sat);
    }
}
