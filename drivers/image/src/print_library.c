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
