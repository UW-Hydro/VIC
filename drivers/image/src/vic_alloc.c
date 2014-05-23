#include <vic_def.h>
#include <vic_run.h>
#include <vic_driver_image.h>

void
vic_alloc()
{
    extern domain_struct       global_domain;
    extern filenames_struct    filenames;
    extern option_struct       options;
    extern soil_con_struct    *soil_con;
    extern veg_con_map_struct *veg_con_map;
    extern veg_con_struct    **veg_con;
    extern veg_lib_struct    **veg_lib;

    double                    *dvar = NULL;

    size_t                     i;
    size_t                     j;
    size_t                     idx;
    size_t                     d2count[2];
    size_t                     d2start[2];

    // TBD: handle decomposed domain

    // allocate memory for number of veg types to be read
    // this is to maintain consistency with the VIC data structures, so that
    // the call to vic_run() is not affected
    dvar = (double *) malloc(global_domain.n_ny * global_domain.n_nx *
                             sizeof(double));
    if (dvar == NULL) {
        nrerror("Memory allocation error in vic_init().");
    }
    // read the number of vegetation types for each grid cell. This result is
    // not stored at this time, but will be read again in vic_init(). Clunky,
    // but so it is
    // overstory
    d2start[0] = 0;
    d2start[1] = 0;
    d2count[0] = global_domain.n_ny;
    d2count[1] = global_domain.n_nx;
    get_nc_field_double(filenames.veglib, "Nveg", d2start, d2count, dvar);

    // allocate memory for soil structure
    soil_con = (soil_con_struct *)
               malloc((size_t) global_domain.ncells_global *
                      sizeof(soil_con_struct));
    if (soil_con == NULL) {
        nrerror("Memory allocation error in vic_alloc().");
    }

    // allocate memory for vegetation mapping structure
    veg_con_map = (veg_con_map_struct *)
                  malloc((size_t) global_domain.ncells_global *
                         sizeof(veg_con_map_struct));
    if (veg_con_map == NULL) {
        nrerror("Memory allocation error in vic_alloc().");
    }

    // allocate memory for vegetation structure
    veg_con = (veg_con_struct **)
              malloc((size_t) global_domain.ncells_global *
                     sizeof(veg_con_struct *));
    if (veg_con == NULL) {
        nrerror("Memory allocation error in vic_alloc().");
    }

    // allocate memory for vegetation structure
    veg_lib = (veg_lib_struct **)
              malloc((size_t) global_domain.ncells_global *
                     sizeof(veg_lib_struct *));
    if (veg_lib == NULL) {
        nrerror("Memory allocation error in vic_alloc().");
    }

    // allocate memory for individual grid cells
    for (i = 0; i < global_domain.ncells_global; i++) {
        // snow band allocation
        soil_con[i].AreaFract = (double *) calloc(options.SNOW_BAND,
                                                  sizeof(double));
        if (soil_con[i].AreaFract == NULL) {
            nrerror("Memory allocation error in vic_alloc().");
        }
        soil_con[i].BandElev = (float *) calloc(options.SNOW_BAND,
                                                sizeof(float));
        if (soil_con[i].BandElev == NULL) {
            nrerror("Memory allocation error in vic_alloc().");
        }
        soil_con[i].Tfactor = (double *) calloc(options.SNOW_BAND,
                                                sizeof(double));
        if (soil_con[i].Tfactor == NULL) {
            nrerror("Memory allocation error in vic_alloc().");
        }
        soil_con[i].Pfactor = (double *) calloc(options.SNOW_BAND,
                                                sizeof(double));
        if (soil_con[i].Pfactor == NULL) {
            nrerror("Memory allocation error in vic_alloc().");
        }
        soil_con[i].AboveTreeLine = (char *) calloc(options.SNOW_BAND,
                                                    sizeof(char));
        if (soil_con[i].AboveTreeLine == NULL) {
            nrerror("Memory allocation error in vic_alloc().");
        }

        // vegetation tile allocation

        veg_con_map[i].nv_types = options.NVEGTYPES;

        veg_con_map[i].vidx = (int *) malloc((size_t) veg_con_map[i].nv_types *
                                             sizeof(int));
        if (veg_con_map[i].vidx == NULL) {
            nrerror("Memory allocation error in vic_alloc().");
        }
        veg_con_map[i].Cv = (double *) malloc((size_t) veg_con_map[i].nv_types *
                                              sizeof(double));
        if (veg_con_map[i].Cv == NULL) {
            nrerror("Memory allocation error in vic_alloc().");
        }

        idx = global_domain.locations[i].global_y_idx * global_domain.n_nx +
              global_domain.locations[i].global_x_idx;
        veg_con_map[i].nv_active = (size_t) dvar[idx] + 1;
        if (options.AboveTreelineVeg >= 0) {
            veg_con_map[i].nv_active += 1;
        }

        veg_con[i] = (veg_con_struct *)
                     malloc((size_t) (veg_con_map[i].nv_active) *
                            sizeof(veg_con_struct));
        if (veg_con[i] == NULL) {
            nrerror("Memory allocation error in vic_alloc().");
        }

        for (j = 0; j < veg_con_map[i].nv_active; j++) {
            veg_con[i][j].zone_depth = calloc(options.ROOT_ZONES,
                                              sizeof(float));
            if (veg_con[i][j].zone_depth == NULL) {
                nrerror("Memory allocation error in vic_alloc().");
            }
            veg_con[i][j].zone_fract = calloc(options.ROOT_ZONES,
                                              sizeof(float));
            if (veg_con[i][j].zone_fract == NULL) {
                nrerror("Memory allocation error in vic_alloc().");
            }
            if (options.CARBON) {
                veg_con[i][j].CanopLayerBnd = calloc(options.Ncanopy,
                                                     sizeof(double));
                if (veg_con[i][j].CanopLayerBnd == NULL) {
                    nrerror("Memory allocation error in vic_alloc().");
                }
            }
        }

        // vegetation library allocation - there is a veg library for each
        // active grid cell

        veg_lib[i] = (veg_lib_struct *) calloc(options.NVEGTYPES,
                                               sizeof(veg_lib_struct));
        if (veg_lib[i] == NULL) {
            nrerror("Memory allocation error in vic_alloc().");
        }
    }
    // cleanup
    free(dvar);
}
