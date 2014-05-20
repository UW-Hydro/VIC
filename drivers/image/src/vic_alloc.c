#include <vic_def.h>
#include <vic_run.h>
#include <vic_driver_image.h>

void
vic_alloc()
{
    extern domain_struct    global_domain;
    extern option_struct    options;
    extern soil_con_struct *soil_con;
    extern veg_con_struct **veg_con;
    extern veg_lib_struct **veg_lib;

    size_t                  i;
    size_t                  j;


    // TBD: handle decomposed domain

    // allocate memory for soil structure
    soil_con = (soil_con_struct *)
               malloc((size_t) global_domain.ncells_global *
                      sizeof(soil_con_struct));
    if (soil_con == NULL) {
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

        veg_con[i] = (veg_con_struct *)
                     malloc((size_t) options.NVEGTYPES *
                            sizeof(veg_con_struct));
        if (veg_con[i] == NULL) {
            nrerror("Memory allocation error in vic_alloc().");
        }

        for (j = 0; j < options.NVEGTYPES; j++) {
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
}
