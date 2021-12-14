/******************************************************************************
 * @section DESCRIPTION
 *
 * Allocate and free memory for the veg_hist data struct
 *****************************************************************************/

#include <vic_driver_classic.h>

/******************************************************************************
 * @brief    Allocate memory for veg his structure.
 *****************************************************************************/
void
alloc_veg_hist(int                nrecs,
               int                nveg,
               veg_hist_struct ***veg_hist)
{
    int i, j;

    *veg_hist = calloc(nrecs, sizeof(*(*veg_hist)));
    check_alloc_status((*veg_hist), "Memory allocation error.");

    for (i = 0; i < nrecs; i++) {
        (*veg_hist)[i] = calloc(nveg + 1, sizeof(*((*veg_hist)[i])));
        check_alloc_status((*veg_hist)[i], "Memory allocation error.");

        for (j = 0; j < nveg + 1; j++) {
            (*veg_hist)[i][j].albedo =
                calloc(NR + 1, sizeof(*((*veg_hist)[i][j].albedo)));
            check_alloc_status((*veg_hist)[i][j].albedo,
                               "Memory allocation error.");
            (*veg_hist)[i][j].displacement = calloc(NR + 1,
                                                    sizeof(*((*veg_hist)[i][j].
                                                             displacement)));
            check_alloc_status((*veg_hist)[i][j].displacement,
                               "Memory allocation error.");
            (*veg_hist)[i][j].fcanopy = calloc(NR + 1,
                                               sizeof(*((*veg_hist)[i][j].
                                                        fcanopy)));
            check_alloc_status((*veg_hist)[i][j].fcanopy,
                               "Memory allocation error.");
            (*veg_hist)[i][j].LAI =
                calloc(NR + 1, sizeof(*((*veg_hist)[i][j].LAI)));
            check_alloc_status((*veg_hist)[i][j].LAI,
                               "Memory allocation error.");
            (*veg_hist)[i][j].roughness = calloc(NR + 1,
                                                 sizeof(*((*veg_hist)[i][j].
                                                          roughness)));
            check_alloc_status((*veg_hist)[i][j].roughness,
                               "Memory allocation error.");
        }
    }
}

/******************************************************************************
 * @brief    Free veg hist structure.
 *****************************************************************************/
void
free_veg_hist(int                nrecs,
              int                nveg,
              veg_hist_struct ***veg_hist)
{
    int i, j;

    if (*veg_hist == NULL) {
        return;
    }

    for (i = 0; i < nrecs; i++) {
        for (j = 0; j < nveg + 1; j++) {
            free((*veg_hist)[i][j].albedo);
            free((*veg_hist)[i][j].displacement);
            free((*veg_hist)[i][j].fcanopy);
            free((*veg_hist)[i][j].LAI);
            free((*veg_hist)[i][j].roughness);
        }
        free((*veg_hist)[i]);
    }

    free(*veg_hist);
}
