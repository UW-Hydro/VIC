/******************************************************************************
 * @section DESCRIPTION
 *
 * Allocate and free memory for the veg_hist data struct
 *****************************************************************************/

 #include <vic_driver_shared_image.h>

/******************************************************************************
 * @brief
 *****************************************************************************/
void
alloc_veg_hist(veg_hist_struct *veg_hist)
{
    veg_hist->albedo = calloc(NR + 1, sizeof(*(veg_hist->albedo)));
    check_alloc_status(veg_hist->albedo, "Memory allocation error.");

    veg_hist->displacement = calloc(NR + 1, sizeof(*(veg_hist->displacement)));
    check_alloc_status(veg_hist->displacement, "Memory allocation error.");

    veg_hist->fcanopy = calloc(NR + 1, sizeof(*(veg_hist->fcanopy)));
    check_alloc_status(veg_hist->fcanopy, "Memory allocation error.");

    veg_hist->LAI = calloc(NR + 1, sizeof(*(veg_hist->LAI)));
    check_alloc_status(veg_hist->LAI, "Memory allocation error.");

    veg_hist->roughness = calloc(NR + 1, sizeof(*(veg_hist->roughness)));
    check_alloc_status(veg_hist->roughness, "Memory allocation error.");
}

/******************************************************************************
 * @brief    Free veg hist structure.
 *****************************************************************************/
void
free_veg_hist(veg_hist_struct *veg_hist)
{
    if (veg_hist == NULL) {
        return;
    }

    free(veg_hist->albedo);
    free(veg_hist->displacement);
    free(veg_hist->fcanopy);
    free(veg_hist->LAI);
    free(veg_hist->roughness);
}
