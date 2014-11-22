#include <vic_def.h>
#include <vic_run.h>
#include <vic_driver_image.h>

void
alloc_veg_hist(veg_hist_struct *veg_hist)
{
    veg_hist->albedo = (double *) calloc(NR + 1, sizeof(double));
    if (veg_hist->albedo == NULL) {
        nrerror("Memory allocation error in alloc_veg_hist().");
    }
    veg_hist->LAI = (double *) calloc(NR + 1, sizeof(double));
    if (veg_hist->LAI == NULL) {
        nrerror("Memory allocation error in alloc_veg_hist().");
    }
    veg_hist->vegcover = (double *) calloc(NR + 1, sizeof(double));
    if (veg_hist->vegcover == NULL) {
        nrerror("Memory allocation error in alloc_veg_hist().");
    }
}

void
free_veg_hist(veg_hist_struct *veg_hist)
{
    if (veg_hist == NULL) {
        return;
    }

    free(veg_hist->albedo);
    free(veg_hist->LAI);
    free(veg_hist->vegcover);
}
