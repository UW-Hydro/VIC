#include <vic_def.h>
#include <vic_run.h>
#include <vic_driver_image.h>

void
alloc_atmos(atmos_data_struct *atmos)
{
    atmos->air_temp = (double *) calloc(NR + 1, sizeof(double));
    if (atmos->air_temp == NULL) {
        nrerror("Memory allocation error in alloc_atmos().");
    }
    atmos->Catm = (double *) calloc(NR + 1, sizeof(double));
    if (atmos->Catm == NULL) {
        nrerror("Memory allocation error in alloc_atmos().");
    }
    atmos->channel_in = (double *) calloc(NR + 1, sizeof(double));
    if (atmos->channel_in == NULL) {
        nrerror("Memory allocation error in alloc_atmos().");
    }
    atmos->coszen = (double *) calloc(NR + 1, sizeof(double));
    if (atmos->coszen == NULL) {
        nrerror("Memory allocation error in alloc_atmos().");
    }
    atmos->density = (double *) calloc(NR + 1, sizeof(double));
    if (atmos->density == NULL) {
        nrerror("Memory allocation error in alloc_atmos().");
    }
    atmos->fdir = (double *) calloc(NR + 1, sizeof(double));
    if (atmos->fdir == NULL) {
        nrerror("Memory allocation error in alloc_atmos().");
    }
    atmos->longwave = (double *) calloc(NR + 1, sizeof(double));
    if (atmos->longwave == NULL) {
        nrerror("Memory allocation error in alloc_atmos().");
    }
    atmos->par = (double *) calloc(NR + 1, sizeof(double));
    if (atmos->par == NULL) {
        nrerror("Memory allocation error in alloc_atmos().");
    }
    atmos->prec = (double *) calloc(NR + 1, sizeof(double));
    if (atmos->prec == NULL) {
        nrerror("Memory allocation error in alloc_atmos().");
    }
    atmos->pressure = (double *) calloc(NR + 1, sizeof(double));
    if (atmos->pressure == NULL) {
        nrerror("Memory allocation error in alloc_atmos().");
    }
    atmos->shortwave = (double *) calloc(NR + 1, sizeof(double));
    if (atmos->shortwave == NULL) {
        nrerror("Memory allocation error in alloc_atmos().");
    }
    atmos->snowflag = (char *) calloc(NR + 1, sizeof(char));
    if (atmos->snowflag == NULL) {
        nrerror("Memory allocation error in alloc_atmos().");
    }
    atmos->tskc = (double *) calloc(NR + 1, sizeof(double));
    if (atmos->tskc == NULL) {
        nrerror("Memory allocation error in alloc_atmos().");
    }
    atmos->vp = (double *) calloc(NR + 1, sizeof(double));
    if (atmos->vp == NULL) {
        nrerror("Memory allocation error in alloc_atmos().");
    }
    atmos->vpd = (double *) calloc(NR + 1, sizeof(double));
    if (atmos->vpd == NULL) {
        nrerror("Memory allocation error in alloc_atmos().");
    }
    atmos->wind = (double *) calloc(NR + 1, sizeof(double));
    if (atmos->wind == NULL) {
        nrerror("Memory allocation error in alloc_atmos().");
    }
}

void
free_atmos(atmos_data_struct *atmos)
{
    if (atmos == NULL) {
        return;
    }

    free(atmos->air_temp);
    free(atmos->Catm);
    free(atmos->channel_in);
    free(atmos->coszen);
    free(atmos->density);
    free(atmos->fdir);
    free(atmos->longwave);
    free(atmos->par);
    free(atmos->prec);
    free(atmos->pressure);
    free(atmos->shortwave);
    free(atmos->snowflag);
    free(atmos->tskc);
    free(atmos->vp);
    free(atmos->vpd);
    free(atmos->wind);
}
