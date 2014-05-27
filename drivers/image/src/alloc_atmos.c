#include <vic_def.h>
#include <vic_run.h>
#include <vic_driver_image.h>

void
alloc_atmos(int                 nrecs,
            atmos_data_struct **atmos)
{
    extern int NR;

    int        i;

    *atmos = (atmos_data_struct *) calloc(nrecs, sizeof(atmos_data_struct));
    if (*atmos == NULL) {
        nrerror("Memory allocation error in alloc_atmos().");
    }

    for (i = 0; i < nrecs; i++) {
        (*atmos)[i].air_temp = (double *) calloc(NR + 1, sizeof(double));
        if ((*atmos)[i].air_temp == NULL) {
            nrerror("Memory allocation error in alloc_atmos().");
        }
        (*atmos)[i].Catm = (double *) calloc(NR + 1, sizeof(double));
        if ((*atmos)[i].Catm == NULL) {
            nrerror("Memory allocation error in alloc_atmos().");
        }
        (*atmos)[i].channel_in = (double *) calloc(NR + 1, sizeof(double));
        if ((*atmos)[i].channel_in == NULL) {
            nrerror("Memory allocation error in alloc_atmos().");
        }
        (*atmos)[i].coszen = (double *) calloc(NR + 1, sizeof(double));
        if ((*atmos)[i].coszen == NULL) {
            nrerror("Memory allocation error in alloc_atmos().");
        }
        (*atmos)[i].density = (double *) calloc(NR + 1, sizeof(double));
        if ((*atmos)[i].density == NULL) {
            nrerror("Memory allocation error in alloc_atmos().");
        }
        (*atmos)[i].fdir = (double *) calloc(NR + 1, sizeof(double));
        if ((*atmos)[i].fdir == NULL) {
            nrerror("Memory allocation error in alloc_atmos().");
        }
        (*atmos)[i].longwave = (double *) calloc(NR + 1, sizeof(double));
        if ((*atmos)[i].longwave == NULL) {
            nrerror("Memory allocation error in alloc_atmos().");
        }
        (*atmos)[i].par = (double *) calloc(NR + 1, sizeof(double));
        if ((*atmos)[i].par == NULL) {
            nrerror("Memory allocation error in alloc_atmos().");
        }
        (*atmos)[i].prec = (double *) calloc(NR + 1, sizeof(double));
        if ((*atmos)[i].prec == NULL) {
            nrerror("Memory allocation error in alloc_atmos().");
        }
        (*atmos)[i].pressure = (double *) calloc(NR + 1, sizeof(double));
        if ((*atmos)[i].pressure == NULL) {
            nrerror("Memory allocation error in alloc_atmos().");
        }
        (*atmos)[i].shortwave = (double *) calloc(NR + 1, sizeof(double));
        if ((*atmos)[i].shortwave == NULL) {
            nrerror("Memory allocation error in alloc_atmos().");
        }
        (*atmos)[i].snowflag = (char *) calloc(NR + 1, sizeof(char));
        if ((*atmos)[i].snowflag == NULL) {
            nrerror("Memory allocation error in alloc_atmos().");
        }
        (*atmos)[i].tskc = (double *) calloc(NR + 1, sizeof(double));
        if ((*atmos)[i].tskc == NULL) {
            nrerror("Memory allocation error in alloc_atmos().");
        }
        (*atmos)[i].vp = (double *) calloc(NR + 1, sizeof(double));
        if ((*atmos)[i].vp == NULL) {
            nrerror("Memory allocation error in alloc_atmos().");
        }
        (*atmos)[i].vpd = (double *) calloc(NR + 1, sizeof(double));
        if ((*atmos)[i].vpd == NULL) {
            nrerror("Memory allocation error in alloc_atmos().");
        }
        (*atmos)[i].wind = (double *) calloc(NR + 1, sizeof(double));
        if ((*atmos)[i].wind == NULL) {
            nrerror("Memory allocation error in alloc_atmos().");
        }
    }
}

void
free_atmos(int                 nrecs,
           atmos_data_struct **atmos)
{
    int i;

    if (*atmos == NULL) {
        return;
    }

    for (i = 0; i < nrecs; i++) {
        free((*atmos)[i].air_temp);
        free((*atmos)[i].Catm);
        free((*atmos)[i].channel_in);
        free((*atmos)[i].coszen);
        free((*atmos)[i].density);
        free((*atmos)[i].fdir);
        free((*atmos)[i].longwave);
        free((*atmos)[i].par);
        free((*atmos)[i].prec);
        free((*atmos)[i].pressure);
        free((*atmos)[i].shortwave);
        free((*atmos)[i].snowflag);
        free((*atmos)[i].tskc);
        free((*atmos)[i].vp);
        free((*atmos)[i].vpd);
        free((*atmos)[i].wind);
    }

    free(*atmos);
}
