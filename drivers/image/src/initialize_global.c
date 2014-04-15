#include <vic_def.h>
#include <vic_run.h>
#include <vic_driver_image.h>

void
initialize_global()
{
/*********************************************************************
   initialize_global              Keith Cherkauer       March 1998

   This subroutine initalizes all global parameters before they are
   called by the model.
*********************************************************************/

    extern option_struct    options;
    extern param_set_struct param_set;

    int                     i, j;

    /** Initialize model option flags **/

    // simulation modes
    options.AboveTreelineVeg = -1;
    options.AERO_RESIST_CANSNOW = AR_406_FULL;
    options.BLOWING = FALSE;
    options.CARBON = FALSE;
    options.CLOSE_ENERGY = FALSE;
    options.COMPUTE_TREELINE = FALSE;
    options.CONTINUEONERROR = TRUE;
    options.CORRPREC = FALSE;
    options.EQUAL_AREA = FALSE;
    options.EXP_TRANS = TRUE;
    options.FROZEN_SOIL = FALSE;
    options.FULL_ENERGY = FALSE;
    options.GRND_FLUX_TYPE = GF_410;
    options.IMPLICIT = TRUE;
    options.LAKES = FALSE;
    options.LAKE_PROFILE = FALSE;
    options.LOG_MATRIC = FALSE;
    options.LW_CLOUD = LW_CLOUD_DEARDORFF;
    options.LW_TYPE = LW_PRATA;
    options.MIN_WIND_SPEED = 0.1;
    options.MTCLIM_SWE_CORR = FALSE;
    options.Ncanopy = 3;
    options.Nfrost = 1;
    options.Nlayer = 3;
    options.Nnode = 3;
    options.NOFLUX = FALSE;
    options.PLAPSE = TRUE;
    options.QUICK_FLUX = TRUE;
    options.QUICK_SOLVE = FALSE;
    options.RC_MODE = RC_JARVIS;
    options.ROOT_ZONES = MISSING;
    options.SHARE_LAYER_MOIST = TRUE;
    options.SNOW_BAND = 1;
    options.SNOW_DENSITY = DENS_BRAS;
    options.SNOW_STEP = 1;
    options.SPATIAL_FROST = FALSE;
    options.SPATIAL_SNOW = FALSE;
    options.SW_PREC_THRESH = 0;
    options.TFALLBACK = TRUE;
    options.VP_INTERP = TRUE;
    options.VP_ITER = VP_ITER_ALWAYS;
    // input options
    options.BASEFLOW = ARNO;
    options.GRID_DECIMAL = 2;
    options.JULY_TAVG_SUPPLIED = FALSE;
    options.ORGANIC_FRACT = FALSE;
    options.VEGLIB_PHOTO = FALSE;
    options.VEGPARAM_LAI = FALSE;
    options.LAI_SRC = LAI_FROM_VEGLIB;
    // state options
    options.BINARY_STATE_FILE = FALSE;
    options.INIT_STATE = FALSE;
    options.SAVE_STATE = FALSE;
    // output options
    options.ALMA_OUTPUT = FALSE;
    options.BINARY_OUTPUT = FALSE;
    options.COMPRESS = FALSE;
    options.MOISTFRACT = FALSE;
    options.Noutfiles = 2;
    options.OUTPUT_FORCE = FALSE;
    options.PRT_HEADER = FALSE;
    options.PRT_SNOW_BAND = FALSE;

    /** Initialize forcing file input controls **/

    for (j = 0; j < N_FORCING_TYPES; j++) {
        param_set.TYPE[j].SUPPLIED = FALSE;
        param_set.TYPE[j].SIGNED = 1;
        param_set.TYPE[j].multiplier = 1;
    }
    for (i = 0; i < 2; i++) {
        param_set.FORCE_DT[i] = MISSING;
        param_set.N_TYPES[i] = MISSING;
        param_set.FORCE_FORMAT[i] = MISSING;
        for (j = 0; j < N_FORCING_TYPES; j++) {
            param_set.FORCE_INDEX[i][j] = MISSING;
        }
    }
}
