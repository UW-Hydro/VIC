/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine determines the current forcing file data type and stores its
 * location in the description of the current forcing file.
 *****************************************************************************/

#include <vic_driver_shared_image.h>

/******************************************************************************
 * @brief    This routine determines the current forcing file data type and
 *           stores its location in the description of the current forcing file.
 *****************************************************************************/
void
set_force_type(char *cmdstr,
               int   file_num,
               int  *field)
{
    extern param_set_struct param_set;

    char                    optstr[MAXSTRING];
    char                    flgstr[MAXSTRING];
    char                    ncvarname[MAXSTRING];
    int                     type = SKIP;

    strcpy(ncvarname, "MISSING");

    /** Initialize flgstr **/
    strcpy(flgstr, "NULL");

    if ((*field) >= (int) param_set.N_TYPES[file_num]) {
        log_err("Too many variables defined for forcing file %i., was "
                "expecting at most %zu and got %d", file_num + 1,
                param_set.N_TYPES[file_num], *field);
    }

    sscanf(cmdstr, "%*s %s %s", optstr, ncvarname);

    /***************************************
       Get meteorological data forcing info
    ***************************************/

    /* type 0: air temperature [C] */
    if (strcasecmp("AIR_TEMP", optstr) == 0) {
        type = AIR_TEMP;
    }
    /* type 1: albedo [fraction] */
    else if (strcasecmp("ALBEDO", optstr) == 0) {
        type = ALBEDO;
    }
    /* type 2: atmospheric CO2 mixing ratio [ppm] */
    else if (strcasecmp("CATM", optstr) == 0) {
        type = CATM;
    }
    /* type 3: incoming channel flow [m3] */
    else if (strcasecmp("CHANNEL_IN", optstr) == 0) {
        type = CHANNEL_IN;
    }
    /* type 4: vegetation cover fraction [fraction] */
    else if (strcasecmp("FCANOPY", optstr) == 0) {
        type = FCANOPY;
    }
    /* type 5: direct fraction of shortwave [fraction] */
    else if (strcasecmp("FDIR", optstr) == 0) {
        type = FDIR;
    }
    /* type 6: LAI [m2/m2] */
    else if (strcasecmp("LAI", optstr) == 0) {
        type = LAI;
    }
    /* type 7: incoming longwave radiation [W/m2] */
    else if (strcasecmp("LWDOWN", optstr) == 0) {
        type = LWDOWN;
    }
    /* type 8: photosynthetically active radiation [uE/m2s] */
    else if (strcasecmp("PAR", optstr) == 0) {
        type = PAR;
    }
    /* type 9: precipitation [mm] */
    else if (strcasecmp("PREC", optstr) == 0) {
        type = PREC;
    }
    /* type 10: air pressure [kPa] */
    else if (strcasecmp("PRESSURE", optstr) == 0) {
        type = PRESSURE;
    }
    /* type 11: vapor pressure [kPa] */
    else if (strcasecmp("VP", optstr) == 0) {
        type = VP;
    }
    /* type 12: incoming shortwave radiation [W/m2] */
    else if (strcasecmp("SWDOWN", optstr) == 0) {
        type = SWDOWN;
    }
    /* type 13: wind speed [m/s] */
    else if (strcasecmp("WIND", optstr) == 0) {
        type = WIND;
    }
    /* type 14: unused (blank) data */
    else if (strcasecmp("SKIP", optstr) == 0) {
        type = SKIP;
    }
    /** Undefined variable type **/
    else {
        log_err("Undefined forcing variable type %s in file %i.",
                optstr, file_num + 1);
    }

    param_set.TYPE[type].SUPPLIED = file_num + 1;
    param_set.FORCE_INDEX[file_num][(*field)] = type;

    if (strcasecmp("MISSING", ncvarname) != 0) {
        strcpy(param_set.TYPE[type].varname, ncvarname);
    }
    else {
        log_err(
            "Must supply netCDF variable name for %s forcing file number %d",
            optstr, file_num + 1);
    }

    param_set.TYPE[type].N_ELEM = 1;

    (*field)++;
}
