/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine reads in a library of vegetation parameters for all vegetation
 * classes used in the model.  The veg class number is used to reference the
 * information in this library.
 *
 * @section LICENSE
 *
 * The Variable Infiltration Capacity (VIC) macroscale hydrological model
 * Copyright (C) 2016 The Computational Hydrology Group, Department of Civil
 * and Environmental Engineering, University of Washington.
 *
 * The VIC model is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
 *****************************************************************************/

#include <vic_driver_classic.h>

/******************************************************************************
 * @brief    Read in a library of vegetation parameters for all vegetation
 *           classes used in the model.
 *****************************************************************************/
veg_lib_struct *
read_veglib(FILE   *veglib,
            size_t *Ntype)
{
    extern option_struct       options;
    extern parameters_struct   param;
    extern global_param_struct global_param;

    veg_lib_struct            *temp;
    size_t                     i, j;
    int                        tmpflag;
    size_t                     Nveg_type;
    char                       str[MAXSTRING];
    double                     maxd;
    char                       tmpstr[MAXSTRING];
    double                     tmp_double;

    rewind(veglib);
    fgets(str, MAXSTRING, veglib);
    Nveg_type = 0;
    while (!feof(veglib)) {
        if (str[0] <= 57 && str[0] >= 48) {
            Nveg_type++;
        }
        fgets(str, MAXSTRING, veglib);
    }
    rewind(veglib);

    // +1 for bare soil
    temp = calloc(Nveg_type + 1, sizeof(*temp));
    options.NVEGTYPES = Nveg_type + 1;

    fscanf(veglib, "%s", str);
    i = 0;
    while (!feof(veglib)) {
        if (str[0] <= 57 && str[0] >= 48) {
            temp[i].NVegLibTypes = Nveg_type;
            temp[i].veg_class = atoi(str);
            fscanf(veglib, "%d", &tmpflag);
            if (tmpflag == 0) {
                temp[i].overstory = false;
            }
            else {
                temp[i].overstory = true;
            }
            fscanf(veglib, "%lf", &temp[i].rarc);
            fscanf(veglib, "%lf", &temp[i].rmin);
            for (j = 0; j < MONTHS_PER_YEAR; j++) {
                fscanf(veglib, "%lf", &temp[i].LAI[j]);
                if (options.LAI_SRC == FROM_VEGLIB && temp[i].overstory &&
                    temp[i].LAI[j] == 0) {
                    log_err("veg library: the specified veg class (%d) "
                            "is listed as an overstory class, but the LAI "
                            "given for this class for month %zu is 0",
                            temp[i].veg_class, j);
                }
                temp[i].Wdmax[j] = param.VEG_LAI_WATER_FACTOR * temp[i].LAI[j];
            }
            /* Default values of fcanopy */
            for (j = 0; j < MONTHS_PER_YEAR; j++) {
                temp[i].fcanopy[j] = 1.00;
            }
            if (options.VEGLIB_FCAN) {
                for (j = 0; j < MONTHS_PER_YEAR; j++) {
                    fscanf(veglib, "%lf", &tmp_double);
                    if (options.FCAN_SRC != FROM_DEFAULT) {
                        temp[i].fcanopy[j] = tmp_double;
                        if (temp[i].fcanopy[j] < 0 ||
                            temp[i].fcanopy[j] > 1) {
                            log_err(
                                "Veg cover fraction must be between 0 and 1 " "(%f)",
                                temp[i].fcanopy[j]);
                        }
                    }
                }
            }
            for (j = 0; j < MONTHS_PER_YEAR; j++) {
                fscanf(veglib, "%lf", &temp[i].albedo[j]);
                if (temp[i].albedo[j] < 0 || temp[i].albedo[j] > 1) {
                    log_err("Albedo must be between 0 and 1 (%f)",
                            temp[i].albedo[j]);
                }
            }
            for (j = 0; j < MONTHS_PER_YEAR; j++) {
                fscanf(veglib, "%lf", &temp[i].roughness[j]);
            }
            temp[i].wind_h = 0.;
            maxd = 0;
            for (j = 0; j < MONTHS_PER_YEAR; j++) {
                fscanf(veglib, "%lf", &temp[i].displacement[j]);
                if (temp[i].displacement[j] > maxd) {
                    maxd = temp[i].displacement[j];
                }
                if (temp[i].LAI[j] > 0 && temp[i].displacement[j] <= 0) {
                    log_err("Vegetation has leaves (LAI = %f), but no "
                            "displacement (%f)",
                            temp[i].LAI[j], temp[i].displacement[j]);
                }
            }
            fscanf(veglib, "%lf", &temp[i].wind_h);
            if (temp[i].wind_h < maxd && temp[i].overstory) {
                log_err("Vegetation reference height (%f) for vegetation "
                        "class %d, must be greater than the maximum "
                        "displacement height (%f) when OVERSTORY has been set "
                        "true.", temp[i].wind_h, temp[i].veg_class, maxd);
            }
            fscanf(veglib, "%lf", &temp[i].RGL);    /* minimum value of incoming
                                                       solar radiation at which there
                                                       will still be transpiration */
            if (temp[i].RGL < 0) {
                log_err("Minimum value of incoming solar radiation at which "
                        "there is transpiration (RGL) must be greater than 0 "
                        "for vegetation class %d.  Check that the vegetation "
                        "library has the correct number of columns.",
                        temp[i].veg_class);
            }
            fscanf(veglib, "%lf", &temp[i].rad_atten); /* vegetation radiation
                                                          attenuation factor */
            if (temp[i].rad_atten < 0 || temp[i].rad_atten > 1) {
                log_err("The vegetation radiation attenuation factor must be "
                        "greater than 0, and less than 1 for vegetation class "
                        "%d.  Check that the vegetation library has the "
                        "correct number of columns.", temp[i].veg_class);
            }
            fscanf(veglib, "%lf", &temp[i].wind_atten); /* canopy wind speed
                                                           attenuation factor */
            fscanf(veglib, "%lf", &temp[i].trunk_ratio); /* ratio of tree height that
                                                            is trunk */
            /* Carbon-cycling parameters */
            if (options.VEGLIB_PHOTO) {
                fscanf(veglib, "%s", tmpstr); /* photosynthetic pathway */
                if (!strcmp(tmpstr, "0") || !strcmp(tmpstr, "C3")) {
                    temp[i].Ctype = PHOTO_C3;
                }
                else if (!strcmp(tmpstr, "1") || !strcmp(tmpstr, "C4")) {
                    temp[i].Ctype = PHOTO_C4;
                }
                if (!strcmp(tmpstr, "C3") || !strcmp(tmpstr, "C4")) {
                    log_warn("Use of strings \"C3\" and \"C4\" as values of "
                             "Ctype is deprecated.  Please replace these with "
                             "\"0\" and \"1\", respectively");
                }
                fscanf(veglib, "%lf", &temp[i].MaxCarboxRate); /* Maximum carboxylation rate at 25 deg C */
                if (temp[i].Ctype == PHOTO_C3) {
                    fscanf(veglib, "%lf", &temp[i].MaxETransport); /* Maximum electron transport rate at 25 deg C */
                    temp[i].CO2Specificity = 0;
                }
                else if (temp[i].Ctype == PHOTO_C4) {
                    fscanf(veglib, "%lf", &temp[i].CO2Specificity); /* CO2 Specificity */
                    temp[i].MaxETransport = 0;
                }
                fscanf(veglib, "%lf", &temp[i].LightUseEff); /* Light-use efficiency */
                fscanf(veglib, "%s", tmpstr); /* Nitrogen-scaling flag */
                temp[i].NscaleFlag = atoi(tmpstr); /* Nitrogen-scaling flag */
                fscanf(veglib, "%lf", &temp[i].Wnpp_inhib); /* Moisture level in top soil layer above which photosynthesis begins experiencing inhibition due to saturation */
                fscanf(veglib, "%lf", &temp[i].NPPfactor_sat); /* photosynthesis multiplier when top soil layer is saturated */
            }
            else {
                temp[i].Wnpp_inhib = 1.0;
                temp[i].NPPfactor_sat = 1.0;
            }

            fgets(str, MAXSTRING, veglib); /* skip over end of line comments */
            i++;
        }
        else {
            fgets(str, MAXSTRING, veglib);
        }
        fscanf(veglib, "%s", str);
    }
    if (i != Nveg_type) {
        log_err("Problem reading vegetation library file - make sure "
                "the file has the right number of columns.");
    }
    *Ntype = Nveg_type;

    // Assign properties of bare soil to default bare soil tile
    temp[i].NVegLibTypes = Nveg_type;
    temp[i].veg_class = Nveg_type + 1;
    temp[i].overstory = false;
    temp[i].rarc = param.SOIL_RARC;
    temp[i].rmin = 0.0;
    for (j = 0; j < MONTHS_PER_YEAR; j++) {
        temp[i].LAI[j] = 0.0;
        temp[i].Wdmax[j] = 0.0;
        temp[i].fcanopy[j] = MIN_FCANOPY;
        temp[i].albedo[j] = param.ALBEDO_BARE_SOIL;
        // These will be assigned in read_soilparam.c
        temp[i].roughness[j] = MISSING;
        temp[i].displacement[j] = MISSING;
    }
    temp[i].wind_h = global_param.wind_h;
    temp[i].RGL = 0.0;
    temp[i].rad_atten = 0.0;
    temp[i].wind_atten = 0.0;
    temp[i].trunk_ratio = 0.0;
    if (options.VEGLIB_PHOTO) {
        temp[i].Ctype = PHOTO_C3;
        temp[i].MaxETransport = 0.0;
        temp[i].CO2Specificity = 0.0;
        temp[i].LightUseEff = 0.0;
        temp[i].NscaleFlag = 0;
        temp[i].Wnpp_inhib = 1.0;
        temp[i].NPPfactor_sat = 1.0;
    }

    return temp;
}

/******************************************************************************
 * @brief    This routine frees the veglib structure.
 *****************************************************************************/
void
free_veglib(veg_lib_struct **veg_lib)
{
    free((char*)(*veg_lib));
}
