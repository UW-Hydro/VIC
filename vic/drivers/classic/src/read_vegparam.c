/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine reads in vegetation parameters for the current grid cell.
 *
 * It also relates each vegetation class in the cell to the appropriate
 * parameters in the vegetation library.
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
 * @brief    Read vegetation parameters.
 *****************************************************************************/
veg_con_struct *
read_vegparam(FILE  *vegparam,
              int    gridcel,
              size_t Nveg_type)
{
    void ttrim(char *string);
    extern veg_lib_struct   *veg_lib;
    extern option_struct     options;
    extern parameters_struct param;

    veg_con_struct          *temp;
    size_t                   j;
    int                      vegetat_type_num;
    int                      vegcel, i, k, skip, veg_class;
    int                      MaxVeg;
    int                      Nfields, NfieldsMax;
    int                      NoOverstory;
    double                   depth_sum;
    double                   sum;
    double                   Cv_sum;
    char                     str[MAX_VEGPARAM_LINE_LENGTH];
    char                     line[MAXSTRING];
    char                     tmpline[MAXSTRING];
    const char               delimiters[] = " \t";
    char                    *token;
    char                    *vegarr[MAX_VEGPARAM_LINE_LENGTH];
    size_t                   length;
    size_t                   cidx;
    double                   tmp;

    skip = 1;
    if (options.VEGPARAM_LAI) {
        skip++;
    }
    if (options.VEGPARAM_FCAN) {
        skip++;
    }
    if (options.VEGPARAM_ALB) {
        skip++;
    }

    NoOverstory = 0;

    while ((fscanf(vegparam, "%d %d", &vegcel,
                   &vegetat_type_num) == 2) && vegcel != gridcel) {
        if (vegetat_type_num < 0) {
            log_err("number of vegetation tiles (%i) given for cell %i "
                    "is < 0.", vegetat_type_num, vegcel);
        }
        for (i = 0; i <= vegetat_type_num * skip; i++) {
            if (fgets(str, MAX_VEGPARAM_LINE_LENGTH, vegparam) == NULL) {
                log_err("unexpected EOF for cell %i while reading root zones "
                        "and LAI", vegcel);
            }
        }
    }
    fgets(str, MAX_VEGPARAM_LINE_LENGTH, vegparam); // read newline at end of veg class line to advance to next line
    if (vegcel != gridcel) {
        log_err("Grid cell %d not found", gridcel);
    }

    // Make sure to allocate extra memory for bare soil tile
    // and optionally an above-treeline veg tile
    MaxVeg = vegetat_type_num + 1;
    if (options.AboveTreelineVeg >= 0) {
        MaxVeg++;
    }

    /** Allocate memory for vegetation grid cell parameters **/
    temp = calloc(MaxVeg, sizeof(*temp));
    Cv_sum = 0.0;

    for (i = 0; i < vegetat_type_num; i++) {
        temp[i].zone_depth = calloc(options.ROOT_ZONES,
                                    sizeof(*(temp[i].zone_depth)));
        temp[i].zone_fract = calloc(options.ROOT_ZONES,
                                    sizeof(*(temp[i].zone_fract)));
        temp[i].vegetat_type_num = vegetat_type_num;

        /* Upper boundaries of canopy layers, expressed in terms of fraction of total LAI  */
        if (options.CARBON) {
            temp[i].CanopLayerBnd = calloc(options.Ncanopy,
                                           sizeof(*(temp[i].CanopLayerBnd)));
            for (cidx = 0; cidx < options.Ncanopy; cidx++) {
                /* apportion LAI equally among layers */
                temp[i].CanopLayerBnd[cidx] =
                    (double) ((cidx + 1)) / (double) (options.Ncanopy);
            }
        }

        // Read the root zones line
        if (fgets(line, MAXSTRING, vegparam) == NULL) {
            log_err("unexpected EOF for cell %i while reading "
                    "vegetat_type_num %d", vegcel, vegetat_type_num);
        }
        strcpy(tmpline, line);
        ttrim(tmpline);
        token = strtok(tmpline, delimiters); /*  token => veg_class, move 'line' pointer to next field */
        Nfields = 0;
        vegarr[Nfields] =
            calloc(MAX_VEGPARAM_LINE_LENGTH, sizeof(*(vegarr[Nfields])));
        strcpy(vegarr[Nfields], token);
        Nfields++;

        token = strtok(NULL, delimiters);
        while (token != NULL && (length = strlen(token)) == 0) {
            token = strtok(NULL, delimiters);
        }
        while (token != NULL) {
            vegarr[Nfields] = calloc(MAX_VEGPARAM_LINE_LENGTH,
                                     sizeof(*(vegarr[Nfields])));
            strcpy(vegarr[Nfields], token);
            Nfields++;
            token = strtok(NULL, delimiters);
            while (token != NULL && (length = strlen(token)) == 0) {
                token = strtok(NULL, delimiters);
            }
        }

        NfieldsMax = 2 + 2 * options.ROOT_ZONES; /* Number of expected fields this line */
        if (options.BLOWING) {
            NfieldsMax += 3;
        }
        if (Nfields != NfieldsMax) {
            log_err("Cell %d - expecting %d fields but found %d in veg line %s",
                    gridcel, NfieldsMax, Nfields, line);
        }

        temp[i].LAKE = 0;
        temp[i].veg_class = atoi(vegarr[0]);
        temp[i].Cv = atof(vegarr[1]);
        depth_sum = 0;
        sum = 0.;
        for (j = 0; j < options.ROOT_ZONES; j++) {
            temp[i].zone_depth[j] = atof(vegarr[2 + j * 2]);
            temp[i].zone_fract[j] = atof(vegarr[3 + j * 2]);
            depth_sum += temp[i].zone_depth[j];
            sum += temp[i].zone_fract[j];
        }
        if (depth_sum <= 0) {
            log_err("Root zone depths must sum to a value greater than 0.");
        }
        if (sum != 1.) {
            log_warn("Root zone fractions sum to more than 1 ( = %f), "
                     "normalizing fractions.  If the sum is large, check that "
                     "your vegetation parameter file is in the form - <zone 1 "
                     "depth> <zone 1 fract> <zone 2 depth> <zone 2 fract>...",
                     sum);
            for (j = 0; j < options.ROOT_ZONES; j++) {
                temp[i].zone_fract[j] /= sum;
            }
        }

        if (options.BLOWING) {
            j = 2 * options.ROOT_ZONES;
            temp[i].sigma_slope = atof(vegarr[2 + j]);
            temp[i].lag_one = atof(vegarr[3 + j]);
            temp[i].fetch = atof(vegarr[4 + j]);
            if (temp[i].sigma_slope <= 0. || temp[i].lag_one <= 0.) {
                log_err("Deviation of terrain slope must be greater than 0.");
            }
            if (temp[i].fetch < 1.0) {
                log_err("BLOWING parameter fetch should be >> 1 but "
                        "cell %i has fetch = %.2f", gridcel, temp[i].fetch);
            }
        }

        veg_class = MISSING;
        for (j = 0; j < Nveg_type; j++) {
            if (temp[i].veg_class == veg_lib[j].veg_class) {
                veg_class = j;
            }
        }
        if (veg_class == MISSING) {
            log_err("The vegetation class id %i in vegetation tile %i from "
                    "cell %i is not defined in the vegetation library file.",
                    temp[i].veg_class, i, gridcel);
        }
        else {
            temp[i].veg_class = veg_class;
        }

        Cv_sum += temp[i].Cv;

        for (k = 0; k < Nfields; k++) {
            free(vegarr[k]);
        }

        for (j = 0; j < MONTHS_PER_YEAR; j++) {
            temp[i].albedo[j] = veg_lib[temp[i].veg_class].albedo[j];
            temp[i].displacement[j] =
                veg_lib[temp[i].veg_class].displacement[j];
            temp[i].fcanopy[j] = veg_lib[temp[i].veg_class].fcanopy[j];
            temp[i].LAI[j] = veg_lib[temp[i].veg_class].LAI[j];
            temp[i].roughness[j] = veg_lib[temp[i].veg_class].roughness[j];
            temp[i].Wdmax[j] = veg_lib[temp[i].veg_class].Wdmax[j];
        }

        if (options.VEGPARAM_LAI) {
            // Read the LAI line
            if (fgets(line, MAXSTRING, vegparam) == NULL) {
                log_err("Unexpected EOF for cell %i while reading LAI for "
                        "vegetat_type_num %d", vegcel, vegetat_type_num);
            }
            Nfields = 0;
            vegarr[Nfields] = calloc(MAX_VEGPARAM_LINE_LENGTH,
                                     sizeof(*(vegarr[Nfields])));
            strcpy(tmpline, line);
            ttrim(tmpline);
            token = strtok(tmpline, delimiters);
            strcpy(vegarr[Nfields], token);
            Nfields++;

            while ((token = strtok(NULL, delimiters)) != NULL) {
                vegarr[Nfields] = calloc(MAX_VEGPARAM_LINE_LENGTH,
                                         sizeof(*(vegarr[Nfields])));
                strcpy(vegarr[Nfields], token);
                Nfields++;
            }
            NfieldsMax = MONTHS_PER_YEAR; /* For LAI */
            if (Nfields != NfieldsMax) {
                log_err("cell %d - expecting %d LAI values but found "
                        "%d in line %s", gridcel, NfieldsMax, Nfields, line);
            }

            if (options.LAI_SRC == FROM_VEGPARAM) {
                for (j = 0; j < MONTHS_PER_YEAR; j++) {
                    tmp = atof(vegarr[j]);
                    if (tmp != NODATA_VH) {
                        temp[i].LAI[j] = tmp;
                    }
                    if (veg_lib[temp[i].veg_class].overstory &&
                        temp[i].LAI[j] == 0) {
                        log_err("cell %d, veg tile %d: the specified "
                                "veg class (%d) is listed as an overstory "
                                "class in the veg LIBRARY, but the LAI given "
                                "in the veg PARAM FILE for this tile for "
                                "month %zu is 0.", gridcel, i + 1,
                                temp[i].veg_class + 1, j + 1);
                    }
                    temp[i].Wdmax[j] =
                        param.VEG_LAI_WATER_FACTOR *
                        temp[i].LAI[j];
                }
            }
            for (k = 0; k < Nfields; k++) {
                free(vegarr[k]);
            }
        }

        if (options.VEGPARAM_FCAN) {
            // Read the fcanopy line
            if (fgets(line, MAXSTRING, vegparam) == NULL) {
                log_err("unexpected EOF for cell %i while reading fcanopy "
                        "for vegetat_type_num %d", vegcel, vegetat_type_num);
            }
            Nfields = 0;
            vegarr[Nfields] = calloc(MAX_VEGPARAM_LINE_LENGTH,
                                     sizeof(*(vegarr[Nfields])));
            strcpy(tmpline, line);
            ttrim(tmpline);
            token = strtok(tmpline, delimiters);
            strcpy(vegarr[Nfields], token);
            Nfields++;

            while ((token = strtok(NULL, delimiters)) != NULL) {
                vegarr[Nfields] = calloc(MAX_VEGPARAM_LINE_LENGTH,
                                         sizeof(*(vegarr[Nfields])));
                strcpy(vegarr[Nfields], token);
                Nfields++;
            }
            NfieldsMax = MONTHS_PER_YEAR; /* For fcanopy */
            if (Nfields != NfieldsMax) {
                log_err("cell %d - expecting %d fcanopy values but found %d "
                        "in line %s", gridcel, NfieldsMax, Nfields, line);
            }

            if (options.FCAN_SRC == FROM_VEGPARAM) {
                for (j = 0; j < MONTHS_PER_YEAR; j++) {
                    tmp = atof(vegarr[j]);
                    if (tmp != NODATA_VH) {
                        temp[i].fcanopy[j] = tmp;
                    }
                }
            }
            for (k = 0; k < Nfields; k++) {
                free(vegarr[k]);
            }
        }

        if (options.VEGPARAM_ALB) {
            // Read the albedo line
            if (fgets(line, MAXSTRING, vegparam) == NULL) {
                log_err("unexpected EOF for cell %i while reading albedo for "
                        "vegetat_type_num %d", vegcel, vegetat_type_num);
            }
            Nfields = 0;
            vegarr[Nfields] = calloc(MAX_VEGPARAM_LINE_LENGTH,
                                     sizeof(*(vegarr[Nfields])));
            strcpy(tmpline, line);
            ttrim(tmpline);
            token = strtok(tmpline, delimiters);
            strcpy(vegarr[Nfields], token);
            Nfields++;

            while ((token = strtok(NULL, delimiters)) != NULL) {
                vegarr[Nfields] = calloc(MAX_VEGPARAM_LINE_LENGTH,
                                         sizeof(*(vegarr[Nfields])));
                strcpy(vegarr[Nfields], token);
                Nfields++;
            }
            NfieldsMax = MONTHS_PER_YEAR; /* For albedo */
            if (Nfields != NfieldsMax) {
                log_err("cell %d - expecting %d albedo values but found %d in "
                        "line %s", gridcel, NfieldsMax, Nfields, line);
            }

            if (options.ALB_SRC == FROM_VEGPARAM) {
                for (j = 0; j < MONTHS_PER_YEAR; j++) {
                    tmp = atof(vegarr[j]);
                    if (tmp != NODATA_VH) {
                        temp[i].albedo[j] = tmp;
                    }
                }
            }
            for (k = 0; k < Nfields; k++) {
                free(vegarr[k]);
            }
        }

        // Determine if cell contains non-overstory vegetation
        if (options.COMPUTE_TREELINE && !veg_lib[temp[i].veg_class].overstory) {
            NoOverstory++;
        }
    }

    // Determine if we have bare soil
    if (Cv_sum > 1.0) {
        log_warn("Cv_sum exceeds 1.0 (%f) at grid cell %d, fractions being "
                 "adjusted to equal 1", Cv_sum, gridcel);
        for (j = 0; j < (size_t)vegetat_type_num; j++) {
            temp[j].Cv = temp[j].Cv / Cv_sum;
        }
        Cv_sum = 1.;
    }
    else if (Cv_sum > 0.99 && Cv_sum < 1.0) {
        log_warn("Cv > 0.99 and Cv < 1.0 at grid cell %d, model "
                 "assuming that bare soil is not to be run - fractions being "
                 "adjusted to equal 1",
                 gridcel);
        for (j = 0; j < (size_t)vegetat_type_num; j++) {
            temp[j].Cv = temp[j].Cv / Cv_sum;
        }
        Cv_sum = 1.;
    }

    // Handle veg above the treeline
    if (options.SNOW_BAND > 1 && options.COMPUTE_TREELINE &&
        (!NoOverstory && Cv_sum == 1.)) {
        // All vegetation in the current cell is defined with overstory.
        // Add default non-overstory vegetation so that snow bands above treeline
        // can be sucessfully simulated.

        if (options.AboveTreelineVeg < 0) {
            // Above treeline snowband should be treated as bare soil
            for (j = 0; j < (size_t)vegetat_type_num; j++) {
                temp[j].Cv -= (0.001 / (double) vegetat_type_num);
            }
            Cv_sum -= 0.001;
        }
        else {
            // Above treeline snowband should use the defined vegetation
            // add vegetation to typenum
            // check that veg type exists in library and does not have overstory
            if (vegetat_type_num > 0) {
                for (j = 0; j < (size_t)vegetat_type_num; j++) {
                    temp[j].Cv -= (0.001 / (double) vegetat_type_num);
                    temp[j].vegetat_type_num++;
                }

                temp[vegetat_type_num].Cv = 0.001;
                temp[vegetat_type_num].veg_class = options.AboveTreelineVeg;
                temp[vegetat_type_num].zone_depth = calloc(options.ROOT_ZONES,
                                                           sizeof(double));
                temp[vegetat_type_num].zone_fract = calloc(options.ROOT_ZONES,
                                                           sizeof(double));
                temp[vegetat_type_num].vegetat_type_num = vegetat_type_num + 1;

                // Since root zones are not defined they are copied from the last
                // vegetation type.
                for (j = 0; j < options.ROOT_ZONES; j++) {
                    temp[vegetat_type_num].zone_depth[j] =
                        temp[vegetat_type_num - 1].zone_depth[j];
                    temp[vegetat_type_num].zone_fract[j] =
                        temp[vegetat_type_num - 1].zone_fract[j];
                }
            }

            // Identify current vegetation class
            veg_class = MISSING;
            for (j = 0; j < Nveg_type; j++) {
                if (temp[vegetat_type_num].veg_class == veg_lib[j].veg_class) {
                    veg_class = j;
                    break;
                }
            }
            if (veg_class == MISSING) {
                log_err("The vegetation class id %i defined for "
                        "above-treeline from cell %i is not defined in the "
                        "vegetation library file.",
                        temp[vegetat_type_num].veg_class, gridcel);
            }
            else {
                temp[vegetat_type_num].veg_class = veg_class;
            }

            if (veg_lib[veg_class].overstory) {
                log_err("Vegetation class %i is defined to have overstory, so "
                        "it cannot be used as the default vegetation type for "
                        "above canopy snow bands.",
                        veg_lib[veg_class].veg_class);
            }
        }
        vegetat_type_num = temp[0].vegetat_type_num;
    }

    // Default bare soil tile - not specified in vegparam file
    i = vegetat_type_num;
    temp[i].veg_class = Nveg_type; // Create a veg_class ID for bare soil, which is not mentioned in the veg library
    temp[i].Cv = 1.0 - Cv_sum;
    if (temp[i].Cv < 0) {
        temp[i].Cv = 0;
    }
    // Don't allocate any root-zone-related arrays
    if (options.BLOWING) {
        if (vegetat_type_num > 0) {
            temp[i].sigma_slope = temp[0].sigma_slope;
            temp[i].lag_one = temp[0].lag_one;
            temp[i].fetch = temp[0].fetch;
        }
        else {
            temp[i].sigma_slope = 0.005;
            temp[i].lag_one = 0.95;
            temp[i].fetch = 2000;
        }
    }
    for (j = 0; j < MONTHS_PER_YEAR; j++) {
        temp[i].albedo[j] = veg_lib[temp[i].veg_class].albedo[j];
        temp[i].displacement[j] =
            veg_lib[temp[i].veg_class].displacement[j];
        temp[i].fcanopy[j] = veg_lib[temp[i].veg_class].fcanopy[j];
        temp[i].LAI[j] = veg_lib[temp[i].veg_class].LAI[j];
        temp[i].roughness[j] = veg_lib[temp[i].veg_class].roughness[j];
        temp[i].Wdmax[j] = veg_lib[temp[i].veg_class].Wdmax[j];
    }

    return temp;
}

/* trim trailing newlines */

#define END '\0'
#define NEW '\n'

/******************************************************************************
 * @brief    trim trailing newlines
 *****************************************************************************/
void
ttrim(char *c)
{
    while ((*c++ != END)) {
        ;
    }
    for (--c; *--c == NEW; *c = END) {
        ;
    }
}
