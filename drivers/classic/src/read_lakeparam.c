/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine reads in lake parameters for the current grid cell.
 *
 * It will either calculate the lake area v. depth profile from a parabolic
 * curve or read in constant values depending on the LAKE_PROFILE flag.
 *
 * @section LICENSE
 *
 * The Variable Infiltration Capacity (VIC) macroscale hydrological model
 * Copyright (C) 2014 The Land Surface Hydrology Group, Department of Civil
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

#include <vic_def.h>
#include <vic_run.h>
#include <vic_driver_classic.h>

/******************************************************************************
 * @brief    This routine reads in lake parameters for the current grid cell.
 *****************************************************************************/
lake_con_struct
read_lakeparam(FILE           *lakeparam,
               soil_con_struct soil_con,
               veg_con_struct *veg_con)
{
    extern option_struct     options;
    extern parameters_struct param;

    size_t                   i;
    unsigned                 lakecel;
    double                   tempdz;
    double                   radius, x;
    char                     tmpstr[MAXSTRING + 1];
    int                      ErrFlag;

    lake_con_struct          temp;

    /*******************************************************************/
    /* Read in general lake parameters.                           */
    /******************************************************************/

    fscanf(lakeparam, "%u %d", &lakecel, &temp.lake_idx);
    while (lakecel != soil_con.gridcel && !feof(lakeparam)) {
        fgets(tmpstr, MAXSTRING, lakeparam); // grid cell number, etc.
        if (temp.lake_idx >= 0) {
            fgets(tmpstr, MAXSTRING, lakeparam); // lake depth-area relationship
        }
        fscanf(lakeparam, "%u %d", &lakecel, &temp.lake_idx);
    }

    // cell number not found
    if (feof(lakeparam)) {
        log_err("Unable to find cell %d in the lake parameter file",
                soil_con.gridcel);
    }

    // read lake parameters from file
    if (temp.lake_idx >= 0) {
        veg_con[temp.lake_idx].LAKE = 1;
        fscanf(lakeparam, "%zu", &temp.numnod);
        if (temp.numnod < 1) {
            log_err("Number of vertical lake nodes (%zu) for cell %d specified "
                    "in the lake parameter file is < 1; increase this number "
                    "to at least 1.", temp.numnod, soil_con.gridcel);
        }
        if (temp.numnod > MAX_LAKE_NODES) {
            log_err("Number of lake nodes (%zu) in cell %d specified in the "
                    "lake parameter file exceeds the maximum allowable (%d), "
                    "edit MAX_LAKE_NODES in user_def.h.", temp.numnod,
                    soil_con.gridcel, MAX_LAKE_NODES);
        }
        fscanf(lakeparam, "%lf", &temp.mindepth);
        if (temp.mindepth < 0) {
            log_err("Minimum lake depth (%f) for cell %d specified in the "
                    "lake parameter file is < 0; increase this number to at "
                    "least 0.", temp.mindepth, soil_con.gridcel);
        }
        fscanf(lakeparam, "%lf", &temp.wfrac);
        if (temp.wfrac < 0 || temp.wfrac > 1) {
            log_err("Lake outlet width fraction (%f) for cell %d specified in "
                    "the lake parameter file falls outside the range 0 to 1.  "
                    "Change wfrac to be between 0 and 1.", temp.wfrac,
                    soil_con.gridcel);
        }
        fscanf(lakeparam, "%lf", &temp.depth_in);
        if (temp.depth_in < 0) {
            log_err("Initial lake depth (%f) for cell %d specified in the "
                    "lake parameter file is < 0; increase this number to at "
                    "least 1.", temp.depth_in, soil_con.gridcel);
        }
        fscanf(lakeparam, "%lf", &temp.rpercent);
        if (temp.rpercent < 0 || temp.rpercent > 1) {
            log_err("Fraction of runoff entering lake catchment (%f) for cell "
                    "%d specified in the lake parameter file falls outside the"
                    " range 0 to 1.  Change rpercent to be between 0 and 1.",
                    temp.rpercent, soil_con.gridcel);
        }
    }
    else { // no lake exists anywhere in this grid cell
        temp.numnod = 0;
        temp.mindepth = 0;
        temp.maxdepth = 0;
        temp.Cl[0] = 0;
        temp.basin[0] = 0;
        temp.z[0] = 0;
        temp.minvolume = 0;
        temp.maxvolume = 0;
        temp.wfrac = 0;
        temp.depth_in = 0;
        temp.rpercent = 0;
        temp.bpercent = 0;
        fgets(tmpstr, MAXSTRING, lakeparam); // skip to end of line
        return temp;
    }
    temp.bpercent = temp.rpercent;

    /*******************************************************************/
    /* Find lake basin area with depth.                           */
    /******************************************************************/

    /* Read in parameters to calculate lake profile. */
    if (!options.LAKE_PROFILE) {
        log_info("LAKE PROFILE being computed.");

        fscanf(lakeparam, "%lf %lf", &temp.z[0], &temp.Cl[0]);
        temp.maxdepth = temp.z[0];
        tempdz = (temp.maxdepth) / ((double)temp.numnod);
        if (temp.Cl[0] < 0.0 || temp.Cl[0] > 1.0) {
            log_err("Lake area fraction (%f) for cell (%d) specified in the "
                    "lake parameter file must be a fraction between 0 and 1.",
                    temp.Cl[0], soil_con.gridcel);
        }

        fgets(tmpstr, MAXSTRING, lakeparam);

        temp.basin[0] = temp.Cl[0] * soil_con.cell_area;

        /**********************************************
           Compute depth area relationship.
        **********************************************/

        radius = sqrt(temp.basin[0] / CONST_PI);

        temp.maxvolume = 0.0;
        for (i = 1; i <= temp.numnod; i++) {
            temp.z[i] = (temp.numnod - i) * tempdz;
            if (temp.z[i] < 0.0) {
                temp.z[i] = 0.0;
            }
            x = pow(temp.z[i] / temp.maxdepth, param.LAKE_BETA) * radius;
            temp.basin[i] = CONST_PI * x * x;
            temp.maxvolume += (temp.basin[i] + temp.basin[i - 1]) * tempdz / 2.;
        }
    }
    /* Read in basin area for each layer depth. */
    /* Assumes that the lake bottom area is zero. */

    else {
        log_info("Reading in the specified lake profile.");
        temp.maxvolume = 0.0;
        temp.Cl[0] = 0; // initialize to 0 in case no lake is defined
        for (i = 0; i < temp.numnod; i++) {
            fscanf(lakeparam, "%lf %lf", &temp.z[i], &temp.Cl[i]);
            temp.basin[i] = temp.Cl[i] * soil_con.cell_area;

            if (i == 0) {
                temp.maxdepth = temp.z[i];
                tempdz = (temp.maxdepth) / ((double) temp.numnod);
            }

            if (temp.Cl[0] < 0.0 || temp.Cl[0] > 1.0) {
                log_err("Lake area fraction (%f) for cell (%d) specified in "
                        "the lake parameter file must be a fraction between 0 "
                        "and 1.", temp.Cl[0], soil_con.gridcel);
            }
        }
        temp.z[temp.numnod] = 0.0;
        temp.basin[temp.numnod] = 0.0;
        temp.Cl[temp.numnod] = 0.0;

        for (i = 1; i <= temp.numnod; i++) {
            temp.maxvolume +=
                (temp.basin[i] +
                 temp.basin[i - 1]) * (temp.z[i - 1] - temp.z[i]) / 2.;
        }
    }

    // Compute volume corresponding to mindepth
    ErrFlag = get_volume(temp, temp.mindepth, &(temp.minvolume));
    if (ErrFlag == ERROR) {
        log_err("problem in get_volume(): depth %f volume %f rec %d",
                temp.mindepth, temp.minvolume, 0);
    }

    // Make sure min < max
    if (temp.mindepth > temp.maxdepth) {
        log_err("Adjusted minimum lake depth %f exceeds the adjusted maximum "
                "lake depth %f.", temp.mindepth, temp.maxdepth);
    }

    // Validate initial conditions
    if (temp.depth_in > temp.maxdepth) {
        log_warn("Initial lake depth %f exceeds the maximum lake depth %f; "
                 "setting initial lake depth equal to max lake depth.",
                 temp.depth_in, temp.maxdepth);
        temp.depth_in = temp.maxdepth;
    }
    else if (temp.depth_in < 0) {
        log_warn("Initial lake depth %f < 0; setting to 0.", temp.depth_in);
        temp.depth_in = 0;
    }

    log_info("Lake plus wetland area = %e km2", temp.basin[0] /
             (M_PER_KM * M_PER_KM));
    return temp;
}
