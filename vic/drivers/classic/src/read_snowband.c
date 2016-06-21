/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine reads snow band median elevaton, and precipitation fraction for
 * use with the snow model.
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
 * @brief    This routine reads snow band median elevaton, and precipitation
 *           fraction for use with the snow model.
 *****************************************************************************/
void
read_snowband(FILE            *snowband,
              soil_con_struct *soil_con)
{
    extern option_struct     options;
    extern parameters_struct param;

    char                     ErrStr[MAXSTRING];
    size_t                   band;
    size_t                   Nbands;
    unsigned int             cell;
    double                   total;
    double                   area_fract;
    double                   prec_frac;
    double                   band_elev;
    double                   avg_elev;

    Nbands = options.SNOW_BAND;

    if (Nbands > 1) {
        /** Find Current Grid Cell in SnowBand File **/
        fscanf(snowband, "%d", &cell);
        while (cell != soil_con->gridcel && !feof(snowband)) {
            fgets(ErrStr, MAXSTRING, snowband);
            fscanf(snowband, "%d", &cell);
        }

        if (feof(snowband)) {
            log_warn("Cannot find current gridcell (%i) in snow band file; "
                     "setting cell to have one elevation band.",
                     soil_con->gridcel);
            /** 1 band is the default; no action necessary **/
            return;
        }

        /** Read Area Fraction **/
        total = 0.;
        for (band = 0; band < Nbands; band++) {
            fscanf(snowband, "%lf", &area_fract);
            if (area_fract < 0) {
                log_err("Negative snow band area fraction (%f) read from file",
                        area_fract);
            }
            soil_con->AreaFract[band] = area_fract;
            total += area_fract;
        }
        if (total != 1.) {
            log_warn("Sum of the snow band area fractions does not equal "
                     "1 (%f), dividing each fraction by the sum", total);
            for (band = 0; band < options.SNOW_BAND; band++) {
                soil_con->AreaFract[band] /= total;
            }
        }

        /** Read Band Elevation **/
        avg_elev = 0;
        for (band = 0; band < Nbands; band++) {
            fscanf(snowband, "%lf", &band_elev);
            if (band_elev < 0) {
                log_err("Negative snow band elevation (%f) read from file",
                        band_elev);
            }
            soil_con->BandElev[band] = band_elev;
            avg_elev += soil_con->BandElev[band] * soil_con->AreaFract[band];
        }
        if (fabs(avg_elev - soil_con->elevation) > 1.0) {
            log_warn("average band elevation %f not equal to grid_cell "
                     "average elevation %f; setting grid cell elevation to "
                     "average band elevation.", avg_elev,
                     soil_con->elevation);
            soil_con->elevation = (double) avg_elev;
        }
        for (band = 0; band < Nbands; band++) {
            soil_con->Tfactor[band] =
                (soil_con->BandElev[band] -
                 soil_con->elevation) * param.LAPSE_RATE;
        }

        /** Read Precipitation Fraction **/
        total = 0.;
        for (band = 0; band < options.SNOW_BAND; band++) {
            fscanf(snowband, "%lf", &prec_frac);
            if (prec_frac < 0) {
                log_err("Snow band precipitation fraction (%f) must be "
                        "between 0 and 1", prec_frac);
            }
            if (prec_frac > 0 && soil_con->AreaFract[band] == 0) {
                log_err("Snow band precipitation fraction (%f) should be 0 "
                        "when the area fraction is 0. (band = %zu)",
                        prec_frac, band);
            }
            soil_con->Pfactor[band] = prec_frac;
            total += prec_frac;
        }
        if (total != 1.) {
            log_warn("Sum of the snow band precipitation fractions "
                     "does not equal %d (%f), dividing each fraction by the "
                     "sum", 1, total);
            for (band = 0; band < options.SNOW_BAND; band++) {
                soil_con->Pfactor[band] /= total;
            }
        }
        for (band = 0; band < options.SNOW_BAND; band++) {
            if (soil_con->AreaFract[band] > 0) {
                soil_con->Pfactor[band] /= soil_con->AreaFract[band];
            }
            else {
                soil_con->Pfactor[band] = 0.;
            }
        }
    }
}
