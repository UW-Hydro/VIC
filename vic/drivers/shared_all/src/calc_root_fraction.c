/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine computes the fraction of roots in each soil layer based on the
 * root zone distribution defined in the vegetation parameter file.  Roots are
 * assumed to be linearly distributed within each root zone.
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

#include <vic_driver_shared_all.h>

/******************************************************************************
 * @brief    This routine computes the fraction of roots in each soil layer.
 *
 *           Overview: take the union of the set of soil layer boundaries and
 *           the set of root zone boundaries, and define depth intervals 'dD'
 *           that are the spaces between these boundaries. For each soil layer,
 *           the total root fraction is the sum of root fractions in all of the
 *           dD[i] that fall within the layer.  To accomplish this sum
 *           efficiently, first compute root density (fraction per unit depth)
 *           as:
 *             root_dens[zone] = fract[zone] / depth[zone]
 *           Assign root densities to each dD[i] as:
 *             root_dens[i] = root_dens[zone] for the zone that dD[i] is in
 *           Then, compute each layer's total root fraction as:
 *             root[layer] = sum_over_i(dD[i] * root_dens[i])
 *****************************************************************************/
void
calc_root_fractions(veg_con_struct  *veg_con,
                    soil_con_struct *soil_con)
{
    extern option_struct options;

    int                  Nveg;
    int                  veg;
    size_t               layer;
    size_t               zone;
    size_t               ltmp;
    double               Lsum;
    double               Lsumprev;
    double               Zsum;
    double               Dsum;
    double               Dsumprev;
    double               dD;
    double               sum_fract;
    double              *root_dens;
    double               dum;

    root_dens = calloc(options.ROOT_ZONES, sizeof(double));

    Nveg = veg_con[0].vegetat_type_num;

    for (veg = 0; veg < Nveg; veg++) {

        // Check that root fractions sum to 1.0
        dum = 0.0;
        for (zone = 0; zone < options.ROOT_ZONES; zone++) {
            dum += veg_con[veg].zone_fract[zone];
        }
        if (!assert_close_double(dum, 1, 0, 1e-4)) {
            log_err("Input root fractions do not sum to 1.0: %f, "
                    "veg class: %d", dum, veg_con[veg].veg_class);
        }

        // Compute density of root fractions in each root zone
        for (zone = 0; zone < options.ROOT_ZONES; zone++) {
            if (veg_con[veg].zone_depth[zone] > 0) {
                root_dens[zone] = veg_con[veg].zone_fract[zone] /
                                  veg_con[veg].zone_depth[zone];
            }
            else {
                root_dens[zone] = 1.0;
            }
        }

        layer = 0;
        zone = 0;
        Lsum = soil_con->depth[layer];
        Lsumprev = 0;
        Zsum = veg_con[veg].zone_depth[zone];
        Dsum = 0;
        Dsumprev = 0;
        sum_fract = 0;

        while (layer < options.Nlayer || zone < options.ROOT_ZONES) {

            // Determine the depth interval for consideration
            // Depth intervals span from the previous layer or zone boundary
            // to the next layer or zone boundary.
            if (Lsum <= Zsum && layer < options.Nlayer) {
                Dsum = Lsum;
            }
            else {
                Dsum = Zsum;
            }
            dD = Dsum - Dsumprev;
            Dsumprev = Dsum;

            // Add root fraction in this interval to the running total
            // for this layer
            if (zone < options.ROOT_ZONES) {
                sum_fract += dD * root_dens[zone];
            }

            // Assign the total root fraction for this soil layer
            // Wait to do this until we've completed a soil layer
            // and either we're not at the final soil layer,
            // or we're at the final layer and we've completed all root zones
            if ( Lsum > Lsumprev && Dsum == Lsum &&
                 ( layer < options.Nlayer - 1 ||
                   ( zone >= options.ROOT_ZONES - 1 && Dsum == Zsum ) ) ) {

                // Assign the total root fraction
                ltmp = layer;
                if (layer >= options.Nlayer - 1) {
                    ltmp = options.Nlayer - 1;
                }
                veg_con[veg].root[ltmp] = sum_fract;

                // Reset running total for next layer
                sum_fract = 0;

            }

            // Decide whether to increment layer or zone
            if (Lsum < Zsum) {
                // zone extends beyond layer, so increment layer
                // it's ok to exceed limit; this tells algorithm we're done
                layer++;
                if (layer < options.Nlayer) {
                    Lsumprev = Lsum;
                    Lsum += soil_con->depth[layer];
                }
                else {
                    // We have reached the bottom of the soil column;
                    // extend Lsum to include all remaining roots within
                    // the final soil layer
                    Lsum = Zsum;
                }
            }
            else if (Lsum > Zsum) {
                // layer extends beyond zone, so increment zone
                // it's ok to exceed limit; this tells algorithm we're done
                zone++;
                if (zone < options.ROOT_ZONES) {
                    Zsum += veg_con[veg].zone_depth[zone];
                }
                else {
                    // We have reached the bottom of the root zone;
                    // extend Zsum to include the remainder of the soil layer
                    Zsum = Lsum;
                }
            }
            else {
                // layer and zone end at the same depth, so increment both
                // it's ok to exceed limit; this tells algorithm we're done
                layer++;
                if (layer < options.Nlayer) {
                    Lsumprev = Lsum;
                    Lsum += soil_con->depth[layer];
                }
                zone++;
                if (zone < options.ROOT_ZONES) {
                    Zsum += veg_con[veg].zone_depth[zone];
                }
            }

        }

        // Final check on root fractions. If they don't sum to 1, throw error
        // Otherwise, rescale by sum to eliminate small rounding errors
        dum = 0.;
        for (layer = 0; layer < options.Nlayer; layer++) {
            if (veg_con[veg].root[layer] < 1.e-4) {
                veg_con[veg].root[layer] = 0.;
            }
            dum += veg_con[veg].root[layer];
        }
        if (!assert_close_double(dum, 1, 0, 1e-4)) {
            log_err("Soil layer root fractions do not sum to 1.0: %f, "
                    "veg class: %d", dum, veg_con[veg].veg_class);
        }
        else {
            for (layer = 0; layer < options.Nlayer; layer++) {
                veg_con[veg].root[layer] /= dum;
            }
        }
    }

    free((char *) root_dens);
}
