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
    size_t               i;
    size_t               n_iter;
    double               sum_fract;
    double               dum;
    double               Zstep;
    double               Zsum;
    double               Lstep;
    double               Lsum;
    double               Zmin_fract;
    double               Zmin_depth;
    double               Zmax;

    Nveg = veg_con[0].vegetat_type_num;

    for (veg = 0; veg < Nveg; veg++) {
        sum_fract = 0;
        layer = 0;
        Lstep = soil_con->depth[layer];
        Lsum = Lstep;
        Zsum = 0;
        zone = 0;

        n_iter = 0;
        while (zone < options.ROOT_ZONES) {
            n_iter++;
            if (n_iter > MAX_ROOT_ITER) {
                log_warn("veg=%d of Nveg=%d", veg, Nveg);
                log_warn("zone %zu of %zu ROOT_ZONES", zone,
                         options.ROOT_ZONES);
                log_err("stuck in an infinite loop");
            }

            Zstep = veg_con[veg].zone_depth[zone];
            if ((Zsum + Zstep) <= Lsum && Zsum >= Lsum - Lstep) {
                /** CASE 1: Root Zone Completely in Soil Layer **/
                sum_fract += veg_con[veg].zone_fract[zone];
            }
            else {
                /** CASE 2: Root Zone Partially in Soil Layer **/
                if (Zsum < Lsum - Lstep) {
                    /** Root zone starts in previous soil layer **/
                    Zmin_depth = Lsum - Lstep;
                    Zmin_fract = linear_interp(Zmin_depth, Zsum, Zsum + Zstep,
                                               0,
                                               veg_con[veg].zone_fract[zone]);
                }
                else {
                    /** Root zone starts in current soil layer **/
                    Zmin_depth = Zsum;
                    Zmin_fract = 0.;
                }
                if (Zsum + Zstep <= Lsum) {
                    /** Root zone ends in current layer **/
                    Zmax = Zsum + Zstep;
                }
                else {
                    /** Root zone extends beyond bottom of current layer **/
                    Zmax = Lsum;
                }
                sum_fract += linear_interp(Zmax, Zsum, Zsum + Zstep, 0,
                                           veg_con[veg].zone_fract[zone]) -
                             Zmin_fract;
            }

            /** Update Root Zone and Soil Layer **/
            if (Zsum + Zstep < Lsum) {
                Zsum += Zstep;
                zone++;
            }
            else if (Zsum + Zstep == Lsum) {
                Zsum += Zstep;
                zone++;
                if (layer < options.Nlayer) {
                    veg_con[veg].root[layer] = sum_fract;
                    sum_fract = 0.;
                }
                layer++;
                if (layer < options.Nlayer) {
                    Lstep = soil_con->depth[layer];
                    Lsum += Lstep;
                }
                else if (layer == options.Nlayer && zone < options.ROOT_ZONES) {
                    Zstep = (double) veg_con[veg].zone_depth[zone];
                    Lstep = Zsum + Zstep - Lsum;
                    if (zone < options.ROOT_ZONES - 1) {
                        for (i = zone + 1; i < options.ROOT_ZONES; i++) {
                            Lstep += veg_con[veg].zone_depth[i];
                        }
                    }
                    Lsum += Lstep;
                }
            }
            else if (Zsum + Zstep > Lsum) {
                zone++;
                if (layer < options.Nlayer) {
                    veg_con[veg].root[layer] = sum_fract;
                    sum_fract = 0.;
                }
                layer++;
                if (layer < options.Nlayer) {
                    Lstep = soil_con->depth[layer];
                    Lsum += Lstep;
                }
                else if (layer == options.Nlayer) {
                    Lstep = Zsum + Zstep - Lsum;
                    if (zone < options.ROOT_ZONES - 1) {
                        for (i = zone + 1; i < options.ROOT_ZONES; i++) {
                            Lstep += veg_con[veg].zone_depth[i];
                        }
                    }
                    Lsum += Lstep;
                }
            }
        }

        if (sum_fract > 0 && layer >= options.Nlayer) {
            veg_con[veg].root[options.Nlayer - 1] += sum_fract;
        }
        else if (sum_fract > 0) {
            veg_con[veg].root[layer] += sum_fract;
        }

        dum = 0.;
        for (layer = 0; layer < options.Nlayer; layer++) {
            if (veg_con[veg].root[layer] < 1.e-4) {
                veg_con[veg].root[layer] = 0.;
            }
            dum += veg_con[veg].root[layer];
        }
        if (dum == 0.0) {
            log_err("Root fractions sum equals zero: %f , Vege Class: %d",
                    dum, veg_con[veg].veg_class);
        }
        for (layer = 0; layer < options.Nlayer; layer++) {
            veg_con[veg].root[layer] /= dum;
        }
    }
}
