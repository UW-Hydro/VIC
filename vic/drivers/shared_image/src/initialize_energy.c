/******************************************************************************
 * @section DESCRIPTION
 *
 * Initialize energy structure.
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

#include <vic_driver_shared_image.h>

/******************************************************************************
 * @brief    Initialize energy structure.
 *****************************************************************************/
void
initialize_energy(energy_bal_struct **energy,
                  size_t              Nveg)
{
    extern option_struct options;

    size_t               band;
    size_t               index;
    size_t               veg;

    // initialize miscellaneous energy balance terms
    for (veg = 0; veg <= Nveg; veg++) {
        for (band = 0; band < options.SNOW_BAND; band++) {
            /* Set fluxes to 0 */
            energy[veg][band].advected_sensible = 0.0;
            energy[veg][band].advection = 0.0;
            energy[veg][band].AtmosError = 0.0;
            energy[veg][band].AtmosLatent = 0.0;
            energy[veg][band].AtmosLatentSub = 0.0;
            energy[veg][band].AtmosSensible = 0.0;
            energy[veg][band].canopy_advection = 0.0;
            energy[veg][band].canopy_latent = 0.0;
            energy[veg][band].canopy_latent_sub = 0.0;
            energy[veg][band].canopy_refreeze = 0.0;
            energy[veg][band].canopy_sensible = 0.0;
            energy[veg][band].deltaCC = 0.0;
            energy[veg][band].deltaH = 0.0;
            energy[veg][band].error = 0.0;
            energy[veg][band].fusion = 0.0;
            energy[veg][band].grnd_flux = 0.0;
            energy[veg][band].latent = 0.0;
            energy[veg][band].latent_sub = 0.0;
            energy[veg][band].longwave = 0.0;
            energy[veg][band].LongOverIn = 0.0;
            energy[veg][band].LongUnderIn = 0.0;
            energy[veg][band].LongUnderOut = 0.0;
            energy[veg][band].melt_energy = 0.0;
            energy[veg][band].NetLongAtmos = 0.0;
            energy[veg][band].NetLongOver = 0.0;
            energy[veg][band].NetLongUnder = 0.0;
            energy[veg][band].NetShortAtmos = 0.0;
            energy[veg][band].NetShortGrnd = 0.0;
            energy[veg][band].NetShortOver = 0.0;
            energy[veg][band].NetShortUnder = 0.0;
            energy[veg][band].out_long_canopy = 0.0;
            energy[veg][band].out_long_surface = 0.0;
            energy[veg][band].refreeze_energy = 0.0;
            energy[veg][band].sensible = 0.0;
            energy[veg][band].shortwave = 0.0;
            energy[veg][band].ShortOverIn = 0.0;
            energy[veg][band].ShortUnderIn = 0.0;
            energy[veg][band].snow_flux = 0.0;
            for (index = 0; index < options.Nnode - 1; index++) {
                energy[veg][band].T[index] = 0.0;
            }
            energy[veg][band].LongUnderOut = 0.0;
            energy[veg][band].Tfoliage = 0.0;
        }
    }

    // initialize Tfallback counters
    for (veg = 0; veg <= Nveg; veg++) {
        for (band = 0; band < options.SNOW_BAND; band++) {
            energy[veg][band].Tfoliage_fbcount = 0;
            energy[veg][band].Tcanopy_fbcount = 0;
            energy[veg][band].Tsurf_fbcount = 0;
            for (index = 0; index < options.Nnode - 1; index++) {
                energy[veg][band].T_fbcount[index] = 0;
            }
        }
    }
}
