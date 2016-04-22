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

#include <vic_driver_shared_all.h>

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
    double               dummy_double;
    bool                 dummy_bool;
    size_t               dummy_size_t;
    unsigned int         dummy_unsigned_int;
    int                  dummy_int;

    // initialize miscellaneous energy balance terms
    for (veg = 0; veg <= Nveg; veg++) {
        for (band = 0; band < options.SNOW_BAND; band++) {
            // Prognostic states
            energy[veg][band].AlbedoLake = dummy_double;
            energy[veg][band].AlbedoOver = dummy_double;
            energy[veg][band].AlbedoUnder = dummy_double;
            energy[veg][band].Cs[0] = dummy_double;
            energy[veg][band].Cs[1] = dummy_double;
            energy[veg][band].frozen = dummy_bool;
            energy[veg][band].kappa[0] = dummy_double;
            energy[veg][band].kappa[1] = dummy_double;
            energy[veg][band].Nfrost = dummy_size_t;
            energy[veg][band].Nthaw = dummy_size_t;
            energy[veg][band].T1_index = dummy_int;
            energy[veg][band].Tcanopy = dummy_int;
            energy[veg][band].Tcanopy_fbflag = dummy_bool;
            energy[veg][band].Tcanopy_fbcount = dummy_unsigned_int;
            energy[veg][band].Tfoliage = 0.0; // dummy_double;
            energy[veg][band].Tfoliage_fbflag = dummy_bool;
            energy[veg][band].Tfoliage_fbcount = dummy_unsigned_int;
            energy[veg][band].Tsurf = dummy_double;
            energy[veg][band].Tsurf_fbflag = dummy_bool;
            energy[veg][band].Tsurf_fbcount = dummy_unsigned_int;
            energy[veg][band].unfrozen = dummy_double;
            for (index = 0; index < options.Nnode - 1; index++) {
                energy[veg][band].Cs_node[index] = dummy_double;
                energy[veg][band].ice[index] = dummy_double;
                energy[veg][band].kappa_node[index] = dummy_double;
                energy[veg][band].moist[index] = dummy_double;
                energy[veg][band].T[index] = dummy_double;
                energy[veg][band].T_fbflag[index] = dummy_bool;
                energy[veg][band].T_fbcount[index] = dummy_unsigned_int;
            }
            for (index = 0; index < MAX_FRONTS - 1; index++) {
                energy[veg][band].fdepth[index] = dummy_double;
                energy[veg][band].tdepth[index] = dummy_double;
            }
            // Fluxes
            energy[veg][band].advected_sensible = dummy_double;
            energy[veg][band].advection = dummy_double;
            energy[veg][band].AtmosError = dummy_double;
            energy[veg][band].AtmosLatent = dummy_double;
            energy[veg][band].AtmosLatentSub = dummy_double;
            energy[veg][band].AtmosSensible = dummy_double;
            energy[veg][band].canopy_advection = dummy_double;
            energy[veg][band].canopy_latent = dummy_double;
            energy[veg][band].canopy_latent_sub = dummy_double;
            energy[veg][band].canopy_refreeze = dummy_double;
            energy[veg][band].canopy_sensible = dummy_double;
            energy[veg][band].deltaCC = dummy_double;
            energy[veg][band].deltaH = dummy_double;
            energy[veg][band].error = dummy_double;
            energy[veg][band].fusion = dummy_double;
            energy[veg][band].grnd_flux = dummy_double;
            energy[veg][band].latent = dummy_double;
            energy[veg][band].latent_sub = dummy_double;
            energy[veg][band].longwave = dummy_double;
            energy[veg][band].LongOverIn = dummy_double;
            energy[veg][band].LongUnderIn = dummy_double;
            energy[veg][band].LongUnderOut = 0.0; //dummy_double;
            energy[veg][band].melt_energy = dummy_double;
            energy[veg][band].NetLongAtmos = dummy_double;
            energy[veg][band].NetLongOver = dummy_double;
            energy[veg][band].NetLongUnder = dummy_double;
            energy[veg][band].NetShortAtmos = dummy_double;
            energy[veg][band].NetShortGrnd = dummy_double;
            energy[veg][band].NetShortOver = dummy_double;
            energy[veg][band].NetShortUnder = dummy_double;
            energy[veg][band].out_long_canopy = dummy_double;
            energy[veg][band].out_long_surface = dummy_double;
            energy[veg][band].refreeze_energy = dummy_double;
            energy[veg][band].sensible = dummy_double;
            energy[veg][band].shortwave = dummy_double;
            energy[veg][band].ShortOverIn = dummy_double;
            energy[veg][band].ShortUnderIn = dummy_double;
            energy[veg][band].snow_flux = 0.0; //dummy_double;
        }
    }
}
