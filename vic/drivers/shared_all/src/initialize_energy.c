/******************************************************************************
 * @section DESCRIPTION
 *
 * Initialize energy structure.
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

    // initialize miscellaneous energy balance terms
    for (veg = 0; veg <= Nveg; veg++) {
        for (band = 0; band < options.SNOW_BAND; band++) {
            // Prognostic states
            energy[veg][band].AlbedoLake = 0.0;
            energy[veg][band].AlbedoOver = 0.0;
            energy[veg][band].AlbedoUnder = 0.0;
            energy[veg][band].Cs[0] = 0.0;
            energy[veg][band].Cs[1] = 0.0;
            energy[veg][band].frozen = false;
            energy[veg][band].kappa[0] = 0.0;
            energy[veg][band].kappa[1] = 0.0;
            energy[veg][band].Nfrost = 0;
            energy[veg][band].Nthaw = 0;
            energy[veg][band].T1_index = 0;
            energy[veg][band].Tcanopy = 0;
            energy[veg][band].Tcanopy_fbflag = false;
            energy[veg][band].Tcanopy_fbcount = 0;
            energy[veg][band].Tfoliage = 0.0;
            energy[veg][band].Tfoliage_fbflag = false;
            energy[veg][band].Tfoliage_fbcount = 0;
            energy[veg][band].Tsurf = 0.0;
            energy[veg][band].Tsurf_fbflag = false;
            energy[veg][band].Tsurf_fbcount = 0;
            energy[veg][band].unfrozen = 0.0;
            for (index = 0; index < options.Nnode - 1; index++) {
                energy[veg][band].Cs_node[index] = 0.0;
                energy[veg][band].ice[index] = 0.0;
                energy[veg][band].kappa_node[index] = 0.0;
                energy[veg][band].moist[index] = 0.0;
                energy[veg][band].T[index] = 0.0;
                energy[veg][band].T_fbflag[index] = false;
                energy[veg][band].T_fbcount[index] = 0;
            }
            for (index = 0; index < MAX_FRONTS - 1; index++) {
                energy[veg][band].fdepth[index] = 0.0;
                energy[veg][band].tdepth[index] = 0.0;
            }
            // Fluxes
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
        }
    }
}
