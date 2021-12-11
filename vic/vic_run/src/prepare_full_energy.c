/******************************************************************************
 * @section DESCRIPTION
 *
 * This subroutine returns the soil thermal properties, moisture and ice
 * contents for the top two layers for use with the QUICK_FLUX ground heat flux
 * solution.
 *****************************************************************************/

#include <vic_run.h>

/******************************************************************************
 * @brief    This subroutine returns the soil thermal properties, moisture and
 *           ice contents for the top two layers for use with the QUICK_FLUX
 *           ground heat flux solution.
 *****************************************************************************/
void
prepare_full_energy(cell_data_struct  *cell,
                    energy_bal_struct *energy,
                    soil_con_struct   *soil_con,
                    double            *moist0,
                    double            *ice0)
{
    extern option_struct options;

    size_t               i;
    layer_data_struct   *layer = NULL;

    layer = calloc(options.Nlayer, sizeof(*layer));
    check_alloc_status(layer, "Memory allocation error.");

    for (i = 0; i < options.Nlayer; i++) {
        layer[i] = cell->layer[i];
    }

    /* Compute top soil layer moisture content (mm/mm) */
    *moist0 = layer[0].moist / (soil_con->depth[0] * MM_PER_M);

    /* Compute top soil layer ice content (mm/mm) */
    if (options.FROZEN_SOIL && soil_con->FS_ACTIVE) {
        if ((energy->T[0] + energy->T[1]) / 2. < 0.) {
            *ice0 = *moist0 -
                    maximum_unfrozen_water((energy->T[0] + energy->T[1]) / 2.,
                                           soil_con->max_moist[0] /
                                           (soil_con->depth[0] * MM_PER_M),
                                           soil_con->bubble[0],
                                           soil_con->expt[0]);
            if (*ice0 < 0.) {
                *ice0 = 0.;
            }
        }
        else {
            *ice0 = 0.;
        }
    }
    else {
        *ice0 = 0.;
    }

    /** Compute Soil Thermal Properties **/
    compute_soil_layer_thermal_properties(layer, soil_con->depth,
                                          soil_con->bulk_dens_min,
                                          soil_con->soil_dens_min,
                                          soil_con->quartz,
                                          soil_con->bulk_density,
                                          soil_con->soil_density,
                                          soil_con->organic,
                                          soil_con->frost_fract,
                                          options.Nlayer);

    /** Save Thermal Conductivities for Energy Balance **/
    energy->kappa[0] = layer[0].kappa;
    energy->Cs[0] = layer[0].Cs;
    energy->kappa[1] = layer[1].kappa;
    energy->Cs[1] = layer[1].Cs;

    free((char *) layer);
}
