/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine resets the values of all output variables to 0.
 *****************************************************************************/

#include <vic_driver_shared_all.h>

/******************************************************************************
 * @brief    This routine resets the values of all output variables to 0.
 *****************************************************************************/
void
zero_output_list(double **out_data)
{
    extern metadata_struct out_metadata[N_OUTVAR_TYPES];

    size_t                 varid, i;

    for (varid = 0; varid < N_OUTVAR_TYPES; varid++) {
        for (i = 0; i < out_metadata[varid].nelem; i++) {
            out_data[varid][i] = 0.;
        }
    }
}
