/******************************************************************************
 * @section DESCRIPTION
 *
 * clean up functions for routing extension
 *****************************************************************************/

#include <rout.h>

/******************************************************************************
 * @brief    Finalize RVIC by freeing memory.
 *****************************************************************************/
void
rout_finalize(void)
{
    extern rout_struct rout;

    free(rout.rout_param.source2outlet_ind);
    free(rout.rout_param.source_time_offset);
    free(rout.rout_param.source_x_ind);
    free(rout.rout_param.source_y_ind);
    free(rout.rout_param.source_lat);
    free(rout.rout_param.source_lon);
    free(rout.rout_param.source_VIC_index);
    free(rout.rout_param.outlet_lat);
    free(rout.rout_param.outlet_lon);
    free(rout.rout_param.outlet_VIC_index);
    free(rout.rout_param.unit_hydrograph);
    free(rout.rout_param.aggrunin);
    free(rout.discharge);
    free(rout.ring);
}
