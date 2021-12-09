/******************************************************************************
 * @section DESCRIPTION
 *
 * This routine initializes the soil variable arrays for each new grid cell.
 *****************************************************************************/

#include <vic_driver_shared_all.h>

/******************************************************************************
 * @brief    This routine initializes the soil variable arrays for each new
 *           grid cell.
 *****************************************************************************/
void
initialize_soil(cell_data_struct **cell,
                size_t             veg_num)
{
    extern option_struct options;

    size_t               veg, band, lindex, frost_area;

    for (veg = 0; veg <= veg_num; veg++) {
        for (band = 0; band < options.SNOW_BAND; band++) {
            // Prognostic states
            cell[veg][band].aero_resist[0] = 0.0;
            cell[veg][band].aero_resist[1] = 0.0;
            cell[veg][band].CLitter = 0.0;
            cell[veg][band].CInter = 0.0;
            cell[veg][band].CSlow = 0.0;
            for (lindex = 0; lindex < options.Nlayer; lindex++) {
                cell[veg][band].layer[lindex].Cs = 0.0;
                cell[veg][band].layer[lindex].T = 0.0;
                for (frost_area = 0; frost_area < options.Nfrost;
                     frost_area++) {
                    cell[veg][band].layer[lindex].ice[frost_area] = 0.0;
                }
                cell[veg][band].layer[lindex].kappa = 0.0;
                cell[veg][band].layer[lindex].moist = 0.0;
                cell[veg][band].layer[lindex].phi = 0.0;
            }
            cell[veg][band].rootmoist = 0.0;
            cell[veg][band].wetness = 0.0;
            // Fluxes
            cell[veg][band].pot_evap = 0.0;
            cell[veg][band].baseflow = 0.0;
            cell[veg][band].runoff = 0.0;
            cell[veg][band].inflow = 0.0;
            cell[veg][band].RhLitter = 0.0;
            cell[veg][band].RhLitter2Atm = 0.0;
            cell[veg][band].RhInter = 0.0;
            cell[veg][band].RhSlow = 0.0;
            cell[veg][band].RhTot = 0.0;
            for (lindex = 0; lindex < options.Nlayer; lindex++) {
                cell[veg][band].layer[lindex].esoil = 0.0;
                cell[veg][band].layer[lindex].transp = 0.0;
                cell[veg][band].layer[lindex].evap = 0.0;
            }
        }
    }
}
