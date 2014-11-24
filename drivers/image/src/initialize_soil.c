#include <vic_def.h>
#include <vic_run.h>
#include <vic_driver_image.h>

void
initialize_soil(cell_data_struct **cell,
                soil_con_struct   *soil_con,
                size_t             veg_num)

/**********************************************************************
        initialize_soil		Keith Cherkauer		July 31, 1996

   This routine initializes the soil variable arrays for each new
   grid cell.
**********************************************************************/
{
    extern option_struct options;

    size_t               veg, band, lindex, frost_area;
    double               tmp_moist[MAX_LAYERS];
    double               tmp_runoff;

    for (veg = 0; veg <= veg_num; veg++) {
        for (band = 0; band < options.SNOW_BAND; band++) {
            cell[veg][band].baseflow = 0;
            cell[veg][band].runoff = 0;
            for (lindex = 0; lindex < options.Nlayer; lindex++) {
                cell[veg][band].layer[lindex].evap = 0;
                cell[veg][band].layer[lindex].moist =
                    soil_con->init_moist[lindex];
                if (cell[veg][band].layer[lindex].moist >
                    soil_con->max_moist[lindex]) {
                    cell[veg][band].layer[lindex].moist =
                        soil_con->max_moist[lindex];
                }
                tmp_moist[lindex] = cell[veg][band].layer[lindex].moist;
                for (frost_area = 0; frost_area < options.Nfrost;
                     frost_area++) {
                    cell[veg][band].layer[lindex].ice[frost_area] = 0;
                }
            }
            compute_runoff_and_asat(soil_con, tmp_moist, 0,
                                    &(cell[veg][band].asat), &tmp_runoff);
            wrap_compute_zwt(soil_con, &(cell[veg][band]));
            cell[veg][band].CLitter = 0;
            cell[veg][band].CInter = 0;
            cell[veg][band].CSlow = 0;
        }
    }
    for (frost_area = 0; frost_area < options.Nfrost; frost_area++) {
        if (options.Nfrost == 1) {
            soil_con->frost_fract[frost_area] = 1.;
        }
        else if (options.Nfrost == 2) {
            soil_con->frost_fract[frost_area] = 0.5;
        }
        else {
            soil_con->frost_fract[frost_area] = 1. / (options.Nfrost - 1);
            if (frost_area == 0 || frost_area == options.Nfrost - 1) {
                soil_con->frost_fract[frost_area] /= 2.;
            }
        }
    }
}
