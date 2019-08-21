#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
#include <math.h>
#include <assert.h>

static char vcid[] =
    "$Id: gl_flow.c,v 5.12.2.20 2013/10/7 03:45:12 vicadmin Exp $";


int
gl_flow(snow_data_struct **snow,
        soil_con_struct   *soil_con,
        veg_con_struct    *veg_con,
        int                Nbands)
/**********************************************************************
        gl_flow      Bibi S. Naz	October 07, 2013

   this routine solve the glacier flow algorithm as part of VIC glacier flow 
   model. see more detail for ......


**********************************************************************/
{
    vicerror("gl_flow routine not complete yet, use GLACIER = SCALING in, your global parameter file");
    return(0);
}

int
gl_volume_area(snow_data_struct **snow,
               soil_con_struct   *soil_con,
               veg_con_struct    *veg_con,
               int                Nbands,
               int                dt)       // hours
/**********************************************************************
    gl_volume_area      Joe Hamman    February 26, 2014

    this routine solves the volume-area scaling glacier algorithm of
    Bahr, D. B., Global Distributions of Glacier Properties: A
    Stochastic Scaling Paradigm, Water Resour. Res., 33,
    1669-1679, 1997a.

**********************************************************************/
{
    extern option_struct options;

    int    iveg;
    int    band;
    int    bot_band;
    double cell_area;                       // km2
    double tile_area;                       // km2
    double ice_vol;                         // km3
    double ice_vol_new;                     // km3
    double max_vol;                         // km3
    double overflow_vol;                    // km3
    double ice_area_old;                    // km2
    double ice_area_temp;                   // km2
    double ice_area_new;                    // km2
    double cum_area;                        // km2
    double iwe_final;                       // mm
    double stepsize;                        // seconds

    cell_area = soil_con->cell_area / METERS_PER_KM / METERS_PER_KM; // m2 --> km2
    stepsize = dt * SECPHOUR;

    for (iveg = 0; iveg <= veg_con[0].vegetat_type_num; iveg++) {
        if (veg_con[iveg].Cv > 0.0) {
            // find total ice volume in veg tile
            tile_area = veg_con[iveg].Cv * cell_area;
            ice_vol = 0.0;
            ice_area_old = 0.0;
            for (band = 0; band < Nbands; band++) {
                snow[iveg][band].gl_overflow = 0.;
                if (snow[iveg][band].iwq > 0.0) {
                    ice_vol += snow[iveg][band].iwq * veg_con[iveg].Cv *
                               soil_con->AreaFract[band];
                    ice_area_old += soil_con->AreaFract[band] *
                                    veg_con[iveg].Cv;
                }
            }
            if (ice_vol > 0.0) {
                // ice_vol in km3 and ice_area in km2
                ice_vol *= cell_area / MMPERMETER / METERS_PER_KM;
                ice_area_old *= cell_area;

                // find new ice area
                ice_area_temp = pow((ice_vol / BAHR_C), (1 / BAHR_LAMBDA));

                // time scaling
                ice_area_new = ice_area_old +
                               (ice_area_temp - ice_area_old) * exp(stepsize / BAHR_T);

                // Check for negative new area after scaling due to ice volume loss
                // Ignore scaling in this edge case
                if (ice_area_new < 0.0) {
                    ice_area_new = ice_area_old;
                }

                // Make sure the area isn't too big
                overflow_vol = 0.;
                if (ice_area_new > tile_area) {
                    ice_area_new = tile_area;
                    if (options.GLACIER_OVERFLOW) {
                        max_vol =  BAHR_C * pow(tile_area, BAHR_LAMBDA);
                        overflow_vol = ice_vol - max_vol;
                        ice_vol = max_vol;
                    }
                }

                // determine where the ice belongs
                cum_area = 0.0;
                bot_band = -1;
                while ((cum_area < ice_area_new) && (bot_band < Nbands)) {
                    bot_band += 1;
                    cum_area += soil_con->AreaFract[bot_band] * cell_area *
                                veg_con[iveg].Cv;
                }

                // spread ice over bands
                iwe_final = ice_vol / (cum_area * METERS_PER_KM * METERS_PER_KM);
                for (band = 0; band < bot_band; band++) {
                    snow[iveg][band].iwq = iwe_final;
                }

                // add gl_runoff to snow structure (lowest band)
                if (options.GLACIER_OVERFLOW) {
                    snow[iveg][bot_band].gl_overflow = overflow_vol / 
                                                       (soil_con->AreaFract[bot_band] * 
                                                       veg_con[iveg].Cv * 
                                                       cell_area * 
                                                       METERS_PER_KM * 
                                                       METERS_PER_KM);
                }
                else {
                    snow[iveg][bot_band].gl_overflow = 0.;
                }
            }
        } /** end current vegetation type **/
    } /** end of vegetation loop **/
    return(0);
}
