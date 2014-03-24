#include <stdio.h>
#include <stdlib.h>
#include <vicNl.h>
#include <math.h>

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
    double kterminus;
    double Hhi;
    double gradhi;
    double areahi;
    double Hlo;
    double gradlo;
    double arealo;
    double Qnet;
    double delH;
    double divQ;
    double gamma = 1.0e-8; // tune this value for the ice dynamics
    double delswe;
    double deliwe;
    int    Lvert = 100; // m - length scale for surface slopes/ice flux
    int    iveg;
    int    band;
    int    Glbands;

    for (iveg = 0; iveg <= veg_con[0].vegetat_type_num; iveg++) {
        if (veg_con[iveg].Cv > 0.0) {
            // find the highest ice
            for (band = 0; band < Nbands; band++) {
                printf("gl_area=%f\n", snow[iveg][band].iwq);
                if ((soil_con->GlAreaFract[band] > 0) &&
                    (soil_con->GlAreaFract[band + 1] == 0)) {
                    Glbands = band;
                    printf("glbands=%d\n", Glbands);
                }
                else {
                    Glbands = Nbands;
                }
            }
            // look for the first occurence and flag the terminus
            for (band = 0; band < Glbands; band++) {
                if (soil_con->AreaFract[band] > 0) {
                    printf("band= %d\n", band);
                    // change in swe in previous month
                    delswe = snow[iveg][band].swq - snow[iveg][band].swqold;
                    snow[iveg][band].swqold = snow[iveg][band].swq;
                    // change in iwe in previous month
                    deliwe = snow[iveg][band].iwq - snow[iveg][band].iwqold;
                    snow[iveg][band].iwqold = snow[iveg][band].iwq;
                    snow[iveg][band].bn = delswe + deliwe; // glacier mass balance

                    if (band > 0) {
                        if ((snow[iveg][band].iwq == 0) &&
                            (snow[iveg][band - 1].iwq > 0)) {
                            kterminus = band - 1;
                        }
                    }

                    // calculate the flux into the cell
                    if (band < Glbands) { // below the top cell
                        if (snow[iveg][band + 1].iwq > 0) { // ice above
                            Hhi =
                                (snow[iveg][band].iwq +
                                 snow[iveg][band + 1].iwq) / 2.;
                            gradhi =
                                (soil_con->BandElev[band +
                                                    1] -
                                 soil_con->BandElev[band]) / Lvert;
                            areahi =
                                (soil_con->GlAreaFract[band] +
                                 soil_con->GlAreaFract[band + 1]) / 2;
                        }
                        else {
                            Hhi = 0;
                            gradhi = 0;
                            areahi = 0;
                        }
                    }
                    else { // top cell
                        if (snow[iveg][band].iwq > 0) { // ice here
                            Hhi = snow[iveg][band].iwq;
                            gradhi = 0; // shuts off incoming ice
                            areahi = 0;
                        }
                        else {
                            Hhi = 0;
                            gradhi = 0;
                            areahi = 0;
                        }
                    }
                    // calculate the flux out of the cell
                    // above the valley bottom
                    if (band > 0) {
                        gradlo =
                            (soil_con->BandElev[band] -
                             soil_con->BandElev[band - 1]) / Lvert;
                        Hlo =
                            (snow[iveg][band].iwq +
                             snow[iveg][band - 1].iwq) / 2;
                        arealo =
                            (soil_con->GlAreaFract[band] +
                             soil_con->GlAreaFract[band - 1]) / 2;
                    }
                    else {
                        Hlo = snow[iveg][band].iwq;
                        gradlo = 1;
                        arealo = soil_con->GlAreaFract[band];
                    }

                    // assign fluxes at the interface
                    snow[iveg][band].Qin = gamma * areahi * pow(Hhi, 5) * pow(
                        gradhi, 3);
                    snow[iveg][band].Qout = gamma * arealo * pow(Hlo, 5) * pow(
                        gradlo, 3);
                    Qnet = snow[iveg][band].Qin - snow[iveg][band].Qout;
                    if ((Qnet > 0) && (soil_con->GlAreaFract[band] == 0)) {
                        // new incoming ice
                        soil_con->GlAreaFract[band] =
                            soil_con->GlAreaFract[band + 1];
                    }

                    if (abs(Qnet) > 0) {
                        if ((Qnet / soil_con->GlAreaFract[band]) < 5) {
                            divQ = Qnet / soil_con->GlAreaFract[band];
                        }
                        else {
                            divQ = 5;
                        }
                    }
                    else {
                        divQ = 0;
                    }

                    delH = snow[iveg][band].bn + divQ;
                    if ((snow[iveg][band].iwq + delH) > 0) {
                        snow[iveg][band].iwq = snow[iveg][band].iwq + delH;
                    }
                    else {
                        snow[iveg][band].iwq = 0;
                    }

                    soil_con->BandElev[band] = soil_con->BandElev[band] + delH;
                    // don't melt into the rock
                    if (snow[iveg][band].iwq < 0) {
                        snow[iveg][band].iwq = 0;
                        soil_con->GlAreaFract[band] = 0;
                    }
                }
            } // End loop through elevation bands
        } /** end current vegetation type **/
    } /** end of vegetation loop **/
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
    int    iveg;
    int    band;
    int    bot_band;
    double cell_area;                       // km2
    double tile_area;                       // km2
    double ice_vol;                         // km3
    double ice_area_old;                    // km2
    double ice_area_temp;                   // km2
    double ice_area_new;                    // km2
    double cum_area;                        // km2
    double iwe_final;                       // mm
    double stepsize;                        // seconds

    cell_area = soil_con->cell_area / MPERKILOMETER / MPERKILOMETER; // m2 --> km2
    stepsize = dt * SECPHOUR;

    for (iveg = 0; iveg <= veg_con[0].vegetat_type_num; iveg++) {
        if (veg_con[iveg].Cv > 0.0) {
            // find total ice volume in veg tile
            tile_area = veg_con[iveg].Cv * cell_area;
            ice_vol = 0.0;
            ice_area_old = 0.0;
            for (band = 0; band < Nbands; band++) {
                if (snow[iveg][band].iwq > 0.0) {
                    ice_vol += snow[iveg][band].iwq * veg_con[iveg].Cv *
                               soil_con->AreaFract[band];
                    ice_area_old += soil_con->AreaFract[band] *
                                    veg_con[iveg].Cv;
                }
            }

            if (ice_vol > 0) {
                // ice_vol in km3 and ice_area in km2
                ice_vol *= cell_area / MMPERMETER / MPERKILOMETER;
                ice_area_old *= veg_con[iveg].Cv * cell_area;

                // find new ice area
                ice_area_temp = pow((ice_vol / BAHR_C), (1 / BAHR_LAMBDA));

                // time scaling
                ice_area_new = ice_area_temp +
                               (ice_area_old - ice_area_temp) * exp(
                    stepsize / BAHR_T);

                // Make sure the area isn't too big
                if (ice_area_new > tile_area) {
                    ice_area_new = tile_area;
                }

                // determine where the ice belongs
                cum_area = 0.0;
                bot_band = -1;
                while (cum_area < ice_area_new && bot_band < Nbands) {
                    bot_band += 1;
                    cum_area += soil_con->AreaFract[bot_band] * cell_area *
                                veg_con[iveg].Cv;
                }

                // spread ice over bands
                iwe_final = ice_vol / cum_area * MPERKILOMETER * MPERKILOMETER;
                for (band = 0; band < bot_band; band++) {
                    snow[iveg][band].iwq = iwe_final;
                }
            }
        } /** end current vegetation type **/
    } /** end of vegetation loop **/
    return(0);
}
