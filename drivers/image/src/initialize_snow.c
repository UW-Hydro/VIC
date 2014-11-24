#include <vic_def.h>
#include <vic_run.h>
#include <vic_driver_image.h>

void
initialize_snow(snow_data_struct **snow,
                size_t             veg_num)

/**********************************************************************
        initialize_snow		Keith Cherkauer		January 22, 1997

   This routine initializes the snow variable arrays for each new
   grid cell.

   VARIABLES INITIALIZED:
    snow[i][j].snow;	          TRUE = snow, FALSE = no snow
    snow[i][j].last_snow;         time steps since last snowfall
    snow[i][j].snow_canopy;       amount of snow on canopy (m)
    snow[i][j].swq;               snow water equivalent of the entire pack (m)
    snow[i][j].surf_water;        liquid water content of the surface
                                  layer (m)
    snow[i][j].pack_water;        liquid water content of the snow pack (m)
    snow[i][j].surf_temp;         depth averaged temperature of the snow pack
                                  surface layer (C)
    snow[i][j].pack_temp;         depth averaged temperature of the snow pack
                                  (C)
    snow[i][j].vapor_flux;        depth of water evaporation, sublimation, or
                                  condensation from snow pack (m)
    snow[i][j].canopy_vapor_flux; depth of water evaporation, sublimation, or
                                  condensation from intercepted snow (m)
    snow[i][j].albedo;            snow surface albedo (fraction)
    snow[i][j].coldcontent;       cold content of snow pack
    snow[i][j].mass_error;        snow mass balance error
    snow[i][j].density;	          snow density (kg/m^3)
    snow[i][j].depth;	          snow depth (m)
    snow[i][j].tmp_int_storage;   temporary canopy storage, used in
                                  snow_canopy
    snow[i][j].Qnet;              Net energy error in snow model
    snow[i][j].band_elev;         median elevation of the current snow band
    snow[i][j].prec_frac;         fracton of precipitation that falls in the
                                                  current snow band

 *********************************************************************/
{
    extern option_struct options;
    size_t               i, j;

    for (i = 0; i <= veg_num; i++) {
        for (j = 0; j < options.SNOW_BAND; j++) {
            // State vars
            snow[i][j].albedo = 0.0;
            snow[i][j].canopy_albedo = 0.0;
            snow[i][j].coldcontent = 0.0;
            snow[i][j].coverage = 0.0;
            snow[i][j].density = 0.0;
            snow[i][j].depth = 0.0;
            snow[i][j].last_snow = 0;
            snow[i][j].max_snow_depth = 0.0;
            snow[i][j].MELTING = FALSE;
            snow[i][j].pack_temp = 0.0;
            snow[i][j].pack_water = 0.0;
            snow[i][j].snow = FALSE;
            snow[i][j].snow_canopy = 0.0;
            snow[i][j].store_coverage = 0.0;
            snow[i][j].store_snow = FALSE;
            snow[i][j].store_swq = 0.0;
            snow[i][j].surf_temp = 0.0;
            snow[i][j].surf_temp_fbflag = FALSE;
            snow[i][j].surf_temp_fbcount = 0;
            snow[i][j].surf_water = 0.0;
            snow[i][j].swq = 0.0;
            snow[i][j].snow_distrib_slope = 0.0;
            snow[i][j].tmp_int_storage = 0.0;
            // Fluxes
            snow[i][j].blowing_flux = 0.0;
            snow[i][j].canopy_vapor_flux = 0.0;
            snow[i][j].mass_error = 0.0;
            snow[i][j].melt = 0.0;
            snow[i][j].Qnet = 0.0;
            snow[i][j].surface_flux = 0.0;
            snow[i][j].transport = 0.0;
            snow[i][j].vapor_flux = 0.0;
        }
    }
}
